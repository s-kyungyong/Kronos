import pandas as pd
import re
from Bio import SeqIO
from Bio.Blast import NCBIXML
from intervaltree import Interval, IntervalTree
import os

#3-mer start sites from Kronos
#and some manually checked NLRs
acceptable_starts = set('''
MAA MAD MAE MAF MAG MAH MAI MAK MAL MAM MAN MAP MAQ MAR MAS MAT MAV MAY
MCA MCP MDA MDE MDF MDG MDH MDI MDL MDM MDP MDQ MDR MDS MDT MDV MDW MDY
MEA MED MEE MEF MEG MEH MEI MEK MEL MEM MEN MEP MEQ MER MES MET MEV MEY
MFL MFR MFS MFV MGA MGD MGE MGF MGG MGM MGP MGR MGS MGT MHD MHE MHR MIA
MIE MIG MIK MKG MKP MLC MLE MLG MLI MLN MLQ MLR MLS MLV MMD MME MMK MMQ
MMS MMW MNA MND MNG MNL MNM MNS MNY MPA MPD MPI MPK MPP MPQ MPR MPW MQA
MQF MQG MQT MQV MRE MRG MRI MRK MRP MRS MRV MRW MRY MSA MSD MSE MSG MSI
MSK MSL MSM MSN MSP MSQ MSR MSS MST MSV MTA MTD MTE MTG MTH MTI MTL MTS
MTT MTV MVA MVD MVE MVG MVK MVL MVM MVP MVR MVS MVT MVV MWS MYA MYD MYT
MYV MYY MPV
'''.split())

def parse_gff3(gff_file):
    gene_cds = {}
    gene_loci = {}
    gene_direction = {}

    with open(gff_file) as gff:
        for line in gff:
            if not 'CDS' in line:
                continue

            fields = line.strip().split('\t')

            if fields[2] != 'CDS':
                continue

            chrom, _, _, start, end, _, direction, _, attributes = fields
            gene_id = re.search('Parent=([^;]+)', attributes).group(1)
            gene_cds.setdefault(gene_id, []).append((int(start), int(end)))
            gene_direction[gene_id] = direction
            gene_loci[gene_id] = (chrom, min(int(start), int(end)), max(int(start), int(end)))

    return gene_cds, gene_loci, gene_direction

def parse_blastp_xml(xml_file):
    results = []
    with open(xml_file) as handle:
        for record in NCBIXML.parse(handle):
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    results.append((
                        record.query, alignment.hit_def, float(hsp.identities), int(hsp.align_length), int(hsp.gaps),
                        hsp.query, hsp.sbjct, int(record.query_length),
                        int(hsp.query_start), int(hsp.query_end), int(alignment.length),
                        int(hsp.sbjct_start), int(hsp.sbjct_end)
                    ))
    return results

def calculate_splice_sites(cds_list, direction):
    cds_sorted = sorted(cds_list, key=lambda x: x[0], reverse=(direction == '-'))
    cumulative = 0
    splice_sites = []
    for i in range(len(cds_sorted) - 1):
        start, end = cds_sorted[i]
        cumulative += abs(end - start) + 1
        splice_sites.append(cumulative // 3)
    return splice_sites

def recount_position(index, aligned_seq):
    count = 0
    for idx, aa in enumerate(aligned_seq):
        if aa != "-": count += 1
        if count == index: return idx
    return None

def map_splice_sites_to_alignment(splice_sites, aligned_seq, query_from):
    mapped = [recount_position(s - query_from + 1, aligned_seq) for s in splice_sites]
    global_gap = aligned_seq.count('-') / len(aligned_seq)
    for s in mapped:
        if s is None: return False
        window = aligned_seq[max(0, s - 1): min(len(aligned_seq), s + 14)]
        local_gap = window.count('-') / len(window)
        if local_gap > global_gap * 1.2:
            return False

    return True

def geneID(transcript, prefix, outtype ):

    if prefix in ['maker_v1', 'maker_v2']:
        gene = transcript.rsplit('-', 1)[0]
        num  = transcript.rsplit('-', 1)[1]

    elif prefix == 'augustus':
        gene = transcript.rsplit('.', 1)[0]
        num  = transcript.rsplit('.', 1)[1].replace('t', '')

    elif prefix == 'snap':
        gene = transcript.rsplit('.', 1)[0]
        num  = transcript.rsplit('.', 1)[1].replace('t', '')

    return gene if outtype == 'gene' else num


def examine_alignment(result, pred_cds, gene_cds, direction, prefix, db):
    qid, sid, identity, aln_len, gaps, qseq, sseq, qlen, q_from, q_to, slen, s_from, s_to = result
    strand = direction[qid]
    query_cds = pred_cds[qid].seq
    splice_sites = calculate_splice_sites(gene_cds[qid], strand) if qid in gene_cds else []

    qcov = (q_to - q_from + 1) / qlen
    scov = (s_to - s_from + 1) / slen
    perid = identity / (q_to - q_from + 1)

    status = {'gene': geneID(qid, prefix, 'gene'), 'transcript': qid,
              'length': qlen, 'reliable': False, 'reason': '', 'selected': False,
              'prefix': prefix, 'DB': db, 'score': 0}

    query_seq = query_cds.translate()
    #sometimes genes may have the correct starting positions
    #but start with some random sequences
    #just check for a stop codon
    if not query_seq.endswith('*'):
        status['reason'] = 'Incomplete'
        return status

    #no glboal matches. discard.
    #if qcov < 0.9 or scov < 0.9:
    #    status['reason'] = 'Partial match'
    #    return status

    #aligned length
    #equals to query and subject length
    if aln_len == qlen == slen:
        status['reliable'] = True
        status['reason'] = 'Global match'
        status['score'] = 2
        return status

    #if sequence identity is low (i.e. there are not good close homolog matches)
    #if the gene is a single exon, and > 750 amino acids, maybe it is an okay divergent homolog
    #just confirm the first some amino acids are aligned to the beginning of the homolog for the presence of CC
    if perid < 0.70 and splice_sites == [] and qlen > 750 and s_from < 150:
        status['reliable'] = True
        status['reason'] = 'Divergent single-CDS'
        status['score'] = 1
        return status

    #regions around splicing sites are gappy
    #likeley mutations around these regions
    #preventing proper translation
    if splice_sites and not map_splice_sites_to_alignment(splice_sites, qseq, q_from):
        status['reason'] = 'Bad splice junction'
        return status

    #examine start sties
    #query and subject start at position 1
    if q_from == s_from == 1:
        status['score'] += 1

    #query and subject start in different positions
    #but query has an accetable start sequences
    #in this case, make sure the hit start is contained within the query
    elif str(query_seq[:3]) in acceptable_starts and s_from == 1:
        status['score'] += 0.75

    #start annotation incorrect but the canonical start is contained in the gene
    elif q_from != 1 and s_from == 1 and str(query_seq[q_from - 1]) == "M":
        status['reason'] += 'Canonical start site contained '
        status['score'] += 0.75

    #allow some mismatches based on sequence divergence
    elif q_from < 4 + min(1 - perid, 0.3) * 30 and s_from < 4+ min(1 - perid, 0.3) * 30:
        status['reason'] += 'Start slightly off '
        status['score'] += 0.75

    #discard other cases
    else:
        status['reason'] = 'Bad start'
        return status

    #examine stop positions
    #query and subject terminate at the end
    if qlen == q_to and slen == s_to:
        status['score'] += 1

    #query sequences were aligned to the end
    #but 50 or less hit sequences were not covered
    #i.e. query ended slightly early
    elif (qlen == q_to and abs(slen - s_to) <= 50):
        status['reason'] += 'Early termination '
        status['score'] += 0.75

    #some termination site may be tricky due to small termination exons.
    #allow some unalignment at the end
    elif (abs(qlen - q_to) < 15 + min(1 - perid, 0.3) * 40 and abs(slen - s_to) < 15 + min(1 - perid, 0.3) * 40):
        status['reason'] += 'End slightly off'
        status['score'] += 0.75

    #otherwise discard
    else:
        status['reason'] = 'Bad stop'
        return status

    status['reliable'] = True
    if status['reason'] == '':
        status['reason'] = 'Good alignment'

    return status

def pick_model(records, prefix):
    # Filter for reliable records
    reliable = [r for r in records if r["reliable"]]
    if not reliable:
        return None

    def sort_key(record):
        # Higher score first, longer transcript, then smaller transcript number
        score = -record.get("score", 0)
        length = -record["length"]
        transcript_num = int(geneID(record["transcript"], prefix, 'number'))
        return (score, length, transcript_num)

    return sorted(reliable, key=sort_key)[0]

def filter_records(records, prefix):
    uniq_ids = list(set(record["gene"] for record in records if record["reliable"]))
    best_models = []
    for uid in uniq_ids:
        transcripts = [record for record in records if record["gene"] == uid]
        best = pick_model(transcripts, prefix)
        if best:
            best_models.append(best)
    return best_models

# Inputs
prefixes = ['maker_v2', 'maker_v1', 'augustus', 'snap']  # from high to low confidence

records_final = []
selected_by_transcript = {}
selected_coordinates = {}
saved_pep = {}
folder = os.getcwd().split('/')[-1]

# Load data
for prefix in prefixes:
    records = []
    gff = f'{prefix}.gff3'
    cds = f'{prefix}.cds.fa'
    dmnd = f'{prefix}.est.dmnd.out'

    gene_cds, gene_loci, gene_direction = parse_gff3(gff)
    pep_dict = SeqIO.to_dict(SeqIO.parse(cds, 'fasta'))
    results = parse_blastp_xml(dmnd)

    # for EST comparison
    for r in results:
        status = examine_alignment(r, pep_dict, gene_cds, gene_direction, prefix, 'long-read transcripts')
        records.append(status)

    # for Kronos comparison
    dmnd = f'{prefix}.kronos.dmnd.out'
    results = parse_blastp_xml(dmnd)
    for r in results:
        status = examine_alignment(r, pep_dict, gene_cds, gene_direction, prefix, 'Kronos reliable')
        records.append(status)

    best_models = filter_records(records, prefix)

    for best_model in best_models:
        transcript = best_model['transcript']
        chrom, start, end = gene_loci[transcript]
        score = best_model["score"]
        better_model = True
        key = f'{prefix}:{transcript}'
        saved_pep[key] = pep_dict[transcript]
        overlapping_loci = list()

        chrom_tree = IntervalTree()
        for item in selected_coordinates:
            c, p, t = item.split(':')

            if c == chrom:
                for start, end in selected_coordinates[item]:
                    if end > start:
                        chrom_tree.addi(start, end, f'{p}:{t}')

        for exon_start, exon_end in gene_cds[transcript]:
            overlaps = chrom_tree.overlap(exon_start, exon_end)
            for interval in overlaps:
                overlapping_loci.append((interval.begin, interval.end, interval.data))

        overlapping_transcripts = set(t for s, e, t in overlapping_loci)

        if len(overlapping_transcripts) == 1:
            prev_transcript = list(overlapping_transcripts)[0]
            prev_model = selected_by_transcript.get(prev_transcript)

            if prev_model and score > prev_model['score']:
                del selected_by_transcript[prev_transcript]
                del selected_coordinates[f'{chrom}:{prev_transcript}']
                selected_by_transcript[key] = best_model
                selected_coordinates[f'{chrom}:{key}'] = gene_cds[transcript]

            else:
                better_model = False

        elif len(overlapping_transcripts) > 1:
            better_model = False

        if better_model:
            selected_by_transcript[key] = best_model
            selected_coordinates[f'{chrom}:{key}'] = gene_cds[transcript]

    for p in prefixes:
        ct = 0

        for key in selected_by_transcript:
            if key.startswith(p):
                ct += 1
        print(f'{p}: {ct}')

    records_final += records

# Finalize selection status
selected_ids = set(selected_by_transcript.keys())
for record in records_final:
    record['selected'] = False
    if f'{record["prefix"]}:{record["transcript"]}' in selected_ids:
        if record['DB'] == selected_by_transcript[f'{record["prefix"]}:{record["transcript"]}']['DB']:
            record['selected'] = True

# Write selected protein sequences to FASTA
with open(folder + ".filtered.putative_nlrs.fa", 'w') as o:
    for i, transcript_id in enumerate(selected_ids, start=1):
        prefix, transcript = transcript_id.split(':')
        seq = saved_pep[transcript_id].seq.translate().replace("*", "")
        o.write(f'>{folder}_{prefix}_{i}    source:{transcript_id}\n{seq}\n')

# Save all records (selected and unselected)
pd.DataFrame(records_final).to_csv('NLR_evaluation.tsv', sep='\t', index=False)
