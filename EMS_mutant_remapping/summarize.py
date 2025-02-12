import argparse
import os

def list_vcf_files(directory):
    """List all VCF files in the given directory."""
    #ignore indel files
    #only get the substitution files
    return [file for file in os.listdir(directory) if file.endswith('maps.vcf')]

def extract_parameter(vcf_filename):
    """Extract parameter from the VCF filename."""
    try:
        return vcf_filename.split('HetMinPer15')[1].split('.reformatted')[0]
    except IndexError:
        raise ValueError(f"Unexpected filename format: {vcf_filename}")

def process_vcf_files(vcfs, params):
    """Process all listed VCF files to calculate mutations."""
    mutations = {}

    for index, vcf in enumerate(vcfs):
        print(f'Processing {vcf}...')
        param = extract_parameter(vcf)

        #make sure the parameter is within the predefined parameters
        if param not in params:
            print(f"Warning: {param} not found in expected params. Skipping...")
            continue

        print(f'Stringency level: {params.index(param) + 1}')
        calculate_mutations(vcf, index, mutations)

    return mutations

def calculate_mutations(vcf, index, mutations):
    """Calculate mutations for a given VCF file."""
    try:
        with open(vcf, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue

                items = line.strip().split('\t')

                if line.startswith('#'):
                    continue

                original, mutated = items[3].upper(), items[4].upper()
                accession = extract_accession(items[-1])

                update_mutation_counts(mutations, accession, original, mutated, index)

    except FileNotFoundError:
        print(f"Error: File {vcf} not found.")

def extract_accession(info_field):
    """Extract the accession from the INFO field safely."""
    for field in info_field.split(';'):
        if field.startswith('seed_avail='):
            return field.split('=')[1]
    return 'UNKNOWN'

def update_mutation_counts(mutations, accession, original, mutated, index):
    """Update mutation counts based on original and mutated values."""
    if accession not in mutations:
        #fields will be formatted as below
        #['# Substitutions', '# EMS-type mutations', '# non-EMS-type mutations', 'G->A', 'C->T', '%EMS']
        mutations[accession] = [0] * (len(vcfs) * 6)

    if original != mutated:
        # Substitution counts
        mutation_types = {'GA': 3, 'CT': 4}
        transition = original + mutated

        if transition in mutation_types:  # EMS-type mutation
            mutations[accession][index * 6 + mutation_types[transition]] += 1
        else:
            mutations[accession][index * 6 + 2] += 1  # Non EMS-type mutation

def calculate_final_params(mutations, vcfs, params, sorted_params):
    """Determine final parameters based on EMS mutation rate thresholds."""
    Fparams = {}

    for key in mutations.keys():
        per_ems = [mutations[key][i * 6 + 5] for i in range(len(vcfs))]  # Extract EMS mutation rates
        param_to_ems = dict(zip(params, per_ems))

        high, medium, low = None, None, None

        #if  > 98% EMS exists, least stringent parameter with this rate
        #is high-confidence
        for param in sorted_params:
            if param_to_ems.get(param, 0) >= 0.98:
                high = param
                break

        if high:
            remaining_params = [x for x in sorted_params if x != high]
        else:
            #if >98% EMS does not exist, predefined parameters are used for high, medium, low
            Fparams[key] = ['HetMinCov5HomMinCov3', 'HetMinCov4HomMinCov3', 'HetMinCov3HomMinCov2']
            continue

        for param in remaining_params:
            #medium: least stringent with >97% EMS
            if param_to_ems.get(param, 0) >= 0.97:
                medium = param
                break

        remaining_params = [x for x in remaining_params if x != medium]

        for param in remaining_params:
            #low: least stringent parameters with >95% EMS
            if param_to_ems.get(param, 0) >= 0.95:
                low = param
                break

        Fparams[key] = [high] + ([medium] if medium else ['N/A']) + ([low] if low else ['N/A'])

    return Fparams

def summarize_mutations(mutations, vcfs, params):
    """Prepare the mutation summary and ensure correct sorting."""
    numeric_keys = [int(k.replace('Kronos', '')) for k in mutations if k.replace('Kronos', '').isdigit()]
    non_numeric_keys = [k for k in mutations if not k.replace('Kronos', '').isdigit()]

    sorted_keys = ['Kronos' + str(k) for k in sorted(numeric_keys)] + sorted(non_numeric_keys)

    for key in sorted_keys:
        for i in range(len(vcfs)):
            base_idx = i * 6
            mutations[key][base_idx] = sum(mutations[key][base_idx + 2: base_idx + 5])
            mutations[key][base_idx + 1] = sum(mutations[key][base_idx + 3: base_idx + 5])
            mutations[key][base_idx + 5] = round(
                mutations[key][base_idx + 1] / mutations[key][base_idx], 4
            ) if mutations[key][base_idx] != 0 else 0

    return sorted_keys

def write_summary_file(mutations, vcfs, params, sorted_keys, Fparams):
    """Write the mutation summary to a file."""
    headers = [''] + [f"{p}\t \t \t \t \t" for p in params]
    descriptions = ['Line'] + ['# Substitutions', '# EMS-type mutations', '# non-EMS-type mutations', 'G->A', 'C->T', '%EMS'] * len(vcfs) + ['High threshold', 'Medium threshold', 'Low threshold']

    try:
        with open('Mutations.summary', 'w') as file:
            file.write('\t'.join(headers) + '\n')
            file.write('\t'.join(descriptions) + '\n')

            for key in sorted_keys:
                mutation_details = '\t'.join(map(str, mutations[key]))
                final_param = '\t'.join(Fparams[key])
                file.write(f'{key}\t{mutation_details}\t{final_param}\n')

        print("Summary file written successfully.")

    except IOError as e:
        print(f"Error writing summary file: {e}")

# Main execution block
if __name__ == '__main__':
    vcfs = list_vcf_files(os.getcwd())

    if not vcfs:
        print("No VCF files found. Exiting...")
        exit(1)

    params = [extract_parameter(v) for v in vcfs]
    #from least stringency to most stringency
    sorted_params = ['HetMinCov3HomMinCov2', 'HetMinCov4HomMinCov3', 'HetMinCov5HomMinCov3', 'HetMinCov6HomMinCov4']

    mutations = process_vcf_files(vcfs, params)
    sorted_keys = summarize_mutations(mutations, vcfs, params)
    Fparams = calculate_final_params(mutations, vcfs, params, sorted_params)

    write_summary_file(mutations, vcfs, params, sorted_keys, Fparams)
