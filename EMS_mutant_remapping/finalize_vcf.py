import os
import sys

def read_parameters(filename):
    """Read parameter assignments from a file."""
    predefined = ['HetMinCov3HomMinCov2', 'HetMinCov4HomMinCov3', 'HetMinCov5HomMinCov3', 'HetMinCov6HomMinCov4', 'N/A']
    params = {}

    try:
        with open(filename, 'r') as file:
            for i, line in enumerate(file):
                if i == 0:
                    continue

                fields = line.strip().split()
                accession = fields[0]
                param = fields[-3:]  # High, Medium, Low

                if any(p not in predefined for p in param):
                    print(f"Warning: {accession} contains an undefined parameter {param}")
                    continue

                params[accession] = param  # Store confidence levels

    except FileNotFoundError:
        print(f"Error: File not found - {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        sys.exit(1)

    return params

def extract_parameter(vcf_filename):
    """Extract parameter from the VCF filename."""
    try:
        return vcf_filename.split('HetMinPer15')[1].split('.reformatted')[0]
    except IndexError:
        raise ValueError(f"Unexpected filename format: {vcf_filename}")

def get_confidence(params, accession, param):
    """Determine the confidence level of a parameter for a given accession."""
    confidence_levels = ['High', 'Medium', 'Low']
    try:
        return confidence_levels[params[accession].index(param)]
    except (ValueError, KeyError):
        return None  # Accession not in parameter file or param not found

def process_vcf(params, vcfs):
    """Process VCF files and store high, medium, and low-confidence calls."""
    variant_dict = {}  # Store calls with highest confidence only

    for vcf, label, desc in vcfs:
        param = extract_parameter(vcf)

        try:
            with open(vcf, 'r') as file:
                for line in file:
                    if line.startswith('#'):
                        continue  # Skip headers

                    fields = line.strip().split('\t')
                    chrom, pos, _, ref, alt, _, _, info_field = fields[:8]

                    try:
                        accession = next((x.split('=')[1] for x in info_field.split(';') if x.startswith('seed_avail=')), None)
                        if accession is None:
                            print(f"Warning: No accession found in {vcf} at {chrom}:{pos}")
                            continue
                        position_key = (chrom, int(pos), accession)
                    except IndexError:
                        print(f"Warning: Could not extract accession from {vcf} at {chrom}:{pos}")
                        continue

                    confidence = get_confidence(params, accession, param)
                    if not confidence:
                        continue

                    # Ensure highest confidence call is retained
                    # Prioritize uniquely mapped reads over multi-mapped reads for each position
                    if position_key not in variant_dict:
                        variant_dict[position_key] = (line, confidence, param, label, desc)

                    else:
                        existing_confidence = variant_dict[position_key][1]
                        existing_label = variant_dict[position_key][3]
                        confidence_priority = {'High': 1, 'Medium': 2, 'Low': 3}

                        if confidence_priority[confidence] == confidence_priority[existing_confidence]:

                            if existing_label == 'Multi' and label == 'Unique':
                                variant_dict[position_key] = (line, confidence, param, label, desc)
                            else:
                                pass

                        elif confidence_priority[confidence] < confidence_priority[existing_confidence]:
                            variant_dict[position_key] = (line, confidence, param, label, desc)

        except FileNotFoundError:
            print(f"Error: File {vcf} not found.")
        except Exception as e:
            print(f"Error processing {vcf}: {e}")

    return variant_dict

def write_final_vcf(output_filename, variant_dict):
    """Write the sorted high-confidence, medium, and low-confidence variants to a VCF file."""
    sorted_variants = sorted(variant_dict.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))  # Sort by chrom, pos, accession

    try:
        with open(output_filename, 'w') as file:
            file.write("##fileformat=VCFv4.0\n")
            file.write("##fileDate=20250211\n")
            file.write("##source=Kronos\n")
            file.write("##reference=Kronosv1.1\n")
            file.write("##phasing=full\n")
            file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            for (chrom, pos, accession), (line, confidence, param, label, desc) in sorted_variants:
                file.write(f"{line.strip()};{desc};Confidence={confidence};Threshold={param};Mapping={label}\n")

        print(f"Final VCF file written: {output_filename}")
    except IOError as e:
        print(f"Error writing output file: {e}")

def main():
    params = read_parameters('TSVs/No_RH/Mutations.summary')

    # Define all output files and corresponding VCF locations
    vcf_configs = [
        # Substitutions, No_RH only
        ('Kronos_v1.1.Exom-capture.corrected.deduped.10kb_bins.RH.byContig.MI.No_RH.maps.substitutions.vcf',
         [('TSVs/No_RH/', 'Unique', 'Type=Substitution;RH=False;Reference=Kronosv1.1'),
          ('TSVs-Multi/No_RH/', 'Multi', 'Type=Substitution;RH=False;Reference=Kronosv1.1')]),

        # Substitutions, RH_only
        ('Kronos_v1.1.Exom-capture.corrected.deduped.10kb_bins.RH.byContig.MI.RH_only.maps.substitutions.vcf',
         [('TSVs/RH_only/', 'Unique', 'Type=Substitution;RH=True;Reference=Kronosv1.1'),
          ('TSVs-Multi/RH_only/', 'Multi', 'Type=Substitution;RH=True;Reference=Kronosv1.1')]),

        # Indels, No_RH only
        ('Kronos_v1.1.Exom-capture.corrected.deduped.10kb_bins.RH.byContig.MI.No_RH.maps.indels.vcf',
         [('TSVs/No_RH/', 'Unique', 'Type=Indel;RH=False;Reference=Kronosv1.1'),
          ('TSVs-Multi/No_RH/', 'Multi', 'Type=Indel;RH=False;Reference=Kronosv1.1')]),

        # Indels, RH_only
        ('Kronos_v1.1.Exom-capture.corrected.deduped.10kb_bins.RH.byContig.MI.RH_only.maps.indels.vcf',
         [('TSVs/RH_only/', 'Unique', 'Type=Indel;RH=True;Reference=Kronosv1.1'),
          ('TSVs-Multi/RH_only/', 'Multi', 'Type=Indel;RH=True;Reference=Kronosv1.1')])
    ]

    # Process each set of VCFs
    for output_vcf, vcf_locations in vcf_configs:
        vcfs = []

        for loc, label, desc in vcf_locations:
            if not os.path.exists(loc):
                print(f"Warning: Directory {loc} does not exist.")
                continue

            # Select only indel files for indel VCF outputs
            if 'indels' in output_vcf:
                vcfs.extend([(os.path.join(loc, v), label, desc) for v in os.listdir(loc) if v.endswith('vcf') and 'indel' in v])
            else:
                vcfs.extend([(os.path.join(loc, v), label, desc) for v in os.listdir(loc) if v.endswith('vcf') and 'indel' not in v])

        if not vcfs:
            print(f"Error: No valid VCF files found for {output_vcf}. Skipping...")
            continue

        print(f"Processing {output_vcf} with {len(vcfs)} VCF files...")

        variant_dict = process_vcf(params, vcfs)
        write_final_vcf(output_vcf, variant_dict)

    print("All VCF files processed successfully.")

if __name__ == '__main__':
    main()
