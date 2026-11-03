import re

def rename_taxa_with_host_location_genotype(input_file, output_file):
    """
    Parses a FASTA file and renames each taxon using host, location, date, and genotype (e.g., H5N1) from the header.
    """
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    with open(output_file, 'w') as outfile:
        for line in lines:
            if line.startswith(">"):
                try:
                    parts = line.strip().split("|")
                    if len(parts) == 3:
                        _, description, date = parts
                    elif len(parts) == 2:
                        description, date = parts
                    else:
                        print(f"Skipping malformed line: {line.strip()}")
                        continue

                    # Extract strain string: e.g., A/duck/Vietnam/NCVD-KA260/2012(H5N1)
                    strain_match = re.search(r"A\/[^|]+?\(\w+\)", description)
                    if not strain_match:
                        print(f"Strain format not found: {line.strip()}")
                        continue

                    strain_full = strain_match.group(0)  # e.g. A/duck/Vietnam/.../2012(H5N1)

                    # Extract genotype: the part in parentheses
                    genotype_match = re.search(r"\((H\d+N\d+)\)", strain_full)
                    genotype = genotype_match.group(1) if genotype_match else "UnknownGenotype"

                    # Remove genotype from strain string
                    strain_clean = re.sub(r"\(H\d+N\d+\)", "", strain_full)

                    # Parse parts from cleaned strain
                    strain_parts = strain_clean.split("/")
                    host = strain_parts[1] if len(strain_parts) > 1 else "UnknownHost"
                    location = strain_parts[2] if len(strain_parts) > 2 else "UnknownLocation"

                    # Parse date from final field
                    date_parts = date.strip().split("-")
                    year = date_parts[0] if len(date_parts) > 0 else "UnknownYear"
                    month = date_parts[1] if len(date_parts) > 1 else "UnknownMonth"
                    day = date_parts[2] if len(date_parts) > 2 else "UnknownDay"

                    new_taxon = f">{host}_{location}_{year}_{month}_{day}_{genotype}"
                    outfile.write(new_taxon + "\n")

                except Exception as e:
                    print(f"Error processing line: {line.strip()} | {e}")
                    continue
            else:
                outfile.write(line)

# Example usage
input_file = 'Dataset1.fasta'
output_file = 'Dataset1_renamed2.fasta'
rename_taxa_with_host_location_genotype(input_file, output_file)
