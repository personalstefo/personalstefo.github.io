import re

input_file = "HA_aligned.fasta"  # Replace with your filename
output_file = "HA_BEAST_ready.fasta"

with open(input_file, "r") as infile:
    lines = infile.readlines()

output_lines = []
for line in lines:
    if line.startswith(">"):
        match = re.search(r"_(\d{4})_(\d{2})_(\d{2})$", line.strip())
        if match:
            year, month, day = match.groups()
            new_date = f"{year}-{month}-{day}"
            new_header = re.sub(r"_(\d{4})_(\d{2})_(\d{2})$", f"_{new_date}", line.strip())
            output_lines.append(f"{new_header}\n")
        else:
            output_lines.append(line)
    else:
        output_lines.append(line)

with open(output_file, "w") as outfile:
    outfile.writelines(output_lines)

print(f"✅ BEAST-ready FASTA written to: {output_file}")
