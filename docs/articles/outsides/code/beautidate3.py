import re

input_file = "HA_BEAST_ready.fasta"         # Replace with your actual FASTA filename
output_file = "HA_with_slash.fasta"

with open(input_file, "r") as infile:
    lines = infile.readlines()

output_lines = []
for line in lines:
    if line.startswith(">"):
        # Replace date format in the last part of the header
        line = re.sub(r"(\d{4})-(\d{2})-(\d{2})", r"\1/\2/\3", line.strip())
        output_lines.append(f"{line}\n")
    else:
        output_lines.append(line)

with open(output_file, "w") as outfile:
    outfile.writelines(output_lines)

print(f"✅ Converted FASTA saved to: {output_file}")
