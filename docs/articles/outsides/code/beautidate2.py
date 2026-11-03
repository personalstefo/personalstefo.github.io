import re

input_file = "HA_BEAST_ready.fasta"  # Change this if needed
output_file = "HA_with_dates.fasta"

with open(input_file, "r") as infile:
    lines = infile.readlines()

output_lines = []
for line in lines:
    if line.startswith(">"):
        parts = line.strip().rsplit("_", 1)  # split once from the right
        if len(parts) == 2 and re.match(r"\d{4}-\d{2}-\d{2}$", parts[1]):
            date = parts[1]
            if not line.strip().endswith(f"|{date}"):
                line = line.strip() + f"|{date}\n"
        output_lines.append(line)
    else:
        output_lines.append(line)

with open(output_file, "w") as outfile:
    outfile.writelines(output_lines)

print(f"✅ Saved BEAST2-compatible FASTA to: {output_file}")
