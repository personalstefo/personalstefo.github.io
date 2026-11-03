def trim_stop_codon(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        seq_id = ""
        sequence = ""
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    # Trim last 3 nucleotides (stop codon)
                    trimmed_seq = sequence[:-3]
                    outfile.write(seq_id + "\n")
                    # Write wrapped lines (optional)
                    for i in range(0, len(trimmed_seq), 60):
                        outfile.write(trimmed_seq[i:i+60] + "\n")
                seq_id = line
                sequence = ""
            else:
                sequence += line
        # Final record
        if seq_id:
            trimmed_seq = sequence[:-3]
            outfile.write(seq_id + "\n")
            for i in range(0, len(trimmed_seq), 60):
                outfile.write(trimmed_seq[i:i+60] + "\n")

# Example usage:
input_fasta = "/Users/nova/Downloads/finalmolevol/PB1_aligned.fasta"
output_fasta = "/Users/nova/Downloads/finalmolevol/PB1_trimmed.fasta"
trim_stop_codon(input_fasta, output_fasta)
