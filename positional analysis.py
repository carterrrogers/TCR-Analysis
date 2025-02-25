# positional analysis

def extract_amino_acid_sequences(file_path):
    """
    Parses a FASTA-like text file to extract amino acid sequences.
    
    :param file_path: Path to the input text file.
    :return: List of amino acid sequences (each sequence is a single concatenated string).
    """
    sequences = []
    current_seq = []

    with open(file_path, 'r', encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line:
                # Ignore blank lines
                continue

            if line.startswith(">"):
                # We've hit a new header line -> save the current sequence if any
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
                # You could store the header text if needed, but here we just skip it
            else:
                # This line is part of the sequence
                current_seq.append(line)

        # Don't forget the last sequence in the file
        if current_seq:
            sequences.append("".join(current_seq))

    return sequences


# Example usage:
file_path = r"C:\Users\carte\Downloads\seqdump (10).txt"
amino_acid_sequences = extract_amino_acid_sequences(file_path)

# Print or inspect the sequences
for seq in amino_acid_sequences:
    print(seq)



