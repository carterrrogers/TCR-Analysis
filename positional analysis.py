# positional analysis

# aligned FASTA txt file to be analyzed
file_path = r"C:\Users\carte\TCR Optimize\TCR Analysis\fasta_aligned_teb_alpha.txt"



def extract_amino_acid_sequences(file_path):
    """
    Parses a FASTA-like text file to extract amino acid sequences.
    
    parameter: file_path - path to the input FASTA txt file
    return: list of amino acid sequences (each sequence is a single concatenated string)
    """
    sequences = []
    current_seq = []

    with open(file_path, 'r', encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line:
                # ignores blank lines
                continue

            if line.startswith(">"):
                # hit a new header line -> save the current sequence if any
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
            else:
                # this line is part of the sequence
                current_seq.append(line)

        # appends the last sequence in the file
        if current_seq:
            sequences.append("".join(current_seq))

    return sequences

from collections import defaultdict

def compute_positionwise_amino_acid_percents(sequences):
    """
    given a list of amino acid sequences, compute the percent composition
    of each amino acid at every position (column).

    param: sequences - list of strings, each a protein (amino acid) sequence.
    return: a list (length = max sequence length) of dictionaries,
             where result[i] = { amino_acid: percent, ... } for position i.
    """
    if not sequences:
        return []

    # determine the maximum length of the sequences
    max_len = max(len(seq) for seq in sequences)

    positionwise_percents = []

    # for each position i in [0..max_len-1],
    # figure out how many sequences have that position,
    # how often each amino acid occurs, etc.
    for i in range(max_len):
        # count how many times each amino acid appears at position i
        aa_counts = defaultdict(int)
        valid_seq_count = 0

        for seq in sequences:
            if len(seq) > i:
                aa_counts[seq[i]] += 1
                valid_seq_count += 1

        if valid_seq_count == 0:
            # no sequences have residue i => skip or store empty
            positionwise_percents.append({})
            continue

        # convert counts to percentages
        position_percent_dict = {}
        for aa, count in aa_counts.items():
            percent = (count / valid_seq_count) * 100
            position_percent_dict[aa] = percent

        positionwise_percents.append(position_percent_dict)

    return positionwise_percents


def print_positionwise_percents(positionwise_percents):
    """
    Pretty-print the position-wise amino acid composition.
    Skips amino acids with 0% composition.
    """
    for i, aa_dict in enumerate(positionwise_percents, start=1):
        if not aa_dict:
            # no sequences had this position
            print(f"Position {i}: (no residues)")
            continue

        # build a string of "AA X%" pairs, sorted descending %
        aa_pairs = sorted(
            aa_dict.items(),
            key=lambda x: x[1],
            reverse=True
        )
        # only display amino acids with > 0%
        aa_str = ", ".join(
            f"{aa}: {percent:.2f}%" for aa, percent in aa_pairs if percent > 0
        )
        print(f"Position {i}: {aa_str}")


def create_frequency_consesus(positionwise_percents):
    """
    given a list of dictionaries (positionwise_percents), where each element i is
    a dictionary mapping amino_acid -> percent at position i, build a consensus
    sequence by choosing the amino acid with the highest percentage at each position.

    param: positionwise_percents - List[Dict[str, float]]
    return: a single consensus sequence as a string.
    """
    consensus = []
    for aa_dict in positionwise_percents:
        if not aa_dict:
            # If there's no data for this position, place 'X' or another placeholder
            consensus.append('X')
        else:
            # Pick the amino acid with the highest percentage
            best_aa, best_percent = max(aa_dict.items(), key=lambda x: x[1])
            consensus.append(best_aa)
    return "".join(consensus)



if __name__ == "__main__":
    # usage:

    # 1) compute the position-wise composition:
    amino_acid_sequences =  extract_amino_acid_sequences(file_path)
    positionwise_percents = compute_positionwise_amino_acid_percents(amino_acid_sequences)

    # 2) print the result:
    print_positionwise_percents(positionwise_percents)
    print(f"Frequency consensus sequence: {create_frequency_consesus(positionwise_percents)}")



