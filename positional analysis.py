# positional analysis

# aligned FASTA txt file to be analyzed
# file_path = r"C:\Users\carte\TCR Optimize\TCR Analysis\teb_alpha_aligned.txt"
# original_seq = "IQKPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDT"


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



def report_differences_with_composition(
    original_seq, 
    consensus_seq, 
    positionwise_percents
):
    """
    Compare each position of original_seq and consensus_seq.
    Return (and print) positions where there's a mismatch,
    along with the percent composition at that index.
    
    :param original_seq: String, the original amino acid sequence
    :param consensus_seq: String, the consensus sequence
    :param positionwise_percents: List[Dict[str, float]] 
                                  from compute_positionwise_amino_acid_percents()
    :return: List of dicts, each capturing mismatch data, e.g.:
             [
               {
                 "position": 3,
                 "original_aa": "K",
                 "consensus_aa": "D",
                 "composition": {"K": 68.29, "T": 21.95, "P": 7.32, "N": 2.44}
               },
               ...
             ]
    """
    differences = []
    length = min(len(original_seq), len(consensus_seq))
    
    for i in range(length):
        orig_aa = original_seq[i]
        cons_aa = consensus_seq[i]
        
        if orig_aa != cons_aa:
            # Get composition dictionary for position i (0-based).
            # Make sure it doesn't go out of range if positionwise_percents is shorter.
            if i < len(positionwise_percents):
                composition_dict = positionwise_percents[i]
            else:
                composition_dict = {}
            
            differences.append({
                "position": i + 1,
                "original_aa": orig_aa,
                "consensus_aa": cons_aa,
                "composition": composition_dict
            })
    
    # Print results in the desired format:
    if not differences:
        print("No mismatches found.")
    else:
        print("Mismatches with composition data:")
        for diff in differences:
            pos = diff["position"]
            o = diff["original_aa"]
            c = diff["consensus_aa"]
            comp_dict = diff["composition"]
            
            # Build a sorted list of "AA: XX.XX%" in descending % order
            comp_sorted = sorted(comp_dict.items(), key=lambda x: x[1], reverse=True)
            comp_str = ", ".join(f"{aa}: {pct:.2f}%" for aa, pct in comp_sorted)
            
            # Print the combined info
            print(
                f"Position {pos} - "
                f"Original: {o}, Consensus: {c}. "
                f"Index Composition: {comp_str}"
            )
    
    return differences


# Example usage in main flow:
if __name__ == "__main__":
    file_path = r"C:\Users\carte\TCR Optimize\TCR Analysis\b57_alpha_aligned.txt"
    original_sequence = "IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDT"

    # 1) Extract sequences & compute
    amino_acid_sequences = extract_amino_acid_sequences(file_path)
    positionwise_percents = compute_positionwise_amino_acid_percents(amino_acid_sequences)

    # 2) Print the composition table (optional)
    print_positionwise_percents(positionwise_percents)

    # 3) Build the consensus
    consensus_sequence = create_frequency_consesus(positionwise_percents)
    print(f"\nFrequency consensus sequence: {consensus_sequence}\n")

    # 4) Compare to original sequence WITH composition details
    diffs = report_differences_with_composition(
        original_sequence, 
        consensus_sequence, 
        positionwise_percents
    )