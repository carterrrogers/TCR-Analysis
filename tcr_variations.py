from itertools import combinations

# Define the beta constant sequence and mismatches
beta_constant = list("FPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQDPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRA")
                      
# Define mismatch positions and alternative residues
mismatch_info = [
    (3, 'K'),  # Position 4 (0-based index 3), mismatch is 'K'
    (5, 'T'),  # Position 6 (0-based index 5), mismatch is 'T'
    (16, 'R'),  # Position 17 (0-based index 16), mismatch is 'R'
    (46, 'Q'),  # Position 47 (0-based index 46), mismatch is 'Q'
    (50, 'S'),  # Position 51 (0-based index 50), mismatch is 'S'
    (56, 'Y'),  # Position 57 (0-based index 56), mismatch is 'Y'
    (68, 'C')   # Position 69 (0-based index 68), mismatch is 'C'
]

# Generate all combinations of changes
variations = []
for num_changes in range(len(mismatch_info) + 1):
    for changes in combinations(mismatch_info, num_changes):
        # Create a copy of the beta constant sequence
        variant = beta_constant[:]
        for position, new_residue in changes:
            variant[position] = new_residue  # Apply the mismatch change
        variations.append("".join(variant))

# Label variations
labeled_variations = [f"Variation {i + 1}: {seq}" for i, seq in enumerate(variations)]

# Write all variations to a file or print them
with open("variations.txt", "w") as f:
    f.write("\n".join(labeled_variations))
