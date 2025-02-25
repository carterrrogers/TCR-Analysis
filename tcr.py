import torch
import esm

# Load a pre-trained ESM model


def read_sequences(seq):
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()  # Set model to evaluation mode

    batch_converter = alphabet.get_batch_converter()

    # Example TCR sequence
    # ("TCR_Example", "CASSLGTDTQYF")
    tcr_sequence = [seq]

    batch_labels, batch_strs, batch_tokens = batch_converter(tcr_sequence)

    # Run the model to get embeddings
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33])

    # Extract and print embeddings for further analysis
    tcr_embeddings = results["representations"][33]
    
    return tcr_embeddings.shape, tcr_embeddings
    
    
test = [("TCR_Example", "CASSLGTDTQYF"),
        ("Beta Chain","FPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQDPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRA"),
        ("Alpha Chain","SILRDNGVVEGVGTSKVRLIRPRDSGVGLKVHRTHIQDTLISNIHVGVLGLGDVSLGVEVSEQANRL"),
        ]

for i in test:
    print(read_sequences(i))  
    # Shape: (sequence_length, embedding_size)
