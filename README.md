# Zig Transformer

Weekend project to implement the encoder-decoder self-attention
transformer model from scratch (i.e. probably never :-p or maybe Mamba
at the rate things are moving).

## Progress

So far I have a slow and complicated byte-pair encoder for
tokenization. It encodes bytes (limited to ASCII at the moment) into
16-bit codes that index into a dictionary. It supports doing the
encoding in blocks because the time complexity is not great.


