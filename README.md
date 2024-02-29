# Zig Transformer

Weekend project to implement the encoder-decoder self-attention
transformer model from scratch (i.e. probably never :-p or maybe Mamba
at the rate things are moving).

## Progress

1. So far I have a slow and complicated byte-pair encoder for
   tokenization. It encodes bytes (limited to ASCII at the moment) into
   16-bit codes that index into a dictionary. It supports doing the
   encoding in blocks because the time complexity is not great.

2. After looking at OpenAI's Tiktokenizer and watching Kaparthy's video on BPE
   I have a much more positive opinion of my own attempt. Here are
   some notes:
       - Training and encoding BPE have poor time complexity which has
         to be worked around. There could be some really clever
         solution, but it's not what people are presently using.
       - Generally a hand written regex is used to split strings to
         avoid the `m` or the `n` in `O(mn)` being large and avoid
         some tokens. I'm not sure whether a sequence of bytes can be
         engineered to cause a timeout with existing libraries and
         encodings. I certainly had this issue with training my own
         BPE with no splitting.
       - UTF-8 is used sort of. The regex works on UTF-8, but then the
         BPE works on bytes. So I observed that the first two bytes in
         a three byte character can be joined into a single token
         while the remaining byte is encoded as a seperate token. Rare
         4-byte UTF-8 chars can end up being 4 tokens where each
         token is represented with 4-bytes.
   Some things I'd like to try
       - [ ] Get my own implementation to work with UTF-8
       - [ ] Get it to encode and decode cl100k_base (GPT4's tokenizer)
       - [ ] Make it fast
       - [ ] Create some training scheme so that it chooses byte pairs
             to split strings on instead of manually specifying a
             regex
       - [ ] Create a CLI
       - [ ] Create a Python binding


