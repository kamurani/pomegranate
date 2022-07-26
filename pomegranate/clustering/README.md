# Clustering ideas


### TODO

- select graphs based on radius 
- get residues within motif 
- do SEQUENCE alignment of that motif region to get similarity score 
- i.e. project the 3d residues into 1d order: maybe just use each node's sequence position?
- maybe test a different way of aligning them? like ordering by distance to psite?
- the idea would be the seq alignment gives the "best" way to overlay two motifs together (could maybe enforce that psite is aligned? but also just let the alignment figure out the best one)
- alternatively, do 3d alignment and get similarity 'energy' score 
- this similarity score is repeated for every combination (or at least, enough pairs)
- this is used to train embeddings of each psite (center node) when graph convolution is applied to it
