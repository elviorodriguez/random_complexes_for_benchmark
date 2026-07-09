# Random sampling of protein complexes for MultimerMapper benchmark
This repository contains the code to randomly sample protein complexes from ComplexTab files (Complex Portal database) to create a benchmark dataset for MultimerMapper. A fixed random seed was introduced in the code to ensure reproducibility.

# Code organization
ComplexTab files correspond to the 2026-01-14 12:01 update of the Complex Portal database.

```
# ComplexTab files
Ce_6239.tsv    # C. elegans
Dm_7227.tsv    # D. melanogaster
Ec_83333.tsv   # E. coli
Hs_9606.tsv    # H. sapiens
Sc_559292.tsv  # S. cerevisiae

# R code
complex_sampling.R

# Directory with sampled complexes ready to run with MultimerMapper
results

# Commands used to execute MultimerMapper
results/commands.sh
```

To reproduce the results, install the dependencies of the R code (libraries) and execute it. This will generate the results directory with the FASTA files used as MultimerMapper input.
