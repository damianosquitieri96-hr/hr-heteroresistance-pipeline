#!/usr/bin/env bash
# Stage 04 — k-mer sketch clonality over the 66 genomes (both arms of all pairs).
# WARNING: k and sketch size are NOT recoverable (docs/PROVENANCE_GAPS.md, item a).
# The two flags below are therefore left as placeholders and must not be presented as
# the parameters used. Re-running this and reproducing
# data/HR_clonality_distance_matrix_66genomes.csv identifies the k that was used.
set -euo pipefail
mash sketch -k <K> -s <SKETCH_SIZE> -o hr66 *.fasta
mash dist hr66.msh hr66.msh > mash_dist.tsv
