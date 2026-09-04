#!/usr/bin/env bash
# Stage 04 — k-mer sketch clonality over the 66 genomes (both arms of all pairs).
#
# k and sketch size were not recorded when the stage was first run; they were recovered
# afterwards by re-running this stage over a grid of k and sketch size and finding the
# combination that reproduces data/HR_clonality_distance_matrix_66genomes.csv. Only
# k=21 -s 200000 does: 490 of the 496 testable pairwise distances are identical to the
# 6 significant digits at which the matrix was stored, the other 6 differ by <= 4.9e-9
# (rounding on intra-clone distances near 1e-4). Grid: data/mash_k_recovery_grid.csv;
# comparison figure: results/mash_k_recovery.png. Mash v2.3, default seed (-S 42).
set -euo pipefail
mash sketch -k 21 -s 200000 -o hr66 *.fasta
mash dist hr66.msh hr66.msh > mash_dist.tsv
