# Stage 06 — predicted variant effect

Missense variants scored with ESM-2 650M (`esm2_t33_650M_UR50D`), log-likelihood ratio
between mutant and wild-type residue at the masked site (masked-marginal). Structural
context (RSA, per-residue pLDDT, distance from the catalytic site) computed on the
AlphaFold DB model of the same protein.

The scripts `esm_score.py` and `struct_analysis.py` are **not present** in the project
archive (`docs/PROVENANCE_GAPS.md`, item f). Surviving outputs: `../../data/vep_missense.csv`
and `../../results/HR_missense_impact.png`.
