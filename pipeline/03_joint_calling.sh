#!/usr/bin/env bash
# Stage 03 — joint two-sample calling, R and S against the REFERENCE ARM's assembly of the
# same pair. bcftools 1.24 + htslib 1.24. Recovered verbatim from the exported VCF headers.
# The reference arm is S in 8 of the 17 called pairs and R in the other 9; pass it as the
# second argument (data/pair_reference_arm.csv, data/joint_vcf_provenance.csv).
set -euo pipefail
PAIR=$1
REFARM=${2:?"secondo argomento: R o S, il braccio di riferimento della coppia (data/pair_reference_arm.csv)"}
REF="${PAIR}${REFARM}.fasta"
bcftools mpileup -q 20 -Q 20 -a AD -f "$REF" "${PAIR}R.bam" "${PAIR}S.bam" \
  | bcftools call -mv -Ou \
  | bcftools norm -f "$REF" -Oz -o "${PAIR}.RS.vcf.gz"
bcftools index -t "${PAIR}.RS.vcf.gz"
