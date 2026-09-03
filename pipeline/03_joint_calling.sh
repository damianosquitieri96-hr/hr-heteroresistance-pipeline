#!/usr/bin/env bash
# Stage 03 — joint two-sample calling, R and S against the S assembly of the same pair.
# bcftools 1.24 + htslib 1.24. Recovered verbatim from the exported VCF headers.
set -euo pipefail
PAIR=$1
bcftools mpileup -q 20 -Q 20 -a AD -f "${PAIR}S.fasta" "${PAIR}R.bam" "${PAIR}S.bam" \
  | bcftools call -mv -Ou \
  | bcftools norm -f "${PAIR}S.fasta" -Oz -o "${PAIR}.RS.vcf.gz"
bcftools index -t "${PAIR}.RS.vcf.gz"
