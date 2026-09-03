#!/usr/bin/env bash
# Stage 02 — long-read mapping of both arms onto the S assembly of the same pair,
# then per-gene and 1-kb-window coverage.
# NOTE: the mapper and version actually used are in the BAM @PG header line and were
# not available when this repo was assembled (docs/PROVENANCE_GAPS.md, item b).
# Recover with: samtools view -H <pair>R.bam | grep '^@PG'   and paste the line below.
set -euo pipefail
PAIR=$1; REF="${PAIR}S.fasta"
for ARM in R S; do
  # <MAPPER> <VERSION> -x map-ont "$REF" "${PAIR}${ARM}.long.fastq.gz" \
  #   | samtools sort -o "${PAIR}${ARM}.bam"; samtools index "${PAIR}${ARM}.bam"
  samtools flagstat "${PAIR}${ARM}.bam" > "${PAIR}${ARM}.long.flagstat.txt"
  samtools bedcov "${PAIR}S.genes.bed" "${PAIR}${ARM}.bam" \
    > "${PAIR}${ARM}.long_on_${PAIR}S.genecov.tsv"
  bedtools makewindows -g <(cut -f1,2 "${REF}.fai") -w 1000 \
    | samtools bedcov - "${PAIR}${ARM}.bam" | gzip \
    > "${PAIR}${ARM}.long_on_${PAIR}S.win1kb.tsv.gz"
done
