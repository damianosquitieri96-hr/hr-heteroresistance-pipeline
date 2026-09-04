#!/usr/bin/env bash
# Stage 02 — long-read mapping of both arms onto the REFERENCE ARM's assembly of the same
# pair, then per-gene and 1-kb-window coverage.
# The reference arm is not the same in every pair: S in 22, R in 11
# (data/pair_reference_arm.csv). Pass it as the second argument.
# NOTE: the mapper and version actually used are in the BAM @PG header line and were
# not available when this repo was assembled (docs/PROVENANCE_GAPS.md, item b).
# Recover them with pipeline/provenance/recover_b_mapper.sh, then paste the line below.
set -euo pipefail
PAIR=$1
REFARM=${2:?"secondo argomento: R o S, il braccio di riferimento della coppia (data/pair_reference_arm.csv)"}
REF="${PAIR}${REFARM}.fasta"
for ARM in R S; do
  # <MAPPER> <VERSION> -x map-ont "$REF" "${PAIR}${ARM}.long.fastq.gz" \
  #   | samtools sort -o "${PAIR}${ARM}.bam"; samtools index "${PAIR}${ARM}.bam"
  samtools flagstat "${PAIR}${ARM}.bam" > "${PAIR}${ARM}.long.flagstat.txt"
  samtools bedcov "${PAIR}${REFARM}.genes.bed" "${PAIR}${ARM}.bam" \
    > "${PAIR}${ARM}.long_on_${PAIR}${REFARM}.genecov.tsv"
  bedtools makewindows -g <(cut -f1,2 "${REF}.fai") -w 1000 \
    | samtools bedcov - "${PAIR}${ARM}.bam" | gzip \
    > "${PAIR}${ARM}.long_on_${PAIR}${REFARM}.win1kb.tsv.gz"
done
