#!/usr/bin/env bash
# Stage 02 — long-read mapping of both arms onto the REFERENCE ARM's assembly of the same
# pair, then per-gene and 1-kb-window coverage.
# The reference arm is not the same in every pair: S in 22, R in 11
# (data/pair_reference_arm.csv). Pass it as the second argument.
# Mapper recovered 2026-09-04 from the @PG records of the BAM headers
# (minimap2 2.31-r1302, samtools 1.24; data/bam_PG_records.tsv). Thread count was 5 or 6
# depending on the run and is the only parameter that varied; it does not affect output.
set -euo pipefail
PAIR=$1
REFARM=${2:?"secondo argomento: R o S, il braccio di riferimento della coppia (data/pair_reference_arm.csv)"}
REF="${PAIR}${REFARM}.fasta"
THREADS=${3:-5}
for ARM in R S; do
  minimap2 -t "$THREADS" -I 1G -K 100M -ax map-ont "$REF" "${PAIR}${ARM}.fastq.gz" \
    | samtools sort -@ 4 -m 512M -o "${PAIR}${ARM}.bam" -
  samtools index "${PAIR}${ARM}.bam"
  samtools flagstat "${PAIR}${ARM}.bam" > "${PAIR}${ARM}.long.flagstat.txt"
  samtools bedcov "${PAIR}${REFARM}.genes.bed" "${PAIR}${ARM}.bam" \
    > "${PAIR}${ARM}.long_on_${PAIR}${REFARM}.genecov.tsv"
  bedtools makewindows -g <(cut -f1,2 "${REF}.fai") -w 1000 \
    | samtools bedcov - "${PAIR}${ARM}.bam" | gzip \
    > "${PAIR}${ARM}.long_on_${PAIR}${REFARM}.win1kb.tsv.gz"
done
