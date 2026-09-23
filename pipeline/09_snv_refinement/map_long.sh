#!/bin/bash
# Map CycloneSEQ long reads of one or more arms onto the polished reference of a pair.
#
#   bash scripts/map_long.sh 254 R S      -> results/254/long/254R.long.bam, 254S.long.bam
#
# Uses the same polished reference as refine_pair.sh, so coordinates match the
# short-read results (candidates_AF.tsv, genome_differential.tsv).
set -euo pipefail

ID="${1:?usage: bash map_long.sh <pair_id> <arm>...}"; shift
[ "$#" -gt 0 ] || { echo "usage: bash map_long.sh <pair_id> <arm>..."; exit 1; }
SSD="${SSD:-/Volumes/PortableSSD}"
PROJ="${PROJ:-$SSD/HR_SNP_refinement}"
THREADS="${THREADS:-8}"
export PATH="${ENVS:-$HOME/miniconda3/envs}/hbss/bin:$PATH"

OUT="$PROJ/results/$ID/long"
LOG="$PROJ/logs/${ID}_long_$(date +%Y%m%d_%H%M%S).log"
TIMING="$PROJ/logs/${ID}_timing.tsv"
mkdir -p "$OUT"
exec > >(tee -a "$LOG") 2>&1

refs=("$PROJ/results/$ID/ref/"*.polished.fasta)
[ "${#refs[@]}" -eq 1 ] && [ -f "${refs[0]}" ] || { echo "need exactly one polished reference for $ID"; exit 1; }
REF="${refs[0]}"

for arm in "$@"; do
    s="${ID}${arm}"
    fq="$SSD/$s/1.Cleandata/$s.filtered_reads.fq.gz"
    [ -f "$fq" ] || { echo "missing $fq"; exit 1; }
    bam="$OUT/$s.long.bam"
    if [ -f "$bam.bai" ]; then echo "$bam exists, skipping"; continue; fi
    t0=$SECONDS
    echo "=== [$ID] long map $s -> $(basename "$REF")  ($(date +%H:%M:%S))"
    minimap2 -ax map-ont -t "$THREADS" -R "@RG\tID:${s}lr\tSM:$s" "$REF" "$fq" \
        | samtools sort -@ 4 -T "$OUT/sort_$s" -o "$bam"
    samtools index "$bam"
    samtools coverage "$bam" > "$OUT/$s.long.coverage.tsv"
    printf "%s\t%s\t%s\n" "$ID" "long_map_$s" "$((SECONDS - t0))" >> "$TIMING"
done
echo "=== [$ID] long mapping done"
