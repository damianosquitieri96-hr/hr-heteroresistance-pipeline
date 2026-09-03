#!/usr/bin/env bash
# Stage 01 — annotation. Bakta v1.12.1, database v6.0 light (doi:10.1099/mgen.0.000685).
set -euo pipefail
PAIR=$1; ARM=$2            # e.g. 202 S
bakta --db "$BAKTA_DB" --output "annot/${PAIR}${ARM}" --prefix "${PAIR}${ARM}" \
      --threads 8 "${PAIR}${ARM}.fasta"
# gene intervals used for coverage aggregation in stage 02
awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1,$2-1,$3,$6}' "annot/${PAIR}${ARM}/${PAIR}${ARM}.tsv" \
  > "${PAIR}${ARM}.genes.bed"
# AMR determinants (version string not recoverable from the exports)
amrfinder -n "${PAIR}${ARM}.fasta" > "${PAIR}${ARM}.amrfinder.tsv"
