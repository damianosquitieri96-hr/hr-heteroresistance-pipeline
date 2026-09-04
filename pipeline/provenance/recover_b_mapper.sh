#!/usr/bin/env bash
# Provenance gap (b) — long-read mapper and version, from the @PG lines of the BAM headers.
#
#   ./recover_b_mapper.sh [BAM_ROOT] [OUT_TSV]
#   ./recover_b_mapper.sh /Volumes/PortableSSD/HR_RESULTS
#
# Read-only on BAM_ROOT: it only reads headers. Writes one TSV (default in the repo's
# data/ directory) with one row per @PG program record found, so the mapper, its version
# and its full command line can be pasted into pipeline/02_mapping_dosage.sh.
set -euo pipefail

ROOT=${1:-/Volumes/PortableSSD/HR_RESULTS}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=${2:-$HERE/../../data/bam_PG_records.tsv}

command -v samtools >/dev/null || { echo "samtools non trovato nel PATH. Attiva un env che lo contiene e rilancia." >&2; exit 1; }
[[ -d $ROOT ]] || { echo "directory non trovata: $ROOT — monta l'SSD o passa il percorso giusto come primo argomento." >&2; exit 1; }

echo "cerco BAM sotto $ROOT ..." >&2
mapfile -t BAMS < <(find "$ROOT" -type f \( -name '*.bam' -o -name '*.cram' \) 2>/dev/null | sort)
echo "BAM/CRAM trovati: ${#BAMS[@]}" >&2
(( ${#BAMS[@]} )) || { echo "nessun BAM trovato: controlla il percorso." >&2; exit 1; }

{
  printf 'file\tpg_id\tpg_pn\tpg_vn\tpg_cl\n'
  for b in "${BAMS[@]}"; do
    samtools view -H "$b" 2>/dev/null | awk -v f="$(basename "$b")" '
      /^@PG/ {
        id=""; pn=""; vn=""; cl="";
        for (i = 2; i <= NF; i++) {
          tag = substr($i, 1, 3); val = substr($i, 4);
          if (tag == "ID:") id = val;
          else if (tag == "PN:") pn = val;
          else if (tag == "VN:") vn = val;
          else if (tag == "CL:") cl = val;
        }
        print f "\t" id "\t" pn "\t" vn "\t" cl;
      }'
  done
} > "$OUT"

echo >&2
echo "=== programmi distinti (PN, VN) ===" >&2
awk -F'\t' 'NR>1 && $3 != "" {print $3"\t"$4}' "$OUT" | sort -u >&2
echo >&2
echo "=== command line distinte ===" >&2
awk -F'\t' 'NR>1 && $5 != "" {print $5}' "$OUT" | sort -u >&2
echo >&2
echo "scritto: $OUT ($(($(wc -l < "$OUT") - 1)) record @PG)" >&2
echo "Incolla la riga del mapper (minimap2/bwa/ngmlr/winnowmap) in pipeline/02_mapping_dosage.sh." >&2
echo "Se un BAM non ha @PG, il mapper non e' recuperabile da quel file: va dichiarato, non dedotto." >&2
