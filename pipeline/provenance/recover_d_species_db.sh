#!/usr/bin/env bash
# Provenance gap (d) — which tool made the species calls in species.tsv, and which
# reference database release it used.
#
#   ./recover_d_species_db.sh [WORK_ROOT] [OUT_TXT]
#   ./recover_d_species_db.sh /Volumes/PortableSSD/HR_WORK
#
# Read-only. species.tsv in the export is three bare columns (pair, genus, species) with no
# header and no tool signature, so the tool must be identified from its own output files.
set -euo pipefail

ROOT=${1:-/Volumes/PortableSSD/HR_WORK}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=${2:-$HERE/../../data/species_db_provenance_scan.txt}
[[ -d $ROOT ]] || { echo "directory non trovata: $ROOT" >&2; exit 1; }

: > "$OUT"
say()   { printf '%s\n' "$*" | tee -a "$OUT"; }
block() { tee -a "$OUT"; }
say "# scan di $ROOT — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
say ""

say "## output per tool di classificazione"
declare -a SIGS=(
  "kraken2:*.kreport" "kraken2:*.k2report" "kraken2:*.kraken2" "kraken2:kraken2.log"
  "bracken:*.bracken" "gtdbtk:gtdbtk.log" "gtdbtk:*.bac120.summary.tsv"
  "gtdbtk:gtdbtk.json" "mlst:*.mlst.tsv" "mlst:mlst.log"
  "kmerfinder:*.KmerFinder*" "kmerfinder:results.txt"
  "bakta:bakta.log" "bakta:*.bakta.json" "prokka:*.prokka.txt"
  "amrfinder:*.amrfinder.tsv" "amrfinder:amrfinder.log"
  "speciesfinder:*.SpeciesFinder*" "rmlst:*.rmlst*" "sylph:*.sylph.tsv"
  "checkm:checkm.log" "ani:*.ani.tsv" "fastani:*.fastani*"
)
FOUND=0
for s in "${SIGS[@]}"; do
  tool=${s%%:*}; pat=${s#*:}
  while IFS= read -r hit; do say "  [$tool] $hit"; FOUND=1; done \
    < <(find "$ROOT" -name "$pat" 2>/dev/null | head -10)
done
(( FOUND )) || say "  nessun output di classificazione trovato sotto $ROOT"
say ""

say "## righe che dichiarano un database o una release"
dbhits=$(find "$ROOT" -type f \( -name '*.log' -o -name '*.tsv' -o -name '*.json' -o -name '*.txt' \) 2>/dev/null \
  | head -400 | while IFS= read -r f; do
      hits=$(grep -m 2 -a -iE 'database[ _-]?(version|release|path|dir)|db[ _-](version|release|path|dir)|--db|GTDB[ -]?(r|release)[0-9]|refseq[ _-]?release|standard[-_][0-9]{8}|k2_[a-z]+_[0-9]{8}|pluspf|amrfinder.*database' \
             "$f" 2>/dev/null || true)
      [[ -n $hits ]] && { printf '  --- %s\n' "$f"; printf '%s\n' "$hits" | sed 's/^/      /'; }
      true
    done)
if [[ -n $dbhits ]]; then printf '%s\n' "$dbhits" | block
else say "  nessuna riga che dichiari un database"; fi
say ""

say "## versione del database AMRFinderPlus (nell'header dei suoi TSV)"
find "$ROOT" -name '*.amrfinder.tsv' 2>/dev/null | head -3 | while IFS= read -r f; do
  say "  $(basename "$f"): $(head -1 "$f" | cut -c1-200)"
done
command -v amrfinder >/dev/null && say "  amrfinder installato ora: $(amrfinder --version 2>/dev/null) / DB $(amrfinder --database_version 2>&1 | tail -1)"
say ""

say "## conda: env che contengono un classificatore, con versioni"
if command -v conda >/dev/null; then
  for e in $(conda env list | awk '!/^#/ && NF {print $1}'); do
    v=$(conda list -n "$e" 2>/dev/null | awk '/^(kraken2|gtdbtk|mlst|kmerfinder|bakta|ncbi-amrfinderplus|sylph|fastani)[ ]/ {print $1"="$2}' | paste -sd' ' -)
    [[ -n $v ]] && say "  $e: $v"
  done
else
  say "  conda non nel PATH"
fi
say ""
say "## esito"
say "Serve la coppia (tool+versione, release del database). Se il tool e' kraken2, la release"
say "e' il nome della directory passata a --db; se e' GTDB-Tk, e' la stringa GTDB rNN nel log;"
say "se e' bakta, e' il campo 'db version' del suo log. Senza una di queste, (d) resta aperto."

echo "scritto: $OUT" >&2
