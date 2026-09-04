#!/usr/bin/env bash
# Provenance gap (c) — hybrid assembler, version, parameters and polishing stages.
#
#   ./recover_c_assembler.sh [WORK_ROOT] [OUT_TXT]
#   ./recover_c_assembler.sh /Volumes/PortableSSD/HR_WORK
#
# Read-only. Every assembler leaves a log or a run-parameter file next to its output; this
# script looks for those, and where the log records the command line it prints it verbatim.
# It reports what it finds and does NOT infer an assembler from the contig naming alone.
set -euo pipefail

ROOT=${1:-/Volumes/PortableSSD/HR_WORK}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=${2:-$HERE/../../data/assembler_provenance_scan.txt}
[[ -d $ROOT ]] || { echo "directory non trovata: $ROOT" >&2; exit 1; }

: > "$OUT"
say()   { printf '%s\n' "$*" | tee -a "$OUT"; }   # one line, to screen and to file
block() { tee -a "$OUT"; }                        # piped block, to screen and to file

say "# scan di $ROOT — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
say ""

# --- 1. files that identify an assembler by their own existence -------------------------
say "## file firma per assemblatore"
declare -a SIGS=(
  "flye:flye.log" "flye:params.json" "flye:assembly_info.txt" "flye:22-plasmids"
  "unicycler:unicycler.log" "unicycler:001_best_spades_graph.gfa"
  "spades:spades.log" "spades:contigs.paths" "spades:input_dataset.yaml"
  "canu:*.report" "canu:canu.out"
  "raven:raven.log" "miniasm:*.paf.gz" "wtdbg2:*.ctg.lay.gz"
  "trycycler:cluster_001" "trycycler:trycycler.log"
  "hybracter:hybracter.log" "dragonflye:dragonflye.log"
  "medaka:medaka.log" "medaka:consensus_probs.hdf"
  "racon:racon.log" "pilon:*.pilon.changes" "pilon:pilon.log"
  "polypolish:polypolish.log" "nextpolish:nextpolish.log"
  "porechop:porechop.log" "filtlong:filtlong.log"
  "guppy:sequencing_summary.txt" "dorado:dorado.log" "nanoq:nanoq.log"
)
FOUND=0
for s in "${SIGS[@]}"; do
  tool=${s%%:*}; pat=${s#*:}
  while IFS= read -r hit; do
    say "  [$tool] $hit"; FOUND=1
  done < <(find "$ROOT" \( -name "$pat" \) 2>/dev/null | head -20)
done
(( FOUND )) || say "  nessun file firma trovato sotto $ROOT"
say ""

# --- 2. command lines recorded inside those logs ---------------------------------------
say "## command line registrate nei log"
while IFS= read -r log; do
  say "  --- $log"
  hits=$(grep -m 8 -a -iE '(^|[^a-z])(flye|unicycler|spades|canu|raven|medaka|racon|pilon|polypolish|nextpolish|trycycler|hybracter|dragonflye)([^a-z]|$).*(--|-[a-z])' \
         "$log" 2>/dev/null || true)
  if [[ -n $hits ]]; then printf '%s\n' "$hits" | sed 's/^/      /' | block
  else say "      (nessuna riga di comando riconoscibile)"; fi
done < <(find "$ROOT" -type f \( -name '*.log' -o -name 'params.json' -o -name '*.report' -o -name 'assembly_info.txt' \) 2>/dev/null | head -40)
say ""

# --- 3. version strings anywhere in the logs -------------------------------------------
say "## stringhe di versione"
vers=$(find "$ROOT" -type f \( -name '*.log' -o -name '*.report' -o -name 'params.json' \) 2>/dev/null \
       | head -60 | while IFS= read -r log; do
           grep -m 3 -a -oiE '(flye|unicycler|spades|canu|raven|medaka|racon|pilon|polypolish|guppy|dorado)[ _v-]*[0-9]+\.[0-9]+(\.[0-9]+)?' \
             "$log" 2>/dev/null || true
         done | sort -u)
if [[ -n $vers ]]; then printf '%s\n' "$vers" | sed 's/^/  /' | block
else say "  nessuna stringa di versione trovata"; fi
say ""

# --- 4. contig naming — reported as a hint, never as the answer -------------------------
say "## naming dei contig negli assembly (indizio, non prova)"
say "   flye -> contig_1 | unicycler -> 1 length=.. | spades -> NODE_1_length_.._cov_.."
say "   canu -> tig00000001 | raven -> Utg.. | wtdbg2 -> ctg1 | trycycler -> cluster_001"
find "$ROOT" -type f \( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \) 2>/dev/null | head -6 |
  while IFS= read -r f; do say "  $(basename "$f"): $(head -1 "$f")"; done
say ""
say "## esito"
say "Se il punto 2 ha stampato una command line, (c) e' chiuso: copiala in pipeline/00_assembly.md."
say "Se non ha stampato nulla, (c) resta un buco dichiarato: il naming del punto 4 e' un indizio"
say "sull'assemblatore, non la sua identificazione, e non va scritto nei Methods come tale."
echo "scritto: $OUT" >&2
