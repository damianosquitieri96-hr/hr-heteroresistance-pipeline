#!/usr/bin/env bash
# Provenance gap (e) — the 16 pairs with no joint two-sample call set.
#
#   ./recover_e_jointvcf.sh --dry-run                 # default: only report what it finds
#   ./recover_e_jointvcf.sh --run [OUTDIR]
#   ./recover_e_jointvcf.sh --run ~/Desktop/hr_jointvcf_16 --short
#
# Reproduces the command recovered from the headers of the 17 existing VCFs, verbatim, so
# the new call sets are comparable with them. Reads from the SSD, writes only into OUTDIR
# (never onto the SSD), and never deletes anything.
#
# WHETHER THIS IS WORTH RUNNING: none of the 16 pairs appears in any VCF-derived output
# (Table 2, intra-pair chromosomal distance, differential sites, variant effect), so no
# published number or figure changes if it is never run. It would extend the intra-pair
# distance from 17 pairs to 33. See docs/PROVENANCE_GAPS.md, item e.
set -euo pipefail

WORK=${WORK:-/Volumes/PortableSSD/HR_WORK}
RES=${RES:-/Volumes/PortableSSD/HR_RESULTS}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REFTAB=$HERE/../../data/pair_reference_arm.csv

MODE=dry
OUTDIR=$HOME/Desktop/hr_jointvcf_16
USE_SHORT=0
while (( $# )); do
  case $1 in
    --dry-run) MODE=dry ;;
    --run) MODE=run ;;
    --short) USE_SHORT=1 ;;
    -*) echo "opzione non riconosciuta: $1" >&2; exit 2 ;;
    *) OUTDIR=$1 ;;
  esac
  shift
done

PAIRS=(70 72 76 87 103 137 155 173 189 197 200 211 219 220 231 257)

[[ -f $REFTAB ]] || { echo "manca $REFTAB" >&2; exit 1; }
if [[ $MODE == run ]]; then
  command -v bcftools >/dev/null || { echo "bcftools non trovato nel PATH." >&2; exit 1; }
  v=$(bcftools --version | head -1)
  echo "bcftools presente: $v" >&2
  grep -q '1\.24' <<<"$v" || echo "ATTENZIONE: i 17 VCF esistenti sono bcftools 1.24+htslib-1.24. Una versione diversa va dichiarata nei Methods." >&2
  mkdir -p "$OUTDIR"
fi
(( USE_SHORT )) && cat >&2 <<'W'
ATTENZIONE --short: i 17 call set esistenti sono stati fatti sui BAM long-read (<pair><arm>.bam).
Usando i BAM short-read le distanze NON sono confrontabili con quelle pubblicate e vanno
riportate come serie separata, non unite alle 17.
W

refarm_of() { awk -F, -v p="$1" 'NR>1 && $1==p {print $2; exit}' "$REFTAB"; }

printf '%-6s %-4s %-42s %-28s %-28s %s\n' pair ref reference_fasta bam_R bam_S stato >&2
for p in "${PAIRS[@]}"; do
  arm=$(refarm_of "$p")
  [[ -n $arm ]] || { printf '%-6s %-4s %s\n' "$p" "?" "braccio di riferimento assente dalla tabella" >&2; continue; }
  ref=$WORK/$p/${p}${arm}.fasta
  if (( USE_SHORT )); then
    bR=$(find "$RES/$p" -maxdepth 1 -name "${p}R.short_on_*.bam" 2>/dev/null | head -1)
    bS=$(find "$RES/$p" -maxdepth 1 -name "${p}S.short_on_*.bam" 2>/dev/null | head -1)
  else
    bR=$RES/$p/${p}R.bam
    bS=$RES/$p/${p}S.bam
  fi
  miss=""
  [[ -f $ref ]] || miss+=" reference"
  [[ -n ${bR:-} && -f ${bR:-} ]] || miss+=" bam_R"
  [[ -n ${bS:-} && -f ${bS:-} ]] || miss+=" bam_S"
  if [[ -n $miss ]]; then
    printf '%-6s %-4s %s\n' "$p" "$arm" "MANCA:$miss" >&2
    continue
  fi
  printf '%-6s %-4s %-42s %-28s %-28s %s\n' \
    "$p" "$arm" "$(basename "$ref")" "$(basename "$bR")" "$(basename "$bS")" \
    "$([[ $MODE == run ]] && echo "chiamo" || echo "pronto")" >&2
  [[ $MODE == run ]] || continue

  out=$OUTDIR/${p}.RS.vcf.gz
  if [[ -f $out ]]; then echo "    esiste gia', salto: $out" >&2; continue; fi
  bcftools mpileup -q 20 -Q 20 -a AD -f "$ref" "$bR" "$bS" \
    | bcftools call -mv -Ou \
    | bcftools norm -f "$ref" -Oz -o "$out"
  bcftools index -t "$out"
  bcftools stats "$out" > "$OUTDIR/${p}.bcfstats.txt"
done

cat >&2 <<EOF

Modo: $MODE${MODE:+ }$([[ $MODE == dry ]] && echo "(nessun file scritto; rilancia con --run per chiamare)")
Output: $OUTDIR
Dopo il run: copia i VCF nella struttura per coppia che recount.py si aspetta, punta la sua
variabile E a quella directory e rilancia lo stage 05. La distanza passa da 17 a 33 coppie:
va ricalcolata per tutte, non concatenata, e i Methods vanno aggiornati con il numero nuovo.
EOF
