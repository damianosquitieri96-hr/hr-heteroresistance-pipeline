#!/usr/bin/env bash
# Fill file_1_name / file_1_md5 / file_2_name / file_2_md5 in ena_runs.tsv from the
# FASTQ files on external storage. Run this on the Mac with the SSD mounted; it is the
# only stage that needs the reads themselves. No credentials are read or written.
#
#   ./fill_ena_runs.sh <READ_DIR> [ena_runs.tsv] [out.tsv]
#
# READ_DIR is searched recursively. A file is assigned to a run by the token
# <pair><arm> (e.g. 70R) appearing in its basename, and to a library by its name:
#   *_R1*/*_1.fastq*  -> Illumina forward      *_R2*/*_2.fastq*  -> Illumina reverse
#   anything matching ont|nanopore|long        -> ONT
#   any remaining single file for that sample  -> ONT
# Samples whose files cannot be resolved are listed on stderr and left as TO_FILL:
# the script never guesses a filename.
set -euo pipefail

READ_DIR=${1:?usage: fill_ena_runs.sh <READ_DIR> [in.tsv] [out.tsv]}
IN=${2:-"$(dirname "$0")/ena_runs.tsv"}
OUT=${3:-"$(dirname "$0")/ena_runs.filled.tsv"}
[[ -d $READ_DIR ]] || { echo "no such directory: $READ_DIR" >&2; exit 1; }

if command -v md5 >/dev/null 2>&1; then md5sum_of() { md5 -q "$1"; }
elif command -v md5sum >/dev/null 2>&1; then md5sum_of() { md5sum "$1" | cut -d' ' -f1; }
else echo "neither md5 nor md5sum found" >&2; exit 1; fi

INDEX=$(mktemp); trap 'rm -f "$INDEX"' EXIT
find "$READ_DIR" -type f \
     \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) \
     -print > "$INDEX"
echo "FASTQ trovati in $READ_DIR: $(wc -l < "$INDEX" | tr -d ' ')" >&2

split_tabs() {                       # tab is IFS-whitespace, so `read` would collapse
  local s=$1; F=()                   # empty fields; split by hand instead.
  while [[ $s == *$'\t'* ]]; do F+=("${s%%$'\t'*}"); s=${s#*$'\t'}; done
  F+=("$s")
}

# pick_file <pair><arm> <ONT|ILMN> <1|2> -> prints one path, or nothing.
# The token must be delimited by a non-alphanumeric (or the start/end of the basename)
# so that token 40R does not also match 240R.
pick_file() {
  local tok=$1 tag=$2 mate=$3 cands
  cands=$(grep -iE "/([^/]*[^0-9A-Za-z])?${tok}([^0-9A-Za-z][^/]*)?\$" "$INDEX" || true)
  [[ -z $cands ]] && return 0
  if [[ $tag == ILMN ]]; then
    printf '%s\n' "$cands" | grep -iE "(_R${mate}[._]|_${mate}\.f(ast)?q)" || true
  else
    local ont
    ont=$(printf '%s\n' "$cands" | grep -iE "(ont|nanopore|long|guppy|dorado)" || true)
    if [[ -n $ont ]]; then printf '%s\n' "$ont"
    else printf '%s\n' "$cands" | grep -ivE "(_R[12][._]|_[12]\.f(ast)?q)" || true
    fi
  fi
}

{
  header_done=0
  while IFS= read -r line; do
    # pass through the leading '#' directives and the column-header row untouched
    if [[ $header_done -eq 0 ]]; then
      printf '%s\n' "$line"
      [[ $line != \#* ]] && header_done=1
      continue
    fi
    split_tabs "$line"
    run_alias=${F[0]}; exp_alias=${F[1]}; fmt=${F[2]}
    f1=${F[3]}; m1=${F[4]}; f2=${F[5]}; m2=${F[6]}; layout=${F[7]}
    tok=${exp_alias#HR}; tag=${tok##*_}; tok=${tok%%_*}   # HR70R_ILMN -> tok=70R tag=ILMN
    for mate in 1 2; do
      [[ $mate == 2 && $layout != PAIRED ]] && continue
      hits=$(pick_file "$tok" "$tag" "$mate"); n=$(printf '%s' "$hits" | grep -c . || true)
      if [[ $n -eq 1 ]]; then
        b=$(basename "$hits"); s=$(md5sum_of "$hits")
        if [[ $mate == 1 ]]; then f1=$b; m1=$s; else f2=$b; m2=$s; fi
      else
        echo "IRRISOLTO ${run_alias} mate${mate}: ${n} candidati per token ${tok}/${tag}" >&2
      fi
    done
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$run_alias" "$exp_alias" "$fmt" "$f1" "$m1" "$f2" "$m2" "$layout"
  done < "$IN"
} > "$OUT"

echo "scritto $OUT" >&2
grep -c 'TO_FILL' "$OUT" >/dev/null 2>&1 && \
  echo "righe con campi ancora TO_FILL: $(grep -c 'TO_FILL' "$OUT")" >&2 || true
