#!/bin/bash
# Refine R-vs-S SNV candidates of one HR pair with DNBSEQ short reads.
#
#   bash scripts/refine_pair.sh 80
#
# 1. copy the reference assembly used by the long-read scan (HR_WORK, read-only)
# 2. polish it with Polypolish using the short reads of the SAME arm
# 3. bwa mem of R and S short reads on the polished reference (+ markdup)
# 4. lift candidate coordinates onto the polished reference
# 5. bcftools mpileup (AD/ADF/ADR) -> AF per arm at candidates + genome-wide
#    differential scan; freebayes --pooled-continuous as cross-check
# Per-step wall times go to logs/<id>_timing.tsv.
set -euo pipefail

ID="${1:?usage: bash refine_pair.sh <pair_id>}"
SSD="${SSD:-/Volumes/PortableSSD}"
PROJ="${PROJ:-$SSD/HR_SNP_refinement}"
SCRIPTS="$(cd "$(dirname "$0")" && pwd)"
THREADS="${THREADS:-8}"
CANDIDATES="${CANDIDATES:?set CANDIDATES to the candidate SNV table (semicolon-separated)}"
ENVS="${ENVS:-$HOME/miniconda3/envs}"
HBSS="$ENVS/hbss/bin"                       # bwa, samtools, bcftools
POLYPOLISH="$ENVS/flye/bin/polypolish"
FREEBAYES="$ENVS/snippy460/bin/freebayes"
PY="$ENVS/bactgenomics/bin/python"          # scipy
MAPQ_MIN=20
BASEQ_MIN=20
MAX_DEPTH=100000
FB_MIN_AF=0.02
FB_MIN_ALT=3
FB_PAD=100

OUT="$PROJ/results/$ID"
TMP="$OUT/tmp"
LOGDIR="$PROJ/logs"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="$LOGDIR/${ID}_${STAMP}.log"
TIMING="$LOGDIR/${ID}_timing.tsv"
mkdir -p "$OUT/ref" "$TMP" "$LOGDIR"
exec > >(tee -a "$LOG") 2>&1

export PATH="$HBSS:$PATH"
STEP="init"
trap 'echo "[ERROR] pair $ID failed during step: $STEP (see $LOG)"' ERR

step() {                                   # step <name> <command...>
    STEP="$1"; shift
    local t0=$SECONDS
    echo "=== [$ID] $STEP  ($(date +%H:%M:%S))"
    "$@"
    printf "%s\t%s\t%s\n" "$ID" "$STEP" "$((SECONDS - t0))" >> "$TIMING"
}

# ---------- inputs ----------
SUMMARY="$SSD/HR_RESULTS/$ID/snp_summary.json"
[ -f "$SUMMARY" ] || { echo "missing $SUMMARY"; exit 1; }
REF_NAME="$(grep -o '"reference_assembly": "[^"]*"' "$SUMMARY" | cut -d'"' -f4)"
REF_ARM="${REF_NAME%.fasta}"               # e.g. 80S
ORIG_REF_SRC="$SSD/HR_WORK/$ID/$REF_NAME"
[ -f "$ORIG_REF_SRC" ] || { echo "missing $ORIG_REF_SRC"; exit 1; }

reads() {                                  # reads <arm> <1|2>
    local f
    f="$(ls "$SSD/${1}short/"*_"$2".fq.gz 2>/dev/null | grep -v '/\._' || true)"
    [ "$(printf '%s\n' "$f" | grep -c .)" -eq 1 ] || { echo "expected one *_$2.fq.gz in $SSD/${1}short: $f" >&2; exit 1; }
    printf '%s' "$f"
}
R_R1="$(reads "${ID}R" 1)"; R_R2="$(reads "${ID}R" 2)"
S_R1="$(reads "${ID}S" 1)"; S_R2="$(reads "${ID}S" 2)"
if [ "$REF_ARM" = "${ID}R" ]; then P_R1="$R_R1"; P_R2="$R_R2"; else P_R1="$S_R1"; P_R2="$S_R2"; fi

ORIG="$OUT/ref/${REF_ARM}.orig.fasta"
POL="$OUT/ref/${REF_ARM}.polished.fasta"
echo "pair $ID | reference $REF_NAME ($REF_ARM) | threads $THREADS"
echo "R: $R_R1 $R_R2"; echo "S: $S_R1 $S_R2"

# ---------- 1-2. polish ----------
polish() {
    cp "$ORIG_REF_SRC" "$ORIG"
    bwa index "$ORIG" 2> "$TMP/bwa_index_orig.log"
    bwa mem -t "$THREADS" -a "$ORIG" "$P_R1" > "$TMP/pol_1.sam" 2> "$TMP/bwa_pol_1.log"
    bwa mem -t "$THREADS" -a "$ORIG" "$P_R2" > "$TMP/pol_2.sam" 2> "$TMP/bwa_pol_2.log"
    "$POLYPOLISH" filter --in1 "$TMP/pol_1.sam" --in2 "$TMP/pol_2.sam" \
        --out1 "$TMP/pol_f1.sam" --out2 "$TMP/pol_f2.sam"
    "$POLYPOLISH" polish "$ORIG" "$TMP/pol_f1.sam" "$TMP/pol_f2.sam" > "$POL"
    rm -f "$TMP"/pol_*.sam
}
if [ -s "$POL" ]; then echo "polished reference exists, reusing"; else step polish polish; fi

index_polished() {
    bwa index "$POL" 2> "$TMP/bwa_index_pol.log"
    samtools faidx "$POL"
}
[ -f "$POL.bwt" ] && [ -f "$POL.fai" ] || step index_polished index_polished

# ---------- 3. map ----------
map_arm() {                                # map_arm <sample> <r1> <r2>
    local bam="$OUT/$1.short.bam"
    bwa mem -t "$THREADS" -R "@RG\tID:$1\tSM:$1\tPL:DNBSEQ" "$POL" "$2" "$3" 2> "$TMP/bwa_$1.log" \
        | samtools fixmate -m -u - - \
        | samtools sort -u -@ 4 -T "$TMP/sort_$1" - \
        | samtools markdup -@ 4 -f "$OUT/$1.markdup.txt" - "$bam"
    samtools index "$bam"
    samtools coverage "$bam" > "$OUT/$1.coverage.tsv"
}
for arm in R S; do
    s="${ID}${arm}"
    [ "$arm" = R ] && r1="$R_R1" r2="$R_R2" || r1="$S_R1" r2="$S_R2"
    [ -f "$OUT/$s.short.bam.bai" ] || step "map_$s" map_arm "$s" "$r1" "$r2"
done
BAMS=("$OUT/${ID}R.short.bam" "$OUT/${ID}S.short.bam")

# ---------- 4. liftover ----------
step liftover "$PY" "$SCRIPTS/liftover_sites.py" --candidates "$CANDIDATES" \
    --pair "$ID" --orig "$ORIG" --polished "$POL" --out "$OUT/candidates_lifted.tsv"

# ---------- 5. calling ----------
MPILEUP_OPTS=(-f "$POL" -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP
              -d "$MAX_DEPTH" -q "$MAPQ_MIN" -Q "$BASEQ_MIN")

candidate_calls() {
    awk -F'\t' 'NR>1 && $11!="NA"{print $10"\t"$11}' "$OUT/candidates_lifted.tsv" \
        | sort -k1,1 -k2,2n > "$TMP/cand_sites.tsv"
    bcftools mpileup "${MPILEUP_OPTS[@]}" -T "$TMP/cand_sites.tsv" "${BAMS[@]}" \
        -Ov -o "$OUT/candidates.mpileup.vcf" 2> "$TMP/mpileup_cand.log"
    awk -v p="$FB_PAD" -F'\t' '{s=$2-1-p; if(s<0)s=0; print $1"\t"s"\t"$2+p}' \
        "$TMP/cand_sites.tsv" > "$TMP/cand_pad.bed"
    if [ -s "$TMP/cand_pad.bed" ]; then
        "$FREEBAYES" -f "$POL" --pooled-continuous -F "$FB_MIN_AF" -C "$FB_MIN_ALT" \
            --min-coverage 20 -m "$MAPQ_MIN" -q "$BASEQ_MIN" -t "$TMP/cand_pad.bed" \
            "${BAMS[@]}" > "$OUT/candidates.freebayes.vcf"
    else                                   # no candidate lifted (e.g. multi-copy gene)
        echo "no lifted candidates: skipping freebayes"
        : > "$OUT/candidates.freebayes.vcf"
    fi
    "$PY" "$SCRIPTS/af_table.py" candidates --vcf "$OUT/candidates.mpileup.vcf" \
        --r-sample "${ID}R" --s-sample "${ID}S" --lifted "$OUT/candidates_lifted.tsv" \
        --freebayes "$OUT/candidates.freebayes.vcf" --out "$OUT/candidates_AF.tsv"
}
step candidates candidate_calls

genome_scan() {                            # contigs in parallel chunks
    local chunk_len n=0
    chunk_len=$(( $(awk '{s+=$2}END{print s}' "$POL.fai") / THREADS + 1 ))
    awk -v L="$chunk_len" '{for(s=1;s<=$2;s+=L){e=s+L-1; if(e>$2)e=$2; print $1":"s"-"e}}' \
        "$POL.fai" > "$TMP/chunks.txt"
    local pids=() pid
    while read -r reg; do
        n=$((n + 1))
        ( set -o pipefail
          bcftools mpileup "${MPILEUP_OPTS[@]}" -r "$reg" "${BAMS[@]}" 2> "$TMP/chunk_$n.log" \
            | bcftools view -i 'N_ALT>1' -Oz -o "$TMP/chunk_$(printf %03d $n).vcf.gz" ) &
        pids+=($!)
    done < "$TMP/chunks.txt"
    for pid in "${pids[@]}"; do
        wait "$pid" || { echo "mpileup chunk failed (see $TMP/chunk_*.log)"; return 1; }
    done
    bcftools concat -Oz -o "$OUT/genome.mpileup.vcf.gz" "$TMP"/chunk_*.vcf.gz
    tabix -p vcf "$OUT/genome.mpileup.vcf.gz"
    rm -f "$TMP"/chunk_*.vcf.gz
    "$PY" "$SCRIPTS/af_table.py" genome --vcf "$OUT/genome.mpileup.vcf.gz" \
        --r-sample "${ID}R" --s-sample "${ID}S" --out "$OUT/genome_differential.tsv"
}
step genome_scan genome_scan

rm -rf "$TMP"
echo "=== [$ID] done. Timing:"
awk -F'\t' -v id="$ID" '$1==id{printf "  %-16s %6.1f min\n", $2, $3/60}' "$TIMING"
