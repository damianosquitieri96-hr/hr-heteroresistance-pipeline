# Stage 09 — short-read refinement of R-vs-S SNV candidates

Re-examines, with DNBSEQ short reads, the within-pair SNV candidates first called from
CycloneSEQ long reads, in the 10 pairs where an SNV was the candidate HR mechanism:
**80, 202, 216, 226, 236, 240, 246, 254, 275, 279**. Pair 40 was not refined: it has no
short reads.

Heteroresistance mutations can be subclonal, so no genotype is called anywhere in this
stage: every site is reported as an **allele fraction per arm** (alt / (ref + alt) reads,
Wilson 95% CI).

Run on 2026-09-23. Outputs: [`../../data/snv_refinement/`](../../data/snv_refinement/).

## Coordinate frame

Each pair keeps the reference arm of the long-read scan (the `reference_assembly` of that
pair, see `data/pair_reference_arm.csv`), polished with the short reads **of the same arm**.
Both arms' short reads are mapped to that one polished assembly. Unlike
`data/HR_intrapair_differential_sites.csv`, the columns `af_R`/`af_S` here are always the
resistant and the susceptible arm, by sample name, whichever arm is the reference.
Coordinates are on the polished assembly (`new_contig`/`new_pos`); `orig_pos` is on the
original Flye assembly.

## Steps

| # | Script | What it does |
|---|---|---|
| 1 | `refine_pair.sh <pair>` | copy reference assembly → Polypolish (filter + polish) with the reference arm's short reads → `bwa mem` of R and S, `samtools fixmate/sort/markdup` → lift candidates (`liftover_sites.py`) → `bcftools mpileup` AD/ADF/ADR at candidates, `freebayes --pooled-continuous` as cross-check → genome-wide two-sample `mpileup` scan → `af_table.py` |
| — | `run_batch.sh <pair>...` | runs step 1 on several pairs in series |
| 2 | `summarize_batch.py` | merges per-pair outputs into `summary_candidates.tsv`, `summary_new_diffs.tsv` (differential sites not among the candidates), `summary_pairs.tsv` |
| 3 | `annotate_diffs.py` | lifts new differential sites back to the original assembly and places them on the Bakta GFF3 of the reference arm (codon effect, or flanking genes if intergenic) |
| 4 | `map_long.sh <pair> <arm>...` | `minimap2 -ax map-ont` of long reads onto the same polished reference (pairs 216, 246, 254) |
| 5 | `long_read_check.py phase\|copies` | haplotype counts of nearby sites on spanning reads; base counts per copy of a repeated gene |
| 6 | `final_table.py` | one table of candidates + new sites + long-read checks → `final_table.tsv` |

Step 5 was run by hand; the BAM and the `--site` arguments of every invocation are the
header lines of its output (`data/snv_refinement/pairs/<id>/long/*.tsv`). The long-read
interpretation of those outputs is written into `final_table.py` (dictionary `CURATED`),
not computed.

## Parameters

| Parameter | Value | Where |
|---|---|---|
| mapping / base quality in `mpileup` | MAPQ ≥ 20, BQ ≥ 20, max depth 100 000 | `refine_pair.sh` |
| freebayes | `--pooled-continuous -F 0.02 -C 3 --min-coverage 20`, ±100 bp around candidates | `refine_pair.sh` |
| minimum depth per arm | 20 (below: *non valutabile*) | `af_table.py` `MIN_DEPTH` |
| alt reads on both strands | ≥ 2 per strand | `MIN_ALT_PER_STRAND` |
| candidate test | Fisher exact R vs S, p < 1e-3 | `P_CANDIDATE` |
| genome scan | ≥ 5 alt reads, \|ΔAF\| ≥ 0.05, p < 1e-8; suggestive tier p < 1e-5 and AF ≥ 0.03 | `af_table.py` |
| fixed vs subclonal | \|ΔAF\| > 0.8 → fixed; `final_table.py`: AF ≥ 0.8 fixed, ≥ 0.05 present | `af_table.py`, `final_table.py` |
| long reads | primary alignments, MAPQ ≥ 20 | `long_read_check.py` `MAPQ_MIN` |

## Tool versions

bwa 0.7.19-r1273 · samtools 1.24 · bcftools 1.24 · minimap2 2.31-r1302 · Polypolish 0.7.1 ·
freebayes v1.3.6 · Python 3.14.6 (pandas 3.0.5, scipy 1.18.0). Annotation: the Bakta GFF3
already produced in stage 01.

## Environment and inputs

Scripts expect the drive layout they were run on (`$SSD/<id>R`, `$SSD/<id>Rshort`,
`$SSD/HR_WORK/<id>`, `$SSD/HR_RESULTS/<id>/snp_summary.json`); none of it is in this
repository (see *Data availability* in the top-level README).

| Variable | Default | Meaning |
|---|---|---|
| `SSD` | `/Volumes/PortableSSD` | root of the raw data (read only) |
| `PROJ` | `$SSD/HR_SNP_refinement` | output root (`results/`, `logs/`) |
| `CANDIDATES` | *(required)* | candidate table; the one used is `data/snv_refinement/candidates_input.csv` |
| `ENVS` | `$HOME/miniconda3/envs` | conda envs: `hbss` (bwa, samtools, bcftools, minimap2), `flye` (polypolish), `snippy460` (freebayes), `bactgenomics` (python) |
| `THREADS` | 8 | |

```bash
export CANDIDATES=data/snv_refinement/candidates_input.csv
bash pipeline/09_snv_refinement/run_batch.sh 80 202 216 226 236 240 246 254 275 279
python pipeline/09_snv_refinement/summarize_batch.py
python pipeline/09_snv_refinement/annotate_diffs.py
bash pipeline/09_snv_refinement/map_long.sh 254 R S   # likewise 216 R, 246 R
python pipeline/09_snv_refinement/final_table.py
```

## Results in brief

- All 10 primary candidates confirmed by short reads: 9 fixed, 1 subclonal
  (246 *dacB* N312K, AF 0.71 in R).
- Every long-read-only subclonal candidate (pairs 80 and 216, 10 sites) is absent from
  the short reads of both arms: CycloneSEQ artefacts.
- AmpC induction / peptidoglycan recycling (*ampR*, *ampD*, *dacB*, *mpl*) is hit in 6 of
  10 pairs.
- Convergent, mutually exclusive subclones within R: 216 *mpl* G29S vs a frameshift
  (short-read haplotypes 63/62/0/0); 246 *dacB* N312K (~62%) vs P175L (~19%), 1 of 418
  spanning long reads carries both.
- 254: R carries 2 × *bla*CTX-M-15 + 1 × CTX-M-178 (chromosomal copy at 2.08 Mb),
  S 3 × CTX-M-15; fixed difference, no amplification.
- New subclonal sites in R: 240 *mpl* H190Y (AF 0.51), 226 *envZ* K136N (AF 0.16).

`final_table.tsv` (37 rows) is written in Italian, as it was produced; column key in
[`../../data/snv_refinement/README.md`](../../data/snv_refinement/README.md).
