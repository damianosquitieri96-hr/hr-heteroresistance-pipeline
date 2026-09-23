# SNV refinement outputs (stage 09)

Produced by [`pipeline/09_snv_refinement/`](../../pipeline/09_snv_refinement/) for pairs
80, 202, 216, 226, 236, 240, 246, 254, 275, 279. Derived tables only; BAMs, VCFs and
polished assemblies are not included. `af_R`/`af_S` are always the resistant/susceptible
arm; coordinates are per pair, on the polished assembly of the reference arm.

| File | Content |
|---|---|
| `candidates_input.csv` | input: long-read SNV candidates per pair (semicolon-separated), as used |
| `final_table.tsv` | final table, 37 rows: candidates + new differential sites + long-read checks |
| `summary_candidates.tsv` | one row per candidate: short-read AF per arm, Fisher p, freebayes, verdict |
| `summary_new_diffs.tsv` | genome-wide R-vs-S differential sites not among the candidates |
| `summary_new_diffs_annotated.tsv` | the same, with Bakta feature, codon effect and aa change |
| `summary_pairs.tsv` | per-pair counts |
| `pairs/<id>/candidates_lifted.tsv` | candidate coordinates lifted to the polished assembly |
| `pairs/<id>/candidates_AF.tsv` | per-candidate read counts per arm and strand |
| `pairs/<id>/genome_differential.tsv` | all differential sites of the genome scan |
| `pairs/<id>/<id>{R,S}.coverage.tsv`, `*.markdup.txt` | `samtools coverage` / `markdup` stats, short reads |
| `pairs/{216,246,254}/long/*.tsv` | long-read coverage, phasing (216 *mpl*, 246 *dacB*), per-copy counts (254 *bla*CTX-M) |

AF columns: `af` point estimate, `lo`/`hi` Wilson 95% CI, `altF`/`altR` alt reads on the
forward/reverse strand, `dp` ref + alt reads.

## `final_table.tsv` column key

The table is in Italian, as produced.

| Column | Meaning | Values |
|---|---|---|
| `origine` | where the site comes from | `candidato long-read` (input candidate), `nuovo (scan short-read)` (genome scan) |
| `prodotto`, `variante`, `effetto` | product, variant, effect | aa change or `contig:pos` for non-coding sites |
| `af_R`, `af_S` | AF (95% CI) | |
| `stato` | state | `fissata in R/S` (fixed), `subclonale in R/S`, `assente` (absent) |
| `esito` | outcome | `confermato` (confirmed; suffix `(short reads)`/`(long reads)` = which reads carry the evidence), `non confermato (artefatto ONT)` (long-read artefact), `debole` (suggestive tier only), `probabile artefatto` (likely artefact) |
| `pathway` | mechanism class | `circuito AmpC / riciclo PG` AmpC induction / peptidoglycan recycling; `beta-lattamasi` beta-lactamase; `efflusso (MexAB-OprM)` efflux; `regolazione porine (EnvZ/OmpR)` porin regulation; `involucro / LPS-ECA` envelope; `RNA polimerasi` RNA polymerase |
| `note` | curated long-read evidence | |

Verdicts in `summary_candidates.tsv`: `confermato` confirmed, `non confermato` not
confirmed, `non valutabile` not assessable (depth < 20, or site absent: 254 *bla*CTX-M, resolved with long reads),
`discordante` opposite direction, `confermato, AF diversa` confirmed with AF shifted > 0.2.
