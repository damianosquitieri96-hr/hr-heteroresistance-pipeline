# Provenance gaps

Parameters of the pipeline as run that could **not** be recovered from the surviving
files. They are listed here rather than reconstructed, and the Methods text points to
this file at each affected stage.

**a. Mash sketch parameters (k, sketch size).** — **CLOSED 2026-09-04.** The intermediate
`mash_dist.tsv` is still missing, but the parameters were identified by re-running
`mash sketch` + `mash dist` over a grid of 11 k (13–31) × 7 sketch sizes (1 000–200 000)
and comparing each result with the surviving matrix. **k = 21, sketch size 200 000**
(Mash v2.3, default seed) is the only combination that reproduces
`data/HR_clonality_distance_matrix_66genomes.csv`: 490 of the 496 testable distances are
identical to the 6 significant digits at which the matrix was stored, and the 6 remaining
residuals are ≤ 4.9e-9, all on intra-clone distances near 1e-4 — rounding at the last
stored digit. No other combination yields a single exact distance.
`pipeline/04_clonality_mash.sh` now carries the flags instead of placeholders. Grid:
`data/mash_k_recovery_grid.csv`; figure: `results/mash_k_recovery.png`.

Caveat, stated because it is part of the evidence: the test used the 32 of 66 assemblies
present on local storage (496 of 2 145 distances). `mash dist` is pairwise and a sketch
does not depend on which other genomes are sketched with it, so this does not weaken the
identification, but the remaining 1 649 distances have not been checked against the
recovered parameters. Re-run the stage on all 66 assemblies when the external drive is
mounted to close that too.

Clonality remains reported as sketch distance only — never as SNP counts and never against
published clonality thresholds. That choice was never a consequence of the missing
parameters; it follows from the thresholds addressing inter-patient transmission rather
than two isolates from one episode.

Each open item below has a recovery script in `pipeline/provenance/`. The scripts read the
external drive and are read-only on it; only the one for item (e) writes, and only into a
directory given on the command line. Run them with the drive mounted:

| Item | Script | What closes it |
|---|---|---|
| b | `recover_b_mapper.sh` | a `@PG` line naming the mapper and its version |
| c | `recover_c_assembler.sh` | a command line in an assembler log |
| d | `recover_d_species_db.sh` | tool + database release, from the classifier's own output |
| e | `recover_e_jointvcf.sh` | 16 new joint call sets, made with the recovered command |

**b. Long-read mapper and version.** Recorded in the `@PG` line of the BAM headers on
external storage, which was not mounted when this repository was assembled. The BAM paths
are known exactly, from the `##bcftoolsCommand` header of the exported VCFs:
`/Volumes/PortableSSD/HR_RESULTS/<pair>/<pair>{R,S}.bam`. Run
`pipeline/provenance/recover_b_mapper.sh`, which writes `data/bam_PG_records.tsv`, then
paste the mapper line into `pipeline/02_mapping_dosage.sh`. If a BAM turns out to carry no
`@PG` record, the mapper is not recoverable from that file and this item stays open.

**c. Hybrid assembly.** Assembler, version, parameters and polishing stages unknown.
The assembly FASTA headers are generic (`>contig_1`), which is consistent with Flye's
default contig naming but is an indication only — Flye is **not** written into the Methods
on that basis. `pipeline/provenance/recover_c_assembler.sh` looks for the log or
run-parameter file that every assembler leaves beside its output (`flye.log`, `params.json`,
`assembly_info.txt`, `unicycler.log`, `spades.log`, a Canu `.report`, and the polisher logs)
and prints any command line it finds verbatim. One line per stage is needed before
`pipeline/00_assembly.md` can be considered complete.

**d. Species-call database version.** `species.tsv` gives the calls but not the tool or the
reference database release used to produce them; the file itself is three bare columns
(pair, genus, species) with no header and no tool signature, so it must be identified from
the classifier's own output. `pipeline/provenance/recover_d_species_db.sh` searches for
those outputs and for any line declaring a database path or release, and lists the conda
environments holding a classifier. Related: the genus-level discordance on **pair 219**
(*E. coli* in `species.tsv`, *C. koseri* in the clinical database) is unresolved and blocks
the ENA deposition of that sample — see `submission/ena/README.md`.

**e. Pairs without a joint VCF.** 16 of 33 pairs have no joint
two-sample VCF (70, 72, 76, 87, 103, 137, 155, 173, 189, 197, 200, 211, 219, 220, 231, 257) and therefore no intra-pair chromosomal distance.

*Scope of the consequence, checked 2026-09-04.* None of these 16 pairs appears in any
VCF-derived output: `data/table2_snv.csv` (5 pairs), `data/HR_intrapair_chromosomal_distance.csv`
(17), `data/HR_intrapair_differential_sites.csv` (14) and `data/vep_missense.csv` (5) contain
none of them. No reported SNV, distance, variant effect or figure is missing a value or
would change if this item is never closed. 12 of the 16 do appear in
`data/HR_betalactam_dosage_ONT_15pairs.csv` and 2 in `data/HR_intraclone_matrices_x1e6.csv`,
but those are derived from read depth and from mash sketches respectively and use no VCF.
Closing this item would extend the intra-pair distance from 17 pairs to 33; it would not
correct anything. `pipeline/provenance/recover_e_jointvcf.sh` reproduces the recovered
bcftools command for these 16 pairs, using each pair's own reference arm and, by default,
the same long-read BAMs the existing 17 used — calling them from the short-read BAMs
instead would give distances that are not comparable with the published 17 and would have
to be reported as a separate series.

The core-genome counts previously reported for the UKE clone 70+76 (values 3 and 52) come
from the superseded windowed BLASTN method and are **not** carried into this repository;
the historical method is documented for transparency in
`docs/legacy_SNP_calling_method_HR163.md`.

## Corrections

**Reference arm — corrected 2026-09-04.** This repository previously stated, in
`README.md`, `docs/methods.md`, `pipeline/00_assembly.md`, `pipeline/02_mapping_dosage.sh`,
`pipeline/03_joint_calling.sh` and the docstring of `recount.py`, that the **susceptible**
arm's assembly was the within-pair reference at every stage. That is wrong for a third of
the cohort. The `##reference` header of the 17 exported VCFs gives the resistant arm for 9
of them, and the `_on_<pair><arm>` infix of the coverage filenames gives the reference arm
for all 33 pairs — the two agree on all 17 where both exist. The reference arm is the
susceptible one in 22 pairs and the resistant one in 11 (76, 87, 212, 226, 236, 240, 254,
265, 275, 279, 283); the table is `data/pair_reference_arm.csv`. Within each pair the same
arm is used for annotation, mapping, dosage and calling, so every pair sits in one
coordinate frame. `recount.py`'s docstring also attributed the choice of the susceptible
arm to a published recommendation; the recommendation is for a *within-pair* reference and
does not specify an arm.

**No result changes.** `recount.py` locates the two samples by their `R`/`S` suffix in the
VCF header and classifies on the maximum and minimum of the two alt fractions, so it is
symmetric in the arms and never assumes which one is the reference; no script in the
repository filters on `af_R` or `af_S`. The five pairs in `data/table2_snv.csv` (40, 80,
202, 216, 246) are all susceptible-arm-referenced, so the allele fractions reported there
as "afR" carry the intended meaning. What changed is the documentation, not a number.

**Interpretive caveat this creates.** In `data/HR_intrapair_differential_sites.csv`, `af_R`
and `af_S` are alt fractions relative to *that pair's* reference arm, so their polarity
flips between pairs: at fixed differences the median is af_R 0.99 / af_S 0.01 in
susceptible-referenced pairs and af_R 0.002 / af_S 0.998 in resistant-referenced ones.
Filtering that table on `af_R ≥ 0.9` to obtain "variants carried by the resistant arm"
would silently drop every resistant-referenced pair. Use the `cls` column together with
`data/pair_reference_arm.csv`.

## Remaining gaps

**f. Variant-effect scripts.** `esm_score.py` and `struct_analysis.py` are not present in
the project archive; only their output table (`data/vep_missense.csv`) and figure survive.
The scoring definition (ESM-2 650M, masked-marginal LLR at the site; RSA, pLDDT and
catalytic-site distance on the AlphaFold DB model) is stated in `docs/methods.md`, but the
exact script is not available.
