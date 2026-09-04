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

**b. Long-read mapper and version.** Recorded in the `@PG` line of the BAM headers on
external storage (`/Volumes/PortableSSD/HR_RESULTS/`), which was not mounted when this
repository was assembled. Recover with
`samtools view -H <pair>R.bam | grep '^@PG'` and paste into `pipeline/02_mapping_dosage.sh`.

**c. Hybrid assembly.** Assembler, version, parameters and polishing stages unknown.
The assembly FASTA headers are generic (`>contig_1`) and carry no assembler signature.
One line per stage is needed before `pipeline/00_assembly.md` can be considered complete.

**d. Species-call database version.** `species.tsv` gives the calls but not the reference
database release used to produce them.

**e. Pairs without a joint VCF.** 16 of 33 pairs have no joint
two-sample VCF (70, 72, 76, 87, 103, 137, 155, 173, 189, 197, 200, 211, 219, 220, 231, 257) and therefore no intra-pair chromosomal distance.
The core-genome counts previously reported for the UKE clone 70+76 (values 3 and 52) come
from the superseded windowed BLASTN method and are **not** carried into this repository;
the historical method is documented for transparency in
`docs/legacy_SNP_calling_method_HR163.md`.

**f. Variant-effect scripts.** `esm_score.py` and `struct_analysis.py` are not present in
the project archive; only their output table (`data/vep_missense.csv`) and figure survive.
The scoring definition (ESM-2 650M, masked-marginal LLR at the site; RSA, pLDDT and
catalytic-site distance on the AlphaFold DB model) is stated in `docs/methods.md`, but the
exact script is not available.
