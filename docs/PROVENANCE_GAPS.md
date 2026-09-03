# Provenance gaps

Parameters of the pipeline as run that could **not** be recovered from the surviving
files. They are listed here rather than reconstructed, and the Methods text points to
this file at each affected stage.

**a. Mash sketch parameters (k, sketch size).** The intermediate `mash_dist.tsv` no
longer exists; the 66-genome distance matrix survives without its command line.
Consequence: clonality is reported as sketch distance only, never as SNP counts and
never against published clonality thresholds. Recoverable by re-running
`mash sketch` + `mash dist` on the 66 assemblies and checking that the distances
reproduce `data/HR_clonality_distance_matrix_66genomes.csv`; the k that reproduces it
is then the k that was used.

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
