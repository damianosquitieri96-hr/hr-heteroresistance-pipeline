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

| Item | Script | Status |
|---|---|---|
| b | `recover_b_mapper.sh` | **closed** — minimap2 2.31-r1302 from `@PG`, 6 BAMs of 3 pairs |
| c | `recover_c_assembler.sh` | **closed** — Flye 2.9.6-b1802 from 20 logs of 10 pairs; corrected "hybrid" |
| c2 | — | open — subsampling tool/depth behind Flye's `sub.fq.gz` input |
| d | `recover_d_species_db.sh` | open — `HR_WORK` scanned, no classifier output found; try `HR_RESULTS` |
| e | `recover_e_jointvcf.sh` | open by choice — would extend coverage, corrects nothing |

**b. Long-read mapper and version.** — **CLOSED 2026-09-04.** Recovered verbatim from the
`@PG` records of the BAM headers: **minimap2 2.31-r1302**,
`minimap2 -t 5 -I 1G -K 100M -ax map-ont <ref>.fasta <reads>.fastq.gz`, piped to
**samtools 1.24** `sort -@ 4 -m 512M`. `--threads` was 5 or 6 depending on the run; nothing
else varied. `pipeline/02_mapping_dosage.sh` now carries the real command instead of a
placeholder. Records: `data/bam_PG_records.tsv` (18 `@PG` rows).

Evidence scope, stated because it is narrower than the cohort: only 6 BAMs, of 3 pairs
(103, 137, 155), are still on the drive — the alignments of the other 30 pairs were deleted
after use. All 6 agree in mapper, version and parameters, and in each recovered command line
the reference file is that pair's reference arm as listed in `data/pair_reference_arm.csv`
(3/3), which is an independent check that these BAMs are the pipeline's own. That the other
60 alignments used the same mapper is an inference from that agreement, not a recovered
fact.

**c. Assembler, version and parameters.** — **CLOSED 2026-09-04**, and it corrected a
substantive claim (see Corrections below). Recovered from the Flye logs on external storage:
**Flye 2.9.6-b1802**, `flye --nano-hq sub.fq.gz --out-dir <asm>/flye --threads 5
--genome-size 5m`, preset config `asm_nano_hq.cfg`; `--threads` 5 or 8 depending on the run.
`--min-ovlp 10000` appears in the internal `flye-modules` calls and comes from the preset,
not the command line. The same `--genome-size 5m` was passed for every isolate irrespective
of species, including the two *P. aeruginosa* pairs; under `--nano-hq` that value feeds only
Flye's internal coverage estimate. `pipeline/00_assembly.md` now carries the command line.
Scan: `data/assembler_provenance_scan.txt`.

Evidence scope: logs survive for 20 arms of 10 pairs (103, 155, 173, 189, 197, 200, 202,
211, 212, 216) spanning four of the five species, and the command line is identical across
all of them. The remaining 46 arms are an inference from that uniformity plus the uniform
`contig_N` naming of all 32 exported assemblies.

**c2. Read subsampling (new, opened by closing c).** Flye's input is `sub.fq.gz`, a
subsampled read set produced by a step upstream of Flye, so neither the tool nor the target
depth or seed appears in the assembler log. This matters more than it looks: subsampling
depth bounds the sensitivity of everything downstream that rests on read support. Any log or
wrapper script beside `HR_WORK/<pair>/.asm_<arm>/sub.fq.gz`, or the read count of that file
against the full library, would close it.

**d. Species-call tool and database release.** — **still open, and now with a negative
result.** `recover_d_species_db.sh` was run against `/Volumes/PortableSSD/HR_WORK` on
2026-09-04 and found **no classification output of any kind** and no line declaring a
database path or release (`data/species_db_provenance_scan.txt`). So `species.tsv` — three
bare columns, no header, no tool signature — has no surviving provenance in the assembly
working directory. Two things are still worth trying before declaring it unrecoverable: run
the same script against `HR_RESULTS` and against any analysis directory outside these two,
and check the shell history of the machine that produced `species.tsv`. Related: the
genus-level discordance on **pair 219** (*E. coli* in `species.tsv`, *C. koseri* in the
clinical database) is unresolved and blocks the ENA deposition of that sample — see
`submission/ena/README.md`. The species names deposited to ENA come from `species.tsv`, so
this item is a deposition dependency, not only a documentation one.

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

**"Hybrid assembly" — corrected 2026-09-04.** This repository previously stated, in
`README.md`, `docs/methods.md`, `pipeline/00_assembly.md`, `manuscript/methods_preprint.md`
and the `design_description` field of `submission/ena/ena_experiments.tsv`, that short and
long reads were combined into a hybrid assembly per isolate. **The assemblies are long-read
only.** Four independent lines of evidence, from closing item (c):

1. the recovered Flye command line takes one ONT read file and no short reads;
2. 7 of the 10 pairs whose logs survive (103, 155, 173, 189, 197, 200, 211) *have*
   short-read libraries that were not used at this stage — so this is not a case of short
   reads being unavailable;
3. no assembler other than Flye and no polisher (medaka, racon, pilon, polypolish,
   nextPolish) leaves a log, run-parameter file or version string anywhere on the drive;
4. all 61 contigs across the 32 exported assemblies carry Flye's `contig_N` naming, with no
   Unicycler or SPAdes header among them.

Point 3 is absence of evidence, so the accurate claim is "no polishing stage is recorded",
not "no polishing was performed". Points 1, 2 and 4 are positive evidence and are what the
correction rests on.

**No result changes**, because no analysis in this repository takes the assembly method as
an input: the assemblies themselves are the input and they are the same files as before.
What changes is what the reads were used for, and there the ENA `design_description` was
wrong in a way worth naming: it claimed the short-read libraries were used for assembly
polishing *and for joint variant calling*. Neither is true. The 17 joint call sets were
made from the long-read BAMs — `data/joint_vcf_provenance.csv`, `bams` column, 0 of 17
name a short-read alignment — and the short-read libraries served as an independent
consistency check on the coverage-based dosage estimates. Both design descriptions have
been rewritten; had they gone to ENA unchanged, the deposited records would have
misdescribed 101 libraries.

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
