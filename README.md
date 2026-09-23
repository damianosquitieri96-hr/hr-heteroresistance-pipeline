# Beta-lactam heteroresistance in bloodstream infections — analysis pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22288189.svg)](https://doi.org/10.5281/zenodo.22288189)

Analysis code, exact command lines, tool versions and thresholds for the two-centre
study of beta-lactam heteroresistance in clonal resistant/susceptible (R/S) isolate
pairs from bloodstream infection episodes (33 pairs, two centres).

This repository documents the pipeline **as it was actually run**. Nothing here is a
re-implementation or a reconstruction: where a command line could not be recovered from
the surviving files, this is stated explicitly in [docs/PROVENANCE_GAPS.md](docs/PROVENANCE_GAPS.md)
rather than filled with plausible parameters.

## Study design (as it bears on the analysis)

Each unit of analysis is a *pair* of isolates recovered from the same bloodstream
infection episode: one arm resistant (R) on the beta-lactam under study, one susceptible (S).
Pairs are identified by number only (40, 70, 72, 76, ...). Two centres contributed
pairs; centre, species and candidate mechanism per pair are in
`data/HR_cohort_mechanism_from_db.csv`.

All within-pair comparisons are **referenced to one arm of the same pair**: that arm's
assembly is the reference for annotation, mapping, gene-dosage estimation and joint variant
calling alike, so each pair sits in a single coordinate frame. The reference arm is the
susceptible one in 22 pairs and the resistant one in 11 (`data/pair_reference_arm.csv`);
consequently the coordinates and the `af_R`/`af_S` polarity in
`data/HR_intrapair_differential_sites.csv` are per-pair and not comparable across pairs.
No inter-patient comparison is made anywhere in the pipeline.

## How to read this repository

| Path | Content |
|---|---|
| `docs/methods.md` | full Methods text; the source the manuscript is cut from |
| `docs/PROVENANCE_GAPS.md` | parameters that could not be recovered, and what that limits |
| `pipeline/00_assembly.md` .. `09_snv_refinement/` | one file per stage; `.sh`/`.py` if scripted, `.md` if run by hand |
| `pipeline/09_snv_refinement/` | short-read refinement of the SNV candidates, 10 pairs, allele fraction per arm |
| `data/` | derived tables only — one row per pair or per variant, no patient-level data |
| `results/` | final figures |
| `submission/ena/` | pre-filled ENA/Webin metadata for the 101 read libraries, not yet submitted |

## Data availability

Raw reads and assemblies are **not** in this repository, and are **not yet deposited**:
they will be submitted to ENA prior to publication, and the BioProject accession will be
added to this section and to `docs/methods.md` at that point. Until then no accession
should be cited for this study. The metadata for that submission — 66 samples,
101 experiments and 101 runs against checklist ERC000028 — is pre-filled in
`submission/ena/`, with every value's basis and every field still to be completed listed
in `submission/ena/README.md`.
`data/` contains only derived tables (CC-BY-4.0); code is MIT.

## Coverage of the cohort

Joint two-sample variant calling — and therefore the intra-pair chromosomal distance —
is available for 17 of 33 pairs. For the remaining 16 pairs
(70, 72, 76, 87, 103, 137, 155, 173, 189, 197, 200, 211, 219, 220, 231, 257) no joint VCF exists and **no intra-pair distance is reported**.
Per-pair status: `data/pairs_vcf_availability.csv`.

## SNV refinement (stage 09)

In the 10 pairs whose candidate mechanism is an SNV (80, 202, 216, 226, 236, 240, 246,
254, 275, 279) the long-read candidates were re-tested with short reads on a
Polypolish-polished reference, reporting allele fraction per arm rather than genotype.
All 10 primary candidates are confirmed (9 fixed, 246 *dacB* N312K subclonal at 0.71);
all 10 long-read-only subclonal calls are artefacts. Method, parameters and versions:
[`pipeline/09_snv_refinement/README.md`](pipeline/09_snv_refinement/README.md);
tables: `data/snv_refinement/`.

## Citation

See `CITATION.cff`. Release v1.0.1 is archived at Zenodo: DOI [10.5281/zenodo.22288189](https://doi.org/10.5281/zenodo.22288189).
