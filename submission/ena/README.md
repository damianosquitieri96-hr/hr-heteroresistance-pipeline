# ENA / Webin submission templates

Pre-filled metadata for the raw reads of the 33 heteroresistant pairs. **The reads are not
yet deposited** — no PRJEB/PRJNA/ERR/SRR accession exists for this project. These templates
are the starting point of that submission; the upload itself is done from the machine that
holds the reads, and no Webin credentials are used or stored anywhere in this repository.

## What is here

| file | level | rows |
|---|---|---|
| `ena_samples.tsv` | SAMPLE | 66 — one per arm of each pair |
| `ena_experiments.tsv` | EXPERIMENT | 101 — one per library |
| `ena_runs.tsv` | RUN | 101 — one per library |
| `fill_ena_runs.sh` | — | fills read filenames and MD5s from the reads |

66 samples but 101 libraries: **every** arm has a long-read library, while a short-read
library exists for both arms of 17 pairs (70, 72, 76, 87, 103, 137, 155, 173, 189, 197,
200, 211, 219, 220, 231, 257, 275) and for arm S only of pair 40 — 35 short-read
libraries in total. This inventory is derived from which `<pair><arm>.{long,short}.flagstat.txt`
files exist in the project export, not from the reads themselves; if a short-read library
exists on storage for a sample listed here as long-read only, that row is missing and must
be added. `library_layout` is likewise evidence-based: the short-read flagstats report
equal read1/read2 counts and 96–99 % properly-paired (PAIRED), the long-read flagstats
report zero paired reads (SINGLE).

`data/ENA_sample_inventory.csv` in the repository root carries the same 66 rows with the
source of every field, plus the mechanism class, assembly availability and joint-VCF
availability per sample.

## Sample checklist

`ena_samples.tsv` targets **ERC000028**, the ENA prokaryotic pathogen minimal sample
checklist, whose field labels and mandatory/optional status were read from
`https://www.ebi.ac.uk/ena/browser/api/xml/ERC000028` rather than transcribed. Its six
mandatory fields are `isolation_source`, `collection date`,
`geographic location (country and/or sea)`, `host health state`, `host scientific name`
and `isolate`.

Filled, with the basis for each value:

- `tax_id` / `scientific_name` — NCBI Taxonomy IDs retrieved by scientific name
  (*E. coli* 562, *K. pneumoniae* 573, *K. aerogenes* 548, *E. cloacae* 550,
  *P. aeruginosa* 287, *C. koseri* 545). Species from `species.tsv` where present
  (24 pairs), otherwise from `data/HR_cohort_mechanism_from_db.csv` (9 pairs).
- `geographic location (country and/or sea)` — **inferred from the collecting centre**,
  UKE → Germany, FPG → Italy. Both are valid INSDC values in the checklist. Verify before
  submitting: the centre is where the isolate was collected, which is not by definition
  the country recorded for the sample.
- `host health state` — `diseased`, a permitted value of that controlled vocabulary,
  for bloodstream-infection isolates.
- `host scientific name` — `Homo sapiens`.
- `isolation_source` — `blood`. This is a project-level assumption from the study design
  (bloodstream infection) and is not recorded per isolate in any surviving file. Check it
  against the clinical records before submitting.
- `isolate`, `sample_alias` — `HR<pair><arm>`, the study-internal pair code already public
  in this repository. **No patient identifier appears in any column of any of these files.**

Left as `TO_FILL`:

- `collection date` — per-isolate, not present in any surviving file. ENA accepts a year
  (`2021`) or a month (`2021-03`) if the exact date cannot be released.

### One discordance to resolve first

Pair **219** is *Escherichia coli* in `species.tsv` and *Citrobacter koseri* in
`HR_cohort_mechanism_from_db.csv` — a genus-level discordance between the genomic species
call and the clinical database record. The templates carry the genomic call (562) for
HR219R and HR219S because that is the call made on these assemblies, but a genus-level
disagreement can also mean a mislabelled or swapped isolate. Resolve it before submitting;
depositing a wrong `tax_id` requires a curator to fix it afterwards.

## Experiment and run levels

Filled because they follow from the study design: `library_strategy=WGS`,
`library_source=GENOMIC`, `library_selection=RANDOM`, `library_layout`, and
`design_description`.

Deliberately empty, to be completed from the sequencing records: `platform`,
`instrument_model`, `library_name`, `nominal_length`, `nominal_sdev`,
`library_construction_protocol`. `study_accession` is `TO_FILL` until the study (project)
is registered in Webin — register it first, then paste the PRJEB accession into that
column.

Note that Webin's own interactive spreadsheet merges the experiment and run levels into a
single read-file sheet. The three-level split here matches the ENA metadata model
(SAMPLE / EXPERIMENT / RUN); when using the interactive route, join
`ena_experiments.tsv` and `ena_runs.tsv` on `experiment_alias`.

## Sequence

1. Register the study in Webin, obtain the PRJEB accession, and put it in the
   `study_accession` column of `ena_experiments.tsv`.
2. Fill `collection date` in `ena_samples.tsv`; confirm `isolation_source`, the country
   values, and the pair-219 species. Submit the samples, and keep the returned ERS
   accessions.
3. Fill `platform` and `instrument_model` in `ena_experiments.tsv` for the ONT and
   Illumina rows.
4. With the SSD mounted, fill filenames and checksums:

   ```
   ./fill_ena_runs.sh /Volumes/PortableSSD/HR_WORK/reads ena_runs.tsv ena_runs.filled.tsv
   ```

   Files are matched to a run by the `<pair><arm>` token in the basename, delimited so
   that `40R` cannot match `240R`, and split into ONT / Illumina by `_R1`/`_R2` markers or
   an `ont|nanopore|long|guppy|dorado` substring. Any run the script cannot resolve to
   exactly one file is reported on stderr and left as `TO_FILL` — it never guesses a
   filename. If the real naming convention defeats the heuristic, adjust `pick_file`.
5. Upload the FASTQ files to the Webin upload area and submit the runs.
6. Once accessions exist, replace the "reads not yet deposited" statement in `README.md`
   and `docs/methods.md` with the PRJEB accession.
