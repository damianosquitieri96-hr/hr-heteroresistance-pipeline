# Provenance recovery scripts

One script per open item in `docs/PROVENANCE_GAPS.md`. They exist because the parameters of
the pipeline as run survive only inside files on the external drive, which is not part of
this repository. Each script reads that drive and reports what it finds; none of them
infers a parameter from something that merely looks consistent with it.

| Script | Item | Reads | Writes |
|---|---|---|---|
| `recover_b_mapper.sh` | b — long-read mapper + version | `HR_RESULTS/*/*.bam` headers | `data/bam_PG_records.tsv` |
| `recover_c_assembler.sh` | c — assembler, version, parameters | `HR_WORK` logs | `data/assembler_provenance_scan.txt` |
| `recover_d_species_db.sh` | d — species-call tool + DB release | `HR_WORK` classifier outputs | `data/species_db_provenance_scan.txt` |
| `recover_e_jointvcf.sh` | e — the 16 pairs with no joint VCF | `HR_WORK` + `HR_RESULTS` | VCFs, into a directory you name |

## Running them

The drive layout is known exactly, from the `##bcftoolsCommand` and `##reference` headers
of the exported VCFs:

```
/Volumes/PortableSSD/HR_WORK/<pair>/<pair>{R,S}.fasta      assemblies, annotation, logs
/Volumes/PortableSSD/HR_RESULTS/<pair>/<pair>{R,S}.bam     alignments
```

With the drive mounted:

```bash
cd pipeline/provenance
./recover_b_mapper.sh                 # defaults to /Volumes/PortableSSD/HR_RESULTS
./recover_c_assembler.sh              # defaults to /Volumes/PortableSSD/HR_WORK
./recover_d_species_db.sh
./recover_e_jointvcf.sh --dry-run     # reports; writes nothing
```

`recover_b_mapper.sh` needs `samtools`; `recover_e_jointvcf.sh --run` needs `bcftools`
(the existing call sets are 1.24 + htslib 1.24, and it warns if yours differs). Pass a
different root as the first argument if the drive is mounted elsewhere.

## Rules these scripts follow

- **Read-only on the drive.** Only `recover_e_jointvcf.sh --run` writes, and only into an
  output directory passed on the command line, never onto the drive. Nothing is deleted.
- **A hint is reported as a hint.** `recover_c_assembler.sh` prints the contig naming and
  says which assembler it is consistent with, and separately says that this is not an
  identification. Only a command line recovered from a log closes item c.
- **Silence is an answer.** If a script finds nothing, the item stays open in
  `docs/PROVENANCE_GAPS.md`. Do not fill it with a plausible parameter: a wrong command line
  in the Methods is worse than a declared gap.

## After a script succeeds

| Item | Where the recovered value goes |
|---|---|
| b | the mapper invocation in `pipeline/02_mapping_dosage.sh`, and the Methods mapping paragraph |
| c | `pipeline/00_assembly.md`, one line per stage (assembler, then polisher) |
| d | the species-call paragraph in `docs/methods.md`, as tool + version + database release |
| e | re-run stage 05 over all 33 pairs; the distance is recomputed for all, not concatenated |

In every case, move the item out of `docs/PROVENANCE_GAPS.md` and state in the commit
message which file the value came from, so the recovery is auditable.
