# Stage 00 — long-read assembly

Run by hand, per arm. Command line recovered 2026-09-04 from the Flye logs on external
storage (`data/assembler_provenance_scan.txt`):

```bash
flye --nano-hq sub.fq.gz --out-dir <asm_dir>/flye --threads 5 --genome-size 5m
```

Flye 2.9.6-b1802, preset config `asm_nano_hq.cfg`. `--threads` was 5 or 8 depending on the
run; nothing else varied. `--genome-size 5m` was passed for every isolate regardless of
species, including the two *Pseudomonas aeruginosa* pairs whose genomes are larger — under
`--nano-hq` this value only feeds Flye's internal coverage estimate and does not set the
assembly length. `--min-ovlp 10000` appears in the internal `flye-modules` invocations; it
comes from the preset config, not from the command line.

**The assemblies are long-read only, not hybrid.** This corrects the earlier statement in
this file. Evidence: the recovered command line takes a single ONT read file and no short
reads; 7 of the 10 pairs whose logs survive *have* short-read libraries that were not used
here; no assembler other than Flye and no polisher (medaka, racon, pilon, polypolish) has
a log, signature file or version string anywhere on the drive; and all 61 contigs across
the 32 exported assemblies carry Flye's `contig_N` naming, with no Unicycler or SPAdes
header among them. Short reads were used for mapping and dosage cross-checks, not for
assembly. Absence of a polisher log is weaker evidence than the command line itself, so
"no polishing step is recorded" is the accurate claim rather than "no polishing was done".

Logs survive for 20 arms of 10 pairs (103, 155, 173, 189, 197, 200, 202, 211, 212, 216),
spanning four of the five species; the command line is identical across them. That the
remaining 46 arms were assembled the same way is an inference from that uniformity plus the
contig naming, not a recovered fact.

**Still open (`docs/PROVENANCE_GAPS.md`, item c2):** the input `sub.fq.gz` is a subsampled
read set, and the subsampling step is upstream of Flye, so neither the tool nor the target
depth appears in its log. Any file named `sub.fq.gz` or a log beside it in
`HR_WORK/<pair>/.asm_<arm>/` would close this.
