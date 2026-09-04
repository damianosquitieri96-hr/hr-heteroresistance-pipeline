# Stage 00 — hybrid assembly

Run by hand. **The command line is not recoverable** (see `docs/PROVENANCE_GAPS.md`, item c).

What is known: each isolate was assembled from paired short reads plus ONT long reads into
a hybrid assembly; one arm's assembly per pair is then used as the reference for all
downstream within-pair steps — the susceptible arm in 22 pairs and the resistant arm in
11, listed in `data/pair_reference_arm.csv`. Assembler, version, parameters and any
polishing round are unknown and are not guessed here.

To complete this file, one line per stage is required (assembler + version + parameters,
then polisher + version if used). Run `pipeline/provenance/recover_c_assembler.sh` against
the working directory on external storage: if it prints a command line from an assembler
log, that line closes this file. The contig headers are `>contig_1`, which is consistent
with Flye's default naming but is an indication, not an identification, and is deliberately
not written into the Methods as the assembler.
