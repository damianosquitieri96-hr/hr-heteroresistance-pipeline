# Stage 00 — hybrid assembly

Run by hand. **The command line is not recoverable** (see `docs/PROVENANCE_GAPS.md`, item c).

What is known: each isolate was assembled from paired short reads plus ONT long reads into
a hybrid assembly; the per-pair susceptible assembly `<pair>S.fasta` is used as the
reference for all downstream within-pair steps. Assembler, version, parameters and any
polishing round are unknown and are not guessed here.

To complete this file, one line per stage is required (assembler + version + parameters,
then polisher + version if used).
