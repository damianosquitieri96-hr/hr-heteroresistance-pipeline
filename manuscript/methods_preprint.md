# Materials and Methods

*Compressed draft for the preprint. Full command lines, thresholds, tool versions and
per-stage provenance are in the public repository cited under Code availability; this
section states what is needed to read the results and defers the rest there.
`PRJEBXXXXXX` is a placeholder and must be replaced after deposition.*

## Isolates and design

Thirty-three clonal resistant/susceptible (R/S) isolate pairs were collected from
bloodstream infection episodes at two centres — 27 pairs in Rome (Fondazione Policlinico
Universitario Agostino Gemelli, Italy) and 6 in Hamburg (Universitätsklinikum
Hamburg-Eppendorf, Germany) — giving 66 isolates: *Escherichia coli* (n = 36),
*Pseudomonas aeruginosa* (12), *Enterobacter cloacae* (8), *Klebsiella pneumoniae* (8) and
*Klebsiella aerogenes* (2). Each pair comprises two isolates from the same episode in the
same patient, differing in phenotype to the beta-lactam under study. Pairs are pseudonymised
as integers; no patient identifiers, admission dates or laboratory accession codes are used
in any analysis or deposited file. **Every comparison in this study is made within a pair**;
no comparison is made between patients or between centres, and no transmission inference is
drawn.

Within each pair, one arm's assembly is the reference for all downstream stages, so that
both arms are analysed in one coordinate frame. That arm is the susceptible isolate in 22
pairs and the resistant isolate in 11 (per-pair table in the repository); coordinates and
reference-relative allele fractions are therefore interpretable within a pair and not across
pairs.

## Sequencing, assembly and annotation

All 66 isolates were sequenced with Oxford Nanopore long reads and 35 of them additionally
with paired-end short reads (both arms of 17 pairs plus one further arm), giving 101 read
libraries. Each isolate was assembled from its long reads alone with Flye 2.9.6
(`--nano-hq`, nominal genome size 5 Mb) from a subsampled read set; assemblies comprise 1–3
contigs. Short reads were not used for assembly or polishing — they serve as an independent
consistency check on the coverage-based dosage estimates below. Genomes were annotated with
Bakta v1.12.1 (database v6.0 light) and resistance determinants inventoried with
AMRFinderPlus. The subsampling depth, the sequencing platform and library chemistry, read
pre-processing, and the AMRFinderPlus and species-call database releases could not be
recovered from the surviving working files and are stated as unrecovered rather than
reconstructed (see Provenance below).

## Gene dosage

Long reads from both arms were mapped onto the pair's reference assembly with minimap2 2.31
(`-ax map-ont`), sorted with samtools 1.24, and coverage was aggregated per gene and in 1-kb
windows. Gene dosage is the R/S coverage ratio after
normalisation to the median chromosomal window coverage of the same alignment, making it a
within-pair, within-alignment quantity; short-read alignments, where available, served as an
independent consistency check.

## Variant calling and intra-pair distance

Variants were called jointly for the two arms in a single two-sample pileup against the
pair's reference assembly (bcftools 1.24/htslib 1.24, `mpileup -q 20 -Q 20 -a AD` →
`call -mv` → `norm`; command line recovered verbatim from the call-set headers). Joint
two-sample calling is what makes an arm-versus-arm allele-fraction comparison meaningful at
a site, and retained allelic depths are what make the classification below possible. Call
sets exist for 17 of the 33 pairs.

Distances are taken from read evidence, not from assembly-versus-assembly comparison over
sliding windows: in that approach paralogue misalignment within repeat families inflates
mismatches by an amount that shifts with the window boundaries, so it is not correctable by
a scaling factor. The superseded counts are not carried forward. Requiring ≥ 10 reads in
both arms, sites were classed as *fixed* (one arm ≥ 0.90 non-reference fraction, the other
≤ 0.10), *subclonal* (≥ 0.20 vs ≤ 0.05, differing by ≥ 0.20) or *shared* — the last
departing from the reference consensus without differing between arms, and excluded. Only
the largest contig was treated as the chromosome; runs of ≥ 5 differential sites spaced
≤ 500 bp were collapsed to one event; indels were counted separately and never added to the
substitution distance; and every discarded call is reported in audit columns, so a distance
reconciles with its raw call set. The reported distance is dispersed fixed substitutions
plus dense tracts, each tract counted once, and is 0–4 per pair (median 1) across the 17
pairs, with two pairs additionally carrying a tract. The classification uses the maximum and
minimum of the two allele fractions and is therefore symmetric in the arms. Permissive
subclonal counts are threshold-dominated — allele fractions pile up immediately above the
0.20 cut-off in depth-depleted arms — so only a filtered count (≥ 20 supporting reads,
between-arm depth ratio ≥ 0.5) is reported. The 16 pairs without a call set receive no
distance, and no value from the superseded method is substituted for them.

## Clonality, predicted variant effect and PAP

All 66 genomes were compared pairwise by k-mer sketch (Mash v2.3, `sketch -k 21 -s 200000`,
then `dist`). These parameters did not survive in any command line and were identified by
re-executing the stage over a grid of 11 k-mer lengths × 7 sketch sizes against the archived
distance matrix: k = 21 with sketch size 200 000 is the only combination that reproduces it,
to the six significant digits at which the matrix was stored. Sketch distances are reported
as such — not converted to SNP counts and not compared with published clonality thresholds,
which address inter-patient transmission rather than two isolates from one episode.

Missense variants were scored with ESM-2 650M (`esm2_t33_650M_UR50D`) as the masked-marginal
log-likelihood ratio between mutant and wild-type residue, with relative solvent
accessibility, per-residue pLDDT and distance to the catalytic site taken from AlphaFold DB
models. Population analysis profile (PAP) metrics were determined per pair and correlated
with the gene-dosage estimates above.

## Provenance

The analysis was reconstructed from surviving intermediate files rather than a
contemporaneous run log, and the repository separates parameters recovered verbatim from
tool-written headers and logs (assembler, read mapper, variant caller and their versions and
command lines, per-pair reference arm) from those recovered by re-execution against archived
output (the Mash parameters above) and those not recovered at all (read subsampling depth,
sequencing platform and library preparation, AMRFinderPlus and species-call database
releases, and the variant-effect scripts, of which only the scoring definition and output
table survive). Where a recovered parameter comes from a subset of the cohort — assembler
logs survive for 20 of 66 arms, alignment headers for 6 — the repository records the subset
and treats the remainder as an inference from their uniformity. Unrecovered parameters are
declared, not replaced by plausible defaults, and each is listed with the specific file that
would close it.

## Data availability

Raw reads for all 101 libraries from the 66 isolates will be deposited in the European
Nucleotide Archive under BioProject `PRJEBXXXXXX` (66 samples against checklist ERC000028;
66 long-read and 35 short-read runs). Sample and run metadata prepared for deposition, and
the derived tables underlying every figure and table in this manuscript, are in the
repository below, so deposited records and analysis inputs can be checked against each
other.

> **Before submission, replace `PRJEBXXXXXX` with the accession returned by Webin and add
> the sample accession range.** No ENA accession exists for this study until deposition and
> none should be cited. Two items must be settled first: per-isolate collection dates, which
> the checklist requires and which survive in no analysis file, and a genus-level species
> discordance on one pair (*E. coli* by genomic call vs *Citrobacter koseri* in the clinical
> database).

## Code availability

All analysis code, with the exact command lines, tool versions and thresholds used at every
stage, and with the per-stage record of which parameters were recovered and which could not
be, is public at <https://github.com/damianosquitieri96-hr/hr-heteroresistance-pipeline>
(release v1.0.1, archived at Zenodo, doi:10.5281/zenodo.22288189).
