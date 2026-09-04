# Materials and Methods

*Draft for the preprint. Every command line, threshold and tool version stated here is the
one in the repository cited under Code availability; parameters that could not be recovered
from the surviving files are stated as unrecovered rather than reconstructed. Placeholders
in the form `PRJEBXXXXXX` must be replaced with real accessions after deposition, and no
accession may be cited until then (see the note under Data availability).*

## Isolates and study design

Thirty-three clonal resistant/susceptible (R/S) isolate pairs were collected from
bloodstream infection episodes at two centres, 27 pairs at Fondazione Policlinico
Universitario Agostino Gemelli (Rome, Italy) and 6 at Universitätsklinikum
Hamburg-Eppendorf (Hamburg, Germany), giving 66 isolates. Each pair comprises two isolates
recovered from the same episode in the same patient, one expressing the resistant and one
the susceptible phenotype to the beta-lactam under study; pairs are pseudonymised as
integers throughout, and no patient identifiers, admission dates or laboratory accession
codes are used in any analysis or deposited file. The species distribution of the 66
isolates is *Escherichia coli* (n = 36), *Pseudomonas aeruginosa* (12), *Enterobacter
cloacae* (8), *Klebsiella pneumoniae* (8) and *Klebsiella aerogenes* (2). All comparisons in
this study are made **within** a pair; no comparison is made between patients or between
centres, and no transmission inference is drawn.

Within each pair, the assembly of one arm serves as the reference genome for every
downstream stage — annotation, read mapping, gene-dosage estimation and variant calling
alike — so that both arms of a pair are analysed in a single coordinate frame. That
reference arm is the susceptible isolate in 22 pairs and the resistant isolate in 11; the
per-pair assignment is tabulated in the repository. Because the reference arm differs
between pairs, genomic coordinates and reference-relative allele fractions are interpretable
within a pair but are not comparable across pairs.

## Sequencing and assembly

All 66 isolates were sequenced with Oxford Nanopore long reads; 35 of them were
additionally sequenced with paired-end short reads (both arms of 17 pairs, plus the
susceptible arm of one further pair), giving 101 read libraries in total. Short-read and
long-read data were combined into one hybrid assembly per isolate. The assembler, its
version, its parameters and any polishing stage **could not be recovered** from the
surviving working files: the assemblies carry generic `contig_N` headers and no assembler
log survives in the project archive. This gap is declared rather than reconstructed, and the
repository documents the recovery procedure that would close it. For the same reason,
sequencing platform and instrument model, library preparation chemistry and read
pre-processing (adapter trimming, length or quality filtering) are not stated here: they
were not recoverable from the exported files and are left unasserted rather than inferred.

Genomes were annotated with Bakta v1.12.1 against database v6.0 (light). Antimicrobial
resistance determinants were inventoried with AMRFinderPlus; the version string and database
release as run are not recoverable from the exports and are therefore not stated.

## Read mapping and gene dosage

Long reads from both arms of each pair were mapped onto that pair's reference-arm assembly.
Per-gene coverage and coverage in 1-kb windows were computed on the resulting alignments;
where short reads were available, they were processed identically and used as an independent
consistency check on the dosage estimates. Gene dosage is expressed as the
resistant/susceptible coverage ratio after normalisation to the median chromosomal window
coverage of the same alignment, so that a ratio is a within-pair, within-alignment quantity.
The mapper and its version are recorded in the `@PG` line of the alignment headers, which
reside on external storage that was not available when the analysis repository was
assembled; this is declared as an open item.

## Joint two-sample variant calling

Variants were called jointly for the two arms of a pair in a single two-sample pileup
against that pair's reference-arm assembly, using bcftools 1.24 with htslib 1.24:

```
bcftools mpileup -q 20 -Q 20 -a AD -f <pair><ref>.fasta <pair>R.bam <pair>S.bam \
  | bcftools call -mv -Ou \
  | bcftools norm -f <pair><ref>.fasta -Oz -o <pair>.RS.vcf.gz
```

This command line was recovered verbatim from the `##bcftoolsCommand`,
`##bcftools_callCommand` and `##bcftools_normCommand` headers of the call sets themselves.
Joint two-sample calling, rather than two independent single-sample calls, is what makes an
arm-versus-arm allele-fraction comparison meaningful at a given site; retaining allelic
depths (`-a AD`) is what makes the subclonal classification below possible. Joint call sets
exist for 17 of the 33 pairs.

## Intra-pair chromosomal distance

Intra-pair distances were computed from read evidence rather than by
assembly-versus-assembly comparison. In the windowed alignment approach used previously,
mismatches contributed by paralogue misalignment within repeat families shift with the
window boundaries, so the resulting overestimate is not constant across pairs and cannot be
removed by a scaling factor; the earlier counts are therefore not carried forward, and the
superseded method is documented separately in the repository for transparency.

From each joint call set, sites were classified using per-sample allelic depths, requiring at
least 10 reads in both arms: *fixed* when one arm carries the non-reference allele at a
fraction ≥ 0.90 and the other ≤ 0.10; *subclonal* when one arm is ≥ 0.20, the other ≤ 0.05
and the two differ by ≥ 0.20; and *shared* otherwise — a shared site departs from the
reference consensus but does not differ between the arms, and is excluded. Only the largest
reference contig was treated as the chromosome; any additional contig ≥ 1 Mb is reported
separately and never pooled, and variants on all other contigs are counted and excluded.
Runs of ≥ 5 differential sites spaced ≤ 500 bp apart were collapsed to a single event, so
that one divergent tract does not enter the distance as tens of independent substitutions.
Indels were counted separately and never added to the substitution distance. Every discarded
call is accounted for in audit columns (shared, low-depth, off-chromosome), so a reported
distance can be reconciled with the raw call set. The reported distance is the number of
dispersed fixed substitutions plus the number of dense tracts, each tract counted once.

Across the 17 pairs with a joint call set, dispersed fixed substitutions number 0–4 per pair
(median 1), and two pairs additionally carry a dense tract. Because the classification rests
on the maximum and minimum of the two allele fractions, it is symmetric in the two arms and
does not depend on which arm is the pair's reference. The permissive subclonal count is
threshold-dominated and is not quoted: in the pairs with large values the allele fractions
accumulate immediately above the 0.20 cut-off and the variant-carrying arm is depth-depleted
relative to the other, the expected signature of mapping noise. A filtered subclonal count,
additionally requiring ≥ 20 supporting reads and a between-arm depth ratio ≥ 0.5, is the
only subclonal figure reported.

Sixteen pairs have no joint call set and therefore no read-based intra-pair distance; no
distance is reported for them, and no value from the superseded method is substituted. These
pairs contribute to the depth-based dosage and sketch-based clonality analyses, which use no
variant calls.

## Population-level clonality

Mash sketches of all 66 genomes were compared pairwise with `mash dist` (Mash v2.3,
`mash sketch -k 21 -s 200000`, default seed). The sketch parameters were not available from
a surviving command line and were instead **identified by re-executing the stage against its
own archived output**: sketching and pairwise distance computation were repeated over a grid
of 11 k-mer lengths (13–31) and 7 sketch sizes (1 000–200 000), and each resulting matrix was
compared with the archived one. Only k = 21 with a sketch size of 200 000 reproduces it —
490 of the 496 distances available for the test are identical to the six significant digits
at which the matrix was stored, and the six residual differences are ≤ 4.9 × 10⁻⁹, all on
intra-clone distances near 10⁻⁴, i.e. rounding at the last stored digit. No other
combination in the grid yields a single exact distance, and the next-best combination's
median residual (4.7 × 10⁻⁴, at sketch size 100 000) is the magnitude of the sampling noise
expected between two independent sketches of that size. The test used the 32 assemblies
available locally, covering 496 of the 2 145 pairwise distances; since `mash dist` is
pairwise and a sketch does not depend on which other genomes are sketched alongside it, this
does not weaken the identification.

Sketch distances are reported as sketch distances. They are not converted into SNP counts
and are not compared with published clonality thresholds, which address inter-patient
transmission rather than two isolates recovered from a single episode.

## Predicted variant effect

Missense variants were scored with ESM-2 650M (`esm2_t33_650M_UR50D`) as the log-likelihood
ratio between the mutant and wild-type residue at the masked site (masked-marginal scoring).
Structural context was taken from AlphaFold DB models of the same proteins: relative solvent
accessibility, per-residue pLDDT and Euclidean distance to the catalytic site were computed
on those models. The scoring scripts are not present in the project archive; the scoring
definition above and the resulting table are, and the absence of the scripts is declared as
an open item.

## Population analysis profiles

Population analysis profile (PAP) metrics were determined per pair and correlated with the
gene-dosage estimates described above; the metrics and the correlation analysis are included
in the repository.

## Provenance and reproducibility statement

The analysis was reconstructed from surviving intermediate files rather than from a
contemporaneous run log. Accordingly, the repository distinguishes three categories, and
this Methods section follows the same discipline: parameters recovered **verbatim** from
file headers (the variant-calling command line and caller version, the reference arm of each
pair); parameters recovered by **re-execution against archived output** (the Mash k-mer
length and sketch size, identified as the unique combination that reproduces the stored
distance matrix); and parameters **not recovered** (assembler and polishing, read mapper and
version, AMRFinderPlus and species-call database releases, sequencing platform and library
preparation, the variant-effect scripts). Unrecovered parameters are stated as such here and
are not replaced by plausible defaults. The repository lists each open item together with the
specific file that would close it.

## Data availability

Raw sequencing reads for all 101 libraries from the 66 isolates will be deposited in the
European Nucleotide Archive under BioProject accession `PRJEBXXXXXX`, with 66 samples
registered against the ENA prokaryotic pathogen minimal sample checklist (ERC000028) and 101
experiment/run pairs (66 long-read, 35 short-read). Sample and run metadata prepared for
deposition — per-isolate species assignment with NCBI taxonomy identifiers, collecting
centre, country and library layout — are included in the analysis repository, so that the
deposited records and the analysis inputs can be checked against one another. Derived tables
underlying every figure and table in this manuscript are included in the same repository.

> **Before submission, replace `PRJEBXXXXXX` with the BioProject accession returned by Webin
> and add the sample accession range.** Until deposition is complete no ENA accession exists
> for this study and none should be cited. Two items must be settled first: the per-isolate
> collection dates, which the checklist requires and which do not survive in any analysis
> file, and a genus-level species discordance on one pair (*Escherichia coli* by genomic call
> versus *Citrobacter koseri* in the clinical database), which must be resolved before that
> sample is registered.

## Code availability

All analysis code, together with the command lines, tool versions and thresholds used at
every stage, is available at
<https://github.com/damianosquitieri96-hr/hr-heteroresistance-pipeline> (release v1.0.1,
archived at Zenodo, doi:10.5281/zenodo.22288189). The repository also documents, per stage,
which parameters were recovered from surviving files and which could not be, and provides
the scripts that would close each remaining gap.
