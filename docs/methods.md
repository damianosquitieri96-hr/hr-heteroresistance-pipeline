# Methods

## Isolates and design

33 clonal R/S isolate pairs from bloodstream infection episodes, two centres.
Pairs are pseudonymised as integers. Within each pair, one arm's assembly serves as the
reference genome for every within-pair comparison (which arm: see below).

## Hybrid assembly

Short-read and long-read (ONT) data were combined into a hybrid assembly per isolate.
Assembler, version, parameters and any polishing stage **could not be recovered** from
the surviving working files (see `docs/PROVENANCE_GAPS.md`, item c). Assembly FASTA
files carry generic `contig_N` headers, which do not identify the assembler.
Each pair has one arm whose assembly serves as the pair's reference for every downstream
stage — annotation, mapping, gene dosage and joint calling all use the same arm, so the
whole pair sits in one coordinate frame. That arm is the susceptible one in 22 pairs and
the resistant one in 11 (76, 87, 212, 226, 236, 240, 254, 265, 275, 279, 283); the full
table is `data/pair_reference_arm.csv`. It was recovered from the `_on_<pair><arm>`
infix of the coverage filenames and, for the 17 pairs with a joint call set, corroborated
by the `##reference` header of the VCF, which agreed in 17/17. The reference arm is a
property of how each pair was processed, not a study-design variable.

## Annotation

Bakta v1.12.1 with database v6.0 (light) (Schwengers et al., doi:10.1099/mgen.0.000685).
Command: `pipeline/01_annotation.sh`. Outputs used downstream: `<pair><ref>.bakta.tsv`
(gene table) and `<pair><ref>.genes.bed` (gene intervals for coverage aggregation).
AMR determinants were called with AMRFinderPlus (`<pair><ref>.amrfinder.tsv`); the version
string is not recoverable from the exports.

## Long-read mapping and gene dosage

Long reads from both arms of a pair were mapped onto the **reference arm's assembly of the
same pair** (`data/pair_reference_arm.csv`), so dosage ratios are within-pair throughout.
Per-gene coverage (`<pair>{R,S}.long_on_<pair><ref>.genecov.tsv`) and 1-kb window coverage
(`...win1kb.tsv.gz`) were computed on those alignments; short-read alignments were
processed identically (`...short_on_...`) and used as a consistency check.
Gene dosage is expressed as the R/S coverage ratio after normalisation to the
median chromosomal window coverage of the same alignment
(`data/HR_betalactam_dosage_ONT_15pairs.csv`). The mapper and its version are recorded
in the `@PG` line of the BAM headers, which reside on external storage and were not
available when this repository was assembled (`docs/PROVENANCE_GAPS.md`, item b).

## Joint two-sample variant calling

bcftools 1.24 with htslib 1.24, run as a single two-sample pileup of the R and S
alignments against the assembly of **one arm** of the same pair:

```
bcftools mpileup -q 20 -Q 20 -a AD -f <pair><ref>.fasta <pair>R.bam <pair>S.bam | bcftools call -mv -Ou | bcftools norm -f <pair><ref>.fasta -Oz -o <pair>.RS.vcf.gz
```

Recovered verbatim from the `##bcftoolsCommand`, `##bcftools_callCommand` and
`##bcftools_normCommand` header lines of the exported VCFs. `-a AD` retains allelic
depths, which the downstream distance recount uses for allele-fraction filtering.

`<ref>` is **not the same arm in every pair**: the `##reference` header lines give the
susceptible arm for 8 pairs (40, 80, 202, 216, 222, 246, 260, 298) and the resistant arm
for 9 (212, 226, 236, 240, 254, 265, 275, 279, 283). The reference arm of each pair is
tabulated in `data/joint_vcf_provenance.csv`, read from the VCF headers. Two consequences
are load-bearing for reading the outputs: the coordinates in
`data/HR_intrapair_differential_sites.csv` are relative to that pair's own reference arm
and are not comparable across pairs, and the callable fraction of each genome is defined by
the arm that happens to be the reference. The recount is symmetric in the two samples and
does not assume which arm is the reference; nothing in it required the susceptible arm.

## Intra-pair chromosomal distance — recomputed from read-based joint calls

### What changed and why

Intra-pair distances previously obtained by assembly-versus-assembly comparison over
sliding windows are not used. In that approach mismatches contributed by paralogue
misalignment inside repeat families depend on where the window boundaries fall, so the
overestimate is not constant across pairs and cannot be corrected by a scaling factor.
The distance is therefore taken from read evidence.

No new caller was installed: the exports already contain, for 17 of the 33 pairs, a
joint two-sample call set (bcftools, samples R and S) made against **one arm's own
assembly from the same pair** — the susceptible arm in 8 pairs, the resistant arm in 9
(`data/joint_vcf_provenance.csv`). Either choice satisfies the recount, which needs a
within-pair reference so that the two arms are compared on the same coordinate frame; it
does not need a particular arm.

### Procedure

Input: `<pair>.RS.vcf.gz` (samples R, S; reference `<pair><ref>.fasta`, the arm recorded in
`data/joint_vcf_provenance.csv`) and `<pair><ref>.contigs.tsv`. Script: `recount.py`.

1. **Reference.** The reference arm of the same pair (`data/pair_reference_arm.csv`).
   Coordinates and REF alleles are on that assembly; a difference is reported when the
   two arms disagree, irrespective of which arm carries the non-reference allele, so the
   procedure is symmetric and does not depend on which arm is the reference. Note that
   `af_R` and `af_S` in the output are therefore polarised by that pair's reference arm.
2. **No masking.** Neither prophage regions (PhiSpy) nor recombination (Gubbins) are
   masked. Recombination inference requires a phylogeny, which two isolates of one
   isogenic pair cannot provide; prophage masking would remove real sites without an
   independent reason to distrust them. This is declared, not applied.
3. **Chromosome only.** The largest reference contig is taken as the chromosome; any
   further contig ≥ 1 Mb is named in `extra_large_contigs` and reported separately, never
   pooled into the count. Variants on all other contigs are counted in
   `off_chromosome` and excluded.
4. **Difference classes** from per-sample allelic depths (FORMAT/AD), requiring ≥ 10
   reads in both arms:
   - *fixed*: one arm ≥ 0.90 non-reference fraction, the other ≤ 0.10;
   - *subclonal*: one arm ≥ 0.20, the other ≤ 0.05, difference ≥ 0.20;
   - *shared*: everything else — both arms carry the same allele mixture, so the site
     is a difference from the reference consensus, not between the arms. Excluded, and
     counted in `shared_sites`.
5. **Dense divergent tracts.** Runs of ≥ 5 differential sites with ≤ 500 bp spacing are
   collapsed to one event (`dense_tracts`, with `sites_in_tracts` giving their extent)
   and excluded from the dispersed counts.
6. **Indels** are counted separately (`indel_differences`) and never added to the
   substitution distance.
7. **Audit columns.** `shared_sites`, `low_depth_sites`, `off_chromosome` account for
   every call discarded, so the reported distance can be reconciled with the raw call
   set.

Reported distance = `fixed_dispersed + dense_tracts` (each tract as a single event).

## Result

Across the 17 recomputable pairs, fixed dispersed substitutions are 0–4 per pair
(median 1); two pairs additionally carry a dense tract (212: two tracts spanning 255
sites; 298: one tract of 7 sites). The fixed sites independently recover the drivers
reported in the variant table, at the expected coordinates — e.g. pair 40
`contig_1:2,810,478` inside *ampD* and pair 80 `contig_1:2,995,285` inside *ampR*.

The permissive `subclonal_dispersed` column is threshold-dominated: for the pairs with
large values the allele fractions pile up just above the 0.20 cut and the variant-carrying
arm is depth-depleted relative to the other, the signature of mapping noise.
`subclonal_filtered` therefore adds ≥ 20 reads supporting the site and a depth ratio
≥ 0.5 between the two arms; only that column should be quoted.

## Pairs that cannot be recomputed

No joint call set is present in the exports for pairs 70, 72, 76, 87, 103, 137, 155,
173, 189, 197, 200, 211, 219, 220, 231, 257. Any distance previously quoted for these —
including the core-genome counts of the UKE 70+76 clone — rests on the earlier method
and should not be reported as a SNP distance until read-based calls are available.

## Flag raised by the recount

Pair 216, the two adjacent variants recorded in the variant table as Q171R (A>G,
af_R 0.53) and Q171H (G>T, af_R 0.50) in an RND efflux membrane-fusion protein, sit at
`contig_3:888,890-888,891`. The recount classes both as *shared*, not as differences,
because the susceptible arm also carries the alternate alleles (af_S 0.077 and 0.071).
Three further observations argue against two convergent subclonal substitutions:
local depth is 17 and 20 reads against a genome-wide median of 47 for this pair; all
alternate reads fall on one strand (DP4 `22,22,0,12` and `24,25,0,13`); and the same
adjacent A>G / G>T motif recurs at `contig_3:2,076,545-2,076,546`, inside a second
RND-family transporter, again at depleted depth (8 and 7 reads). The two rows should be
treated as candidate mapping artefacts pending targeted confirmation.

## Limitations

The fraction of the reference not interrogated at ≥ 10 reads in both arms is not
quantified here; only the number of discarded low-depth call positions is reported.
No canonical clonality threshold exists for these species, and the published thresholds
address inter-patient transmission rather than an intra-patient pair sampled from one
episode, so the distances are reported as such and not compared with a cut-off.


## k-mer sketch clonality

Mash sketches of all 66 genomes (both arms of all pairs) were compared pairwise to give
`data/HR_clonality_distance_matrix_66genomes.csv`; intra-clone submatrices (values x1e6)
are in `data/HR_intraclone_matrices_x1e6.csv`, and Figure 3 (dendrograms plus intra-clone
heat maps) is produced by `pipeline/08_figures/hr_clonality_fig.py`.

The sketch parameters were **recovered by re-running the stage against its own output**
rather than from a surviving command line: `mash sketch -k 21 -s 200000` (Mash v2.3,
default seed), followed by `mash dist`. `mash sketch` and `mash dist` were re-run over a
grid of 11 k-mer lengths (13–31) × 7 sketch sizes (1 000–200 000), and each resulting
matrix was compared with the archived one. Only k = 21 with a sketch size of 200 000
reproduces it: 490 of the 496 pairwise distances available for the test are identical to
the 6 significant digits at which the matrix was stored, and the 6 remaining residuals are
≤ 4.9e-9 and all fall on intra-clone distances near 1e-4, i.e. they are rounding at the
last stored digit. No other combination in the grid produces a single exact distance; the
next-best median residual, at k = 21 and sketch size 100 000, is 4.7e-4, which is the size
of the sampling noise expected between two independent 100 000-hash sketches. The grid is
in `data/mash_k_recovery_grid.csv` and the comparison in `results/mash_k_recovery.png`.

The test used the 32 of 66 assemblies available on local storage — 496 of the 2 145
pairwise distances — because the remaining assemblies are on external storage. This is not
a limitation of the identification: `mash dist` is pairwise and a sketch does not depend on
which other genomes are sketched alongside it, so a combination that reproduces 496
distances exactly reproduces the matrix.

Sketch distances are nonetheless still reported as *sketch distances only*. They are not
converted into SNP counts and are not compared with published clonality thresholds: those
thresholds address inter-patient transmission, not two isolates from the same episode.

Figure 3 was deliberately **not** recomputed with the corrected distance method: it rests
on sketch distances, which do not pass through a windowed assembly-versus-assembly
alignment and are therefore not affected by the inflation that motivated the recount
(see `docs/legacy_SNP_calling_method_HR163.md` and
`data/spec163_SNP_method_comparison.csv`).

## Predicted variant effect

Missense variants were scored with ESM-2 650M (`esm2_t33_650M_UR50D`) as the
log-likelihood ratio between mutant and wild-type residue at the masked site
(masked-marginal scoring). Structural context was taken from AlphaFold DB models of the
same proteins: relative solvent accessibility (RSA), per-residue pLDDT and Euclidean
distance from the catalytic site were computed on that model. Scores and structural
metrics are in `data/vep_missense.csv`; the figure is `results/HR_missense_impact.png`.
The scoring scripts themselves are **not present** in the project archive
(`docs/PROVENANCE_GAPS.md`, item f).

## PAP and gene dosage

Population analysis profile (PAP) metrics per pair are in the `PAP_metrics` sheet of
`data/tablesSNPeCNV_completed.xlsx` and were correlated with the gene-dosage estimates
above; see `pipeline/07_pap_correlation/`.

## Code availability

All analysis code, together with the exact command lines, tool versions and thresholds
used at every stage — hybrid assembly, annotation, long-read mapping and gene-dosage
estimation, joint two-sample variant calling, intra-pair chromosomal distance, k-mer
sketch clonality analysis, predicted variant effect and PAP correlation — is available at
`https://github.com/damianosquitieri96-hr/hr-heteroresistance-pipeline` (release v1.0.1, archived at
Zenodo, DOI [10.5281/zenodo.22288189](https://doi.org/10.5281/zenodo.22288189)). Derived tables underlying every figure and table are included in
the same repository. Sequencing reads and assemblies will be deposited in ENA prior to publication; the
BioProject accession will be inserted here on submission.
