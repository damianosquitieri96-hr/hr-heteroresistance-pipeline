# Derived tables

Derived tables only: one row per pair, per variant or per gene. No patient identifiers,
no admission dates, no laboratory accession codes. Pairs are integers.
Licensed CC-BY-4.0 (see `../LICENSE`).

| File | Content |
|---|---|
| `HR_cohort_mechanism_from_db.csv` | pair -> centre, species, candidate mechanism |
| `pairs_vcf_availability.csv` | per-pair availability of the joint VCF / intra-pair distance |
| `pair_reference_arm.csv` | which arm's assembly is each pair's reference (S in 22, R in 11), with the evidence for it |
| `bam_PG_records.tsv` | `@PG` program records of the surviving BAMs — closes gap (b): minimap2 2.31-r1302 |
| `assembler_provenance_scan.txt` | assembler-log scan of the working drive — closes gap (c): Flye 2.9.6-b1802, long-read only |
| `species_db_provenance_scan.txt` | classifier-output scan of the working drive — negative result, gap (d) stays open |
| `joint_vcf_provenance.csv` | the 17 joint call sets: reference arm, BAMs, bcftools version, mpileup flags, from the VCF headers |
| `HR_intrapair_chromosomal_distance.csv` | one row per pair, with audit columns |
| `HR_intrapair_differential_sites.csv` | differential sites: allele fraction, depth, class |
| `table2_snv.csv` | Table 2 of the manuscript (variants per pair) |
| `HR_betalactam_dosage_ONT_15pairs.csv` | gene dosage, R/S coverage ratios |
| `HR_clonality_distance_matrix_66genomes.csv` | mash sketch distances, 66 genomes (k=21, sketch 200 000) |
| `mash_k_recovery_grid.csv` | k / sketch-size grid vs the archived matrix; identifies the parameters above |
| `ENA_sample_inventory.csv` | 66 read-set anagraphics with the source of every ENA field |
| `HR_intraclone_matrices_x1e6.csv` | intra-clone submatrices (values x1e6) |
| `vep_missense.csv` | missense variants: ESM-2 LLR, RSA, pLDDT, catalytic-site distance |
| `spec163_SNP_method_comparison.csv` | pair 163: superseded windowed method vs joint calling |
| `snv_refinement/` | stage 09: short-read allele fractions of the SNV candidates, 10 pairs (own README) |
| `tablesSNPeCNV_completed.xlsx` | workbook: snp, cnv, PAP_metrics, amp_units, methods, vep_missense |
