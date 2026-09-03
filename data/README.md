# Derived tables

Derived tables only: one row per pair, per variant or per gene. No patient identifiers,
no admission dates, no laboratory accession codes. Pairs are integers.
Licensed CC-BY-4.0 (see `../LICENSE`).

| File | Content |
|---|---|
| `HR_cohort_mechanism_from_db.csv` | pair -> centre, species, candidate mechanism |
| `pairs_vcf_availability.csv` | per-pair availability of the joint VCF / intra-pair distance |
| `HR_intrapair_chromosomal_distance.csv` | one row per pair, with audit columns |
| `HR_intrapair_differential_sites.csv` | differential sites: allele fraction, depth, class |
| `table2_snv.csv` | Table 2 of the manuscript (variants per pair) |
| `HR_betalactam_dosage_ONT_15pairs.csv` | gene dosage, R/S coverage ratios |
| `HR_clonality_distance_matrix_66genomes.csv` | mash sketch distances, 66 genomes |
| `HR_intraclone_matrices_x1e6.csv` | intra-clone submatrices (values x1e6) |
| `vep_missense.csv` | missense variants: ESM-2 LLR, RSA, pLDDT, catalytic-site distance |
| `spec163_SNP_method_comparison.csv` | pair 163: superseded windowed method vs joint calling |
| `tablesSNPeCNV_completed.xlsx` | workbook: snp, cnv, PAP_metrics, amp_units, methods, vep_missense |
