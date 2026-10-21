Overview

This R script (homology_recovery.R) evaluates the performance of homology-based prediction methods for microRNA (miRNA) loci by comparing predicted copy numbers against the manually curated miRNA annotations in MirGeneDB (v3) across 19 mammalian species.

Key Features

The script performs three main analyses:

Comparative Analysis: Compares homology-based predicted miRNA-like locus counts with curated miRNA family copy numbers from MirGeneDB (v3)

Recovery Rate Calculation: Quantifies the recovery rate of homology-based predictions relative to MirGeneDB annotations for each species

Correlation Analysis: Performs Spearman's rank correlation test to assess the relationship between locus copy numbers from both methods

Important: All analyses are performed twice:

Including shared zero copy counts (where both methods predict 0 copies)
Excluding shared zero copy counts

This dual approach provides a more nuanced understanding of method performance, particularly for species or miRNA families with absent loci.

Input:
1) miRNA_copy_number_homology.csv # homology-based predicted locus copy numbers across 19 species
2) copy_number_table_19_mirgenedb # curated copy numbers from MirGeneDB (v3) across the same 19 species

Output:

The script generates:

1) Statistical summaries of recovery rates per species
2) Spearman correlation coefficients and p-values
3) Visualization plots (if implemented)
4) Tabular outputs comparing predicted vs. curated copy numbers
