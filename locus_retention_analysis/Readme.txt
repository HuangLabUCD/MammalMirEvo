Description

This R script (retention_analysis.R) compares genomic-locus retention between intronic and intergenic Eutherian-origin miRNA-like families. It validates the input data, calculates descriptive statistics, evaluates overdispersion, fits binomial and quasibinomial generalized linear models, performs a Wilcoxon rank-sum sensitivity analysis, and generates a visualization of the observed retention distributions.

Features

Data Preparation and Validation:
1) Imports data from the "Sheet1" worksheet of eutheria_intronic.xlsx
2) Retains loci classified unambiguously as intronic or intergenic
3) Uses intergenic loci as the reference category
4) Calculates the number of species in which each locus was not retained
5) Recalculates retention rates and checks them against the supplied values
6) Verifies that every family has 58 assessable species and that locus counts are valid

Descriptive Statistics:
1) Number of miRNA-like families in each genomic-location category
2) Mean, standard deviation, median, minimum, and maximum retained-locus counts
3) Mean, standard deviation, median, and interquartile range of retention rates
4) Mean and median retention percentages

Statistical Analysis:
1) Ordinary binomial generalized linear model for evaluating overdispersion
2) Pearson chi-square dispersion estimate
3) Quasibinomial generalized linear model comparing intronic with intergenic loci
4) Intronic-versus-intergenic odds ratio with a 95% confidence interval
5) Model-predicted retention probabilities for both genomic-location categories
6) Wilcoxon rank-sum sensitivity analysis of family-level retention rates
7) Difference between the mean intronic and intergenic retention rates

Visualization Output:
1) Boxplot and jittered observations comparing retention rates between intronic and intergenic loci
2) Retention rates displayed as percentages from 0% to 100%

Output Files:
1) eutheria_retention_intronic_intergenic.pdf
2) eutheria_retention_descriptive_statistics.csv
3) eutheria_retention_quasibinomial_results.csv
4) eutheria_retention_predicted_probabilities.csv

Dependencies:
1) R
2) readxl
3) dplyr
4) ggplot2
5) broom
6) scales

Usage:
Place eutheria_intronic.xlsx in the working directory, ensure that the required R packages are installed, and run:

Rscript retention_analysis.R

#input: eutheria_intronic.xlsx, with data in a worksheet named "Sheet1"

Required input columns:
1) location: genomic-location category; only "intronic" and "intergenic" are analyzed
2) all_synteny_species: number of species in which the locus can be assessed; expected to equal 58 for every family
3) retention_locus: number of species retaining the locus
4) retention_rate: supplied retention proportion used for validation

The script writes all output files to the current working directory.
