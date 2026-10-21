Description

This R script (miRNA_locus_MRCA_inference.R) identifies the most recent common ancestor (MRCA) node of each synteny-verified miRNA family. For each family, it identifies all species with a copy number of at least one and reports their MRCA node in a rooted phylogenetic tree.

The method assumes a single ancestral origin followed by possible lineage-specific losses. It reports the MRCA node only and does not estimate evolutionary age.

Features

Input Validation:
1) Verifies that the tree file and Excel workbook exist
2) Requires a rooted tree with unique tip labels
3) Checks that the expected Excel worksheet and miRNA_family column are present
4) Checks for missing or duplicated miRNA family identifiers
5) Requires numeric, non-negative, integer copy numbers without missing values
6) Confirms that the species columns in the workbook exactly match the tree tips

Presence and Absence Classification:
1) Copy number greater than or equal to one is classified as present
2) Copy number equal to zero is classified as absent

MRCA Inference:
1) Families present in two or more species are assigned to the MRCA of all species containing the family
2) Families present in one species are assigned to that species tip
3) Families absent from all species receive a missing origin-node value

Input Files

1) pruned_tree.nwk
   A rooted Newick tree for the 62 species included in the analysis.

2) Copy_table.xlsx
   An Excel workbook containing miRNA family copy numbers. The script reads the worksheet named copy_number_table_62.

Required Workbook Structure:
1) miRNA_family: unique miRNA family identifier
2) One numeric copy-number column for each species in the tree
3) Species column names must match the Newick tip labels exactly

Output File

miRNA_ancestral_origins.csv

Output Columns:
1) miRNA: miRNA family identifier
2) origin_node: numerical tip or internal-node identifier used by the ape phylo object
3) n_present: number of species containing the miRNA family
4) n_total: total number of species assessed
5) present_proportion: proportion of assessed species containing the family
6) present_species: comma-separated list of species containing the family

Dependencies

1) R
2) ape
3) readxl

Usage

Place the script and both default input files in the same working directory, then run:

Rscript miRNA_locus_MRCA_inference.R

The script also accepts optional command-line paths in the following order:

Rscript miRNA_locus_MRCA_inference.R <tree_file> <copy_number_workbook> <output_csv>

Example:

Rscript miRNA_locus_MRCA_inference.R pruned_tree.nwk Copy_table.xlsx miRNA_ancestral_origins.csv

Interpretation Notes

1) The MRCA assignment assumes one origin of each miRNA family. Independent gains could place the reported MRCA deeper in the tree than the true origin of each independently gained locus.
2) Numerical node identifiers depend on the topology and ordering of the input tree and may change if the tree is pruned, rerooted, or reordered.
