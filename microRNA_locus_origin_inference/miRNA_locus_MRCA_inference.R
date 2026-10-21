# ============================================================================
# MRCA NODE INFERENCE FOR SYNTENY-VERIFIED miRNAs
# Method: MRCA of all species with the miRNA
# Note: The analysis reports MRCA nodes, not evolutionary age estimates
# ============================================================================

library(ape)
library(readxl)

# ----------------------------------------------------------------------------
# 1. Input and output files
# ----------------------------------------------------------------------------
# Optional command-line arguments:
#   1) tree file
#   2) Excel copy-number table
#   3) output CSV file

args <- commandArgs(trailingOnly = TRUE)

tree_file <- if (length(args) >= 1) args[1] else "pruned_tree.nwk"
copy_table_file <- if (length(args) >= 2) args[2] else "Copy_table.xlsx"
output_file <- if (length(args) >= 3) args[3] else "miRNA_ancestral_origins.csv"
sheet_name <- "copy_number_table_62"

if (!file.exists(tree_file)) {
  stop("Tree file not found: ", tree_file)
}

if (!file.exists(copy_table_file)) {
  stop("Copy-number table not found: ", copy_table_file)
}


# ----------------------------------------------------------------------------
# 2. Read and validate the phylogenetic tree
# ----------------------------------------------------------------------------

tree <- read.tree(tree_file)

if (is.null(tree) || Ntip(tree) < 2) {
  stop("The tree could not be read or contains fewer than two tips.")
}

if (!is.rooted(tree)) {
  stop("The input tree must be rooted for MRCA inference.")
}

if (anyDuplicated(tree$tip.label)) {
  stop("The tree contains duplicated tip labels.")
}


# ----------------------------------------------------------------------------
# 3. Read and validate the copy-number table
# ----------------------------------------------------------------------------

available_sheets <- excel_sheets(copy_table_file)

if (!(sheet_name %in% available_sheets)) {
  stop(
    "Worksheet '", sheet_name, "' was not found. Available worksheets: ",
    paste(available_sheets, collapse = ", ")
  )
}

input_table <- read_excel(
  copy_table_file,
  sheet = sheet_name
)

if (!("miRNA_family" %in% names(input_table))) {
  stop("The worksheet must contain a 'miRNA_family' column.")
}

miRNA_ids <- input_table$miRNA_family

if (any(is.na(miRNA_ids) | trimws(miRNA_ids) == "")) {
  stop("The 'miRNA_family' column contains missing or empty identifiers.")
}

if (anyDuplicated(miRNA_ids)) {
  duplicated_ids <- unique(miRNA_ids[duplicated(miRNA_ids)])
  stop(
    "Duplicated miRNA family identifiers: ",
    paste(duplicated_ids, collapse = ", ")
  )
}

copy_table <- as.data.frame(
  input_table[, setdiff(names(input_table), "miRNA_family"), drop = FALSE],
  check.names = FALSE
)

if (!all(vapply(copy_table, is.numeric, logical(1)))) {
  non_numeric <- names(copy_table)[
    !vapply(copy_table, is.numeric, logical(1))
  ]
  stop(
    "Non-numeric species columns: ",
    paste(non_numeric, collapse = ", ")
  )
}

copy_matrix <- as.matrix(copy_table)

if (anyNA(copy_matrix)) {
  stop("The copy-number table contains missing values.")
}

if (any(copy_matrix < 0)) {
  stop("The copy-number table contains negative values.")
}

if (any(abs(copy_matrix - round(copy_matrix)) > sqrt(.Machine$double.eps))) {
  stop("The copy-number table contains non-integer values.")
}


# ----------------------------------------------------------------------------
# 4. Match the copy-number table to the tree
# ----------------------------------------------------------------------------

missing_species <- setdiff(tree$tip.label, colnames(copy_matrix))
extra_species <- setdiff(colnames(copy_matrix), tree$tip.label)

if (length(missing_species) > 0) {
  stop(
    "Species present in the tree but missing from the table: ",
    paste(missing_species, collapse = ", ")
  )
}

if (length(extra_species) > 0) {
  stop(
    "Species present in the table but missing from the tree: ",
    paste(extra_species, collapse = ", ")
  )
}

copy_matrix <- copy_matrix[, tree$tip.label, drop = FALSE]
rownames(copy_matrix) <- miRNA_ids

# Copy number >= 1 indicates presence.
presence_matrix <- (copy_matrix > 0) * 1L


# ----------------------------------------------------------------------------
# 5. Infer the MRCA node for each miRNA family
# ----------------------------------------------------------------------------

results <- data.frame(
  miRNA = rownames(presence_matrix),
  origin_node = NA_integer_,
  n_present = NA_integer_,
  n_total = ncol(presence_matrix),
  present_proportion = NA_real_,
  present_species = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(presence_matrix))) {
  presence_vec <- presence_matrix[i, ]
  present_species <- colnames(presence_matrix)[presence_vec == 1L]
  n_present <- length(present_species)

  results$n_present[i] <- n_present
  results$present_proportion[i] <- n_present / ncol(presence_matrix)

  if (n_present > 0) {
    results$present_species[i] <- paste(present_species, collapse = ", ")
  }

  if (n_present == 0) {
    next
  }

  if (n_present == 1) {
    # With one observed species, report that species' tip number.
    results$origin_node[i] <- match(present_species, tree$tip.label)
    next
  }

  present_nodes <- match(present_species, tree$tip.label)
  mrca_node <- getMRCA(tree, present_nodes)

  results$origin_node[i] <- mrca_node
}


# ----------------------------------------------------------------------------
# 6. Export results
# ----------------------------------------------------------------------------

write.csv(
  results,
  output_file,
  row.names = FALSE,
  na = ""
)

cat(
  "Analyzed", nrow(results), "miRNA families across",
  ncol(presence_matrix), "species.\n",
  "Results written to:", output_file, "\n"
)
