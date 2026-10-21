# Comparison of MirGeneDB and homology-based miRNA copy numbers
#
# Analyses performed:
#
# 1. All families
#    Includes all miRNA families, including shared zeroes.
#
# 2. Nonzero-union families
#    Excludes only cases where both MirGeneDB and prediction
#    equal zero.
#
# 3. Recovery of curated loci
#    Measures how well the homology-based prediction recovers
#    MirGeneDB-supported loci without penalising additional
#    predicted loci.
#
# Copy-level recovery:
#
#   sum(min(Prediction, MirGeneDB)) / sum(MirGeneDB)
#
# Presence recovery:
#
#   MirGeneDB-positive family/species comparisons with at least
#   one predicted locus / all MirGeneDB-positive comparisons
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

library(ggplot2)
library(reshape2)


# ============================================================
# 2. Set working directory
# ============================================================

#setwd("SET YOUR OWN WORKING PATH")


# ============================================================
# 3. Read input files
# ============================================================

mirgenedb <- read.csv(
  "miRNA_copy_number_table_homology.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

prediction <- read.delim(
  "copy_number_table_19_mirgenedb.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 4. Standardise family-name columns
# ============================================================

colnames(mirgenedb)[1] <- "miRNA_family"
colnames(prediction)[1] <- "miRNA_family"

mirgenedb$miRNA_family <- trimws(
  mirgenedb$miRNA_family
)

prediction$miRNA_family <- trimws(
  prediction$miRNA_family
)


# Check for duplicated family names
if (anyDuplicated(mirgenedb$miRNA_family) > 0) {
  
  duplicated_names <- unique(
    mirgenedb$miRNA_family[
      duplicated(mirgenedb$miRNA_family)
    ]
  )
  
  stop(
    paste(
      "Duplicated miRNA family names in MirGeneDB table:",
      paste(duplicated_names, collapse = ", ")
    )
  )
}


if (anyDuplicated(prediction$miRNA_family) > 0) {
  
  duplicated_names <- unique(
    prediction$miRNA_family[
      duplicated(prediction$miRNA_family)
    ]
  )
  
  stop(
    paste(
      "Duplicated miRNA family names in prediction table:",
      paste(duplicated_names, collapse = ", ")
    )
  )
}


# ============================================================
# 5. Identify species present in both tables
# ============================================================

mirgenedb_species <- setdiff(
  colnames(mirgenedb),
  "miRNA_family"
)

prediction_species <- setdiff(
  colnames(prediction),
  "miRNA_family"
)

common_species <- intersect(
  mirgenedb_species,
  prediction_species
)


if (length(common_species) == 0) {
  stop(
    paste(
      "No matching species-column names were found between",
      "the MirGeneDB and prediction tables."
    )
  )
}


# Use species order file if available
if (file.exists("19_species_abbr.txt")) {
  
  species_order <- readLines(
    "19_species_abbr.txt",
    warn = FALSE
  )
  
  species_order <- trimws(
    species_order
  )
  
  species_order <- unique(
    species_order[nzchar(species_order)]
  )
  
  missing_species <- setdiff(
    species_order,
    common_species
  )
  
  if (length(missing_species) > 0) {
    warning(
      paste(
        "The following species in 19_species_abbr.txt",
        "were not present in both input tables:",
        paste(missing_species, collapse = ", ")
      )
    )
  }
  
  species <- species_order[
    species_order %in% common_species
  ]
  
} else {
  
  warning(
    paste(
      "19_species_abbr.txt was not found.",
      "Species will be analysed in the order present",
      "in the input tables."
    )
  )
  
  species <- common_species
}


if (length(species) == 0) {
  stop(
    "No species remained after matching the input tables."
  )
}


cat(
  "Species included in the analysis:",
  length(species),
  "\n"
)

print(species)


# ============================================================
# 6. Convert copy-number columns to numeric
# ============================================================

for (sp in species) {
  
  mirgenedb[[sp]] <- suppressWarnings(
    as.numeric(mirgenedb[[sp]])
  )
  
  prediction[[sp]] <- suppressWarnings(
    as.numeric(prediction[[sp]])
  )
}


# ============================================================
# 7. Rename species columns before merging
# ============================================================

mirgenedb_renamed <- mirgenedb[
  ,
  c("miRNA_family", species),
  drop = FALSE
]

prediction_renamed <- prediction[
  ,
  c("miRNA_family", species),
  drop = FALSE
]


colnames(mirgenedb_renamed)[
  colnames(mirgenedb_renamed) != "miRNA_family"
] <- paste0(
  species,
  "_MirGeneDB"
)


colnames(prediction_renamed)[
  colnames(prediction_renamed) != "miRNA_family"
] <- paste0(
  species,
  "_Prediction"
)


# ============================================================
# 8. Merge the two copy-number tables
# ============================================================

comparison_data <- merge(
  mirgenedb_renamed,
  prediction_renamed,
  by = "miRNA_family",
  all = TRUE,
  sort = FALSE
)


# Replace missing values with zero
copy_number_columns <- setdiff(
  colnames(comparison_data),
  "miRNA_family"
)

for (column_name in copy_number_columns) {
  
  comparison_data[[column_name]] <- suppressWarnings(
    as.numeric(comparison_data[[column_name]])
  )
  
  comparison_data[[column_name]][
    is.na(comparison_data[[column_name]])
  ] <- 0
}


cat(
  "Number of miRNA families after merging:",
  nrow(comparison_data),
  "\n"
)


# ============================================================
# 9. Function for copy-number comparison statistics
# ============================================================

calculate_statistics <- function(
    mirgenedb_counts,
    prediction_counts,
    remove_shared_zeros = FALSE
) {
  
  mirgenedb_counts <- as.numeric(
    mirgenedb_counts
  )
  
  prediction_counts <- as.numeric(
    prediction_counts
  )
  
  valid <- !is.na(mirgenedb_counts) &
    !is.na(prediction_counts)
  
  x <- mirgenedb_counts[valid]
  y <- prediction_counts[valid]
  
  
  # Remove only shared-zero observations
  if (remove_shared_zeros) {
    
    keep <- x > 0 | y > 0
    
    x <- x[keep]
    y <- y[keep]
  }
  
  
  total <- length(x)
  
  
  if (total == 0) {
    
    return(
      data.frame(
        Total = 0,
        Agreement_count = NA,
        Agreement_percent = NA,
        Difference_1_count = NA,
        Difference_1_percent = NA,
        Difference_2_5_count = NA,
        Difference_2_5_percent = NA,
        Difference_over_5_count = NA,
        Difference_over_5_percent = NA,
        Prediction_higher_count = NA,
        Prediction_higher_percent_total = NA,
        Prediction_lower_count = NA,
        Prediction_lower_percent_total = NA,
        Prediction_higher_percent_disagreements = NA,
        Prediction_lower_percent_disagreements = NA,
        Spearman_rho = NA,
        Spearman_p_value = NA,
        Pearson_log2 = NA,
        Pearson_log2_p_value = NA
      )
    )
  }
  
  
  difference <- y - x
  
  absolute_difference <- abs(
    difference
  )
  
  
  agreement_count <- sum(
    difference == 0
  )
  
  difference_1_count <- sum(
    absolute_difference == 1
  )
  
  difference_2_5_count <- sum(
    absolute_difference >= 2 &
      absolute_difference <= 5
  )
  
  difference_over_5_count <- sum(
    absolute_difference > 5
  )
  
  prediction_higher_count <- sum(
    difference > 0
  )
  
  prediction_lower_count <- sum(
    difference < 0
  )
  
  disagreement_count <- sum(
    difference != 0
  )
  
  
  # Spearman correlation
  if (
    length(x) > 1 &&
    length(unique(x)) > 1 &&
    length(unique(y)) > 1
  ) {
    
    spearman_test <- suppressWarnings(
      cor.test(
        x,
        y,
        method = "spearman",
        exact = FALSE
      )
    )
    
    spearman_rho <- unname(
      spearman_test$estimate
    )
    
    spearman_p_value <- spearman_test$p.value
    
  } else {
    
    spearman_rho <- NA
    spearman_p_value <- NA
  }
  
  
  # Pearson correlation after log2(count + 1)
  log_x <- log2(
    x + 1
  )
  
  log_y <- log2(
    y + 1
  )
  
  
  if (
    length(log_x) > 1 &&
    length(unique(log_x)) > 1 &&
    length(unique(log_y)) > 1
  ) {
    
    pearson_test <- cor.test(
      log_x,
      log_y,
      method = "pearson"
    )
    
    pearson_log2 <- unname(
      pearson_test$estimate
    )
    
    pearson_log2_p_value <- pearson_test$p.value
    
  } else {
    
    pearson_log2 <- NA
    pearson_log2_p_value <- NA
  }
  
  
  data.frame(
    Total = total,
    
    Agreement_count = agreement_count,
    
    Agreement_percent =
      100 * agreement_count / total,
    
    Difference_1_count =
      difference_1_count,
    
    Difference_1_percent =
      100 * difference_1_count / total,
    
    Difference_2_5_count =
      difference_2_5_count,
    
    Difference_2_5_percent =
      100 * difference_2_5_count / total,
    
    Difference_over_5_count =
      difference_over_5_count,
    
    Difference_over_5_percent =
      100 * difference_over_5_count / total,
    
    Prediction_higher_count =
      prediction_higher_count,
    
    Prediction_higher_percent_total =
      100 * prediction_higher_count / total,
    
    Prediction_lower_count =
      prediction_lower_count,
    
    Prediction_lower_percent_total =
      100 * prediction_lower_count / total,
    
    Prediction_higher_percent_disagreements =
      ifelse(
        disagreement_count > 0,
        100 * prediction_higher_count /
          disagreement_count,
        NA
      ),
    
    Prediction_lower_percent_disagreements =
      ifelse(
        disagreement_count > 0,
        100 * prediction_lower_count /
          disagreement_count,
        NA
      ),
    
    Spearman_rho =
      spearman_rho,
    
    Spearman_p_value =
      spearman_p_value,
    
    Pearson_log2 =
      pearson_log2,
    
    Pearson_log2_p_value =
      pearson_log2_p_value
  )
}


# ============================================================
# 10. Function for recovery statistics
# ============================================================

calculate_recovery <- function(
    mirgenedb_counts,
    prediction_counts
) {
  
  mirgenedb_counts <- as.numeric(
    mirgenedb_counts
  )
  
  prediction_counts <- as.numeric(
    prediction_counts
  )
  
  valid <- !is.na(mirgenedb_counts) &
    !is.na(prediction_counts)
  
  x <- mirgenedb_counts[valid]
  y <- prediction_counts[valid]
  
  
  # Copy-level recovery
  total_curated_loci <- sum(
    x
  )
  
  recovered_curated_loci <- sum(
    pmin(x, y)
  )
  
  missing_curated_loci <-
    total_curated_loci -
    recovered_curated_loci
  
  copy_recovery_percent <- ifelse(
    total_curated_loci > 0,
    100 * recovered_curated_loci /
      total_curated_loci,
    NA
  )
  
  
  # Presence-level recovery
  mirgenedb_positive <- x > 0
  
  mirgenedb_positive_comparisons <- sum(
    mirgenedb_positive
  )
  
  presence_recovered_comparisons <- sum(
    mirgenedb_positive &
      y > 0
  )
  
  presence_missed_comparisons <-
    mirgenedb_positive_comparisons -
    presence_recovered_comparisons
  
  presence_recovery_percent <- ifelse(
    mirgenedb_positive_comparisons > 0,
    100 * presence_recovered_comparisons /
      mirgenedb_positive_comparisons,
    NA
  )
  
  
  # Full-copy recovery at comparison level
  fully_recovered_comparisons <- sum(
    mirgenedb_positive &
      y >= x
  )
  
  partially_recovered_comparisons <- sum(
    mirgenedb_positive &
      y > 0 &
      y < x
  )
  
  completely_missed_comparisons <- sum(
    mirgenedb_positive &
      y == 0
  )
  
  full_copy_recovery_percent <- ifelse(
    mirgenedb_positive_comparisons > 0,
    100 * fully_recovered_comparisons /
      mirgenedb_positive_comparisons,
    NA
  )
  
  
  data.frame(
    Total_curated_loci =
      total_curated_loci,
    
    Recovered_curated_loci =
      recovered_curated_loci,
    
    Missing_curated_loci =
      missing_curated_loci,
    
    Copy_recovery_percent =
      copy_recovery_percent,
    
    MirGeneDB_positive_comparisons =
      mirgenedb_positive_comparisons,
    
    Presence_recovered_comparisons =
      presence_recovered_comparisons,
    
    Presence_missed_comparisons =
      presence_missed_comparisons,
    
    Presence_recovery_percent =
      presence_recovery_percent,
    
    Fully_recovered_comparisons =
      fully_recovered_comparisons,
    
    Partially_recovered_comparisons =
      partially_recovered_comparisons,
    
    Completely_missed_comparisons =
      completely_missed_comparisons,
    
    Full_copy_recovery_percent =
      full_copy_recovery_percent
  )
}


# ============================================================
# 11. Species-level analyses
# ============================================================

species_results <- list()

recovery_results <- list()

family_level_results <- list()


for (sp in species) {
  
  mirgenedb_column <- paste0(
    sp,
    "_MirGeneDB"
  )
  
  prediction_column <- paste0(
    sp,
    "_Prediction"
  )
  
  
  mirgenedb_counts <- as.numeric(
    comparison_data[[mirgenedb_column]]
  )
  
  prediction_counts <- as.numeric(
    comparison_data[[prediction_column]]
  )
  
  
  valid <- !is.na(mirgenedb_counts) &
    !is.na(prediction_counts)
  
  
  shared_zero <- valid &
    mirgenedb_counts == 0 &
    prediction_counts == 0
  
  
  nonzero_union <- valid &
    (
      mirgenedb_counts > 0 |
        prediction_counts > 0
    )
  
  
  # Analysis 1: all families
  all_statistics <- calculate_statistics(
    mirgenedb_counts = mirgenedb_counts,
    prediction_counts = prediction_counts,
    remove_shared_zeros = FALSE
  )
  
  all_statistics$Species <- sp
  
  all_statistics$Analysis <-
    "All families"
  
  
  # Analysis 2: excluding shared zeroes
  nonzero_statistics <- calculate_statistics(
    mirgenedb_counts = mirgenedb_counts,
    prediction_counts = prediction_counts,
    remove_shared_zeros = TRUE
  )
  
  nonzero_statistics$Species <- sp
  
  nonzero_statistics$Analysis <-
    "At least one method nonzero"
  
  
  species_results[[paste0(sp, "_all")]] <- all_statistics
  
  
  species_results[[paste0(sp, "_nonzero")]] <- nonzero_statistics
  
  
  # Recovery statistics
  recovery_statistics <- calculate_recovery(
    mirgenedb_counts = mirgenedb_counts,
    prediction_counts = prediction_counts
  )
  
  recovery_statistics$Species <- sp
  
  recovery_results[[sp]] <-
    recovery_statistics
  
  
  # Family-level output
  family_level_results[[sp]] <- data.frame(
    Species = sp,
    
    miRNA_family =
      comparison_data$miRNA_family,
    
    MirGeneDB =
      mirgenedb_counts,
    
    Prediction =
      prediction_counts,
    
    Difference =
      prediction_counts -
      mirgenedb_counts,
    
    Absolute_difference =
      abs(
        prediction_counts -
          mirgenedb_counts
      ),
    
    Shared_zero =
      shared_zero,
    
    Included_in_nonzero_analysis =
      nonzero_union,
    
    MirGeneDB_positive =
      mirgenedb_counts > 0,
    
    Presence_recovered =
      mirgenedb_counts > 0 &
      prediction_counts > 0,
    
    Full_copy_number_recovered =
      mirgenedb_counts > 0 &
      prediction_counts >= mirgenedb_counts,
    
    Curated_copies_recovered =
      pmin(
        mirgenedb_counts,
        prediction_counts
      ),
    
    Curated_copies_missing =
      pmax(
        mirgenedb_counts -
          prediction_counts,
        0
      ),
    
    Additional_predicted_copies =
      pmax(
        prediction_counts -
          mirgenedb_counts,
        0
      ),
    
    Comparison = ifelse(
      prediction_counts == mirgenedb_counts,
      "Identical copy number",
      ifelse(
        prediction_counts > mirgenedb_counts,
        "Additional predicted loci",
        "Fewer predicted loci than curated"
      )
    ),
    
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 12. Save species-level comparison statistics
# ============================================================

species_summary <- do.call(
  rbind,
  species_results
)

rownames(species_summary) <- NULL


species_summary <- species_summary[
  ,
  c(
    "Species",
    "Analysis",
    "Total",
    "Agreement_count",
    "Agreement_percent",
    "Difference_1_count",
    "Difference_1_percent",
    "Difference_2_5_count",
    "Difference_2_5_percent",
    "Difference_over_5_count",
    "Difference_over_5_percent",
    "Prediction_higher_count",
    "Prediction_higher_percent_total",
    "Prediction_lower_count",
    "Prediction_lower_percent_total",
    "Prediction_higher_percent_disagreements",
    "Prediction_lower_percent_disagreements",
    "Spearman_rho",
    "Spearman_p_value",
    "Pearson_log2",
    "Pearson_log2_p_value"
  )
]


write.csv(
  species_summary,
  "copy_number_comparison_by_species.csv",
  row.names = FALSE
)


# ============================================================
# 13. Save species-level recovery statistics
# ============================================================

recovery_summary <- do.call(
  rbind,
  recovery_results
)

rownames(recovery_summary) <- NULL


recovery_summary <- recovery_summary[
  ,
  c(
    "Species",
    "Total_curated_loci",
    "Recovered_curated_loci",
    "Missing_curated_loci",
    "Copy_recovery_percent",
    "MirGeneDB_positive_comparisons",
    "Presence_recovered_comparisons",
    "Presence_missed_comparisons",
    "Presence_recovery_percent",
    "Fully_recovered_comparisons",
    "Partially_recovered_comparisons",
    "Completely_missed_comparisons",
    "Full_copy_recovery_percent"
  )
]


write.csv(
  recovery_summary,
  "curated_locus_recovery_by_species.csv",
  row.names = FALSE
)


# ============================================================
# 14. Save family-level comparison
# ============================================================

family_level_output <- do.call(
  rbind,
  family_level_results
)

rownames(family_level_output) <- NULL


write.csv(
  family_level_output,
  "copy_number_comparison_family_level.csv",
  row.names = FALSE
)


# ============================================================
# 15. Shared-zero summary by species
# ============================================================

shared_zero_summary <- do.call(
  rbind,
  lapply(
    species,
    function(sp) {
      
      mirgenedb_counts <- as.numeric(
        comparison_data[[paste0(sp, "_MirGeneDB")]]
      )
      
      prediction_counts <- as.numeric(
        comparison_data[[paste0(sp, "_Prediction")]]
      )
      
      
      valid <- !is.na(mirgenedb_counts) &
        !is.na(prediction_counts)
      
      
      total_count <- sum(
        valid
      )
      
      
      shared_zero_count <- sum(
        valid &
          mirgenedb_counts == 0 &
          prediction_counts == 0
      )
      
      
      nonzero_union_count <- sum(
        valid &
          (
            mirgenedb_counts > 0 |
              prediction_counts > 0
          )
      )
      
      
      data.frame(
        Species = sp,
        
        Total_families =
          total_count,
        
        Shared_zero_count =
          shared_zero_count,
        
        Shared_zero_percent =
          ifelse(
            total_count > 0,
            100 * shared_zero_count /
              total_count,
            NA
          ),
        
        Nonzero_union_count =
          nonzero_union_count,
        
        Nonzero_union_percent =
          ifelse(
            total_count > 0,
            100 * nonzero_union_count /
              total_count,
            NA
          )
      )
    }
  )
)


rownames(shared_zero_summary) <- NULL


write.csv(
  shared_zero_summary,
  "shared_zero_summary_by_species.csv",
  row.names = FALSE
)


# ============================================================
# 16. Pooled vectors across all species
# ============================================================

pooled_mirgenedb <- unlist(
  lapply(
    species,
    function(sp) {
      
      as.numeric(
        comparison_data[[paste0(sp, "_MirGeneDB")]]
      )
    }
  ),
  use.names = FALSE
)


pooled_prediction <- unlist(
  lapply(
    species,
    function(sp) {
      
      as.numeric(
        comparison_data[[paste0(sp, "_Prediction")]]
      )
    }
  ),
  use.names = FALSE
)


# ============================================================
# 17. Pooled comparison statistics
# ============================================================

pooled_all <- calculate_statistics(
  mirgenedb_counts =
    pooled_mirgenedb,
  
  prediction_counts =
    pooled_prediction,
  
  remove_shared_zeros =
    FALSE
)

pooled_all$Analysis <-
  "All species-family comparisons"


pooled_nonzero <- calculate_statistics(
  mirgenedb_counts =
    pooled_mirgenedb,
  
  prediction_counts =
    pooled_prediction,
  
  remove_shared_zeros =
    TRUE
)

pooled_nonzero$Analysis <-
  "Species-family comparisons with at least one nonzero value"


pooled_summary <- rbind(
  pooled_all,
  pooled_nonzero
)


pooled_total <- sum(
  !is.na(pooled_mirgenedb) &
    !is.na(pooled_prediction)
)


pooled_shared_zero_count <- sum(
  pooled_mirgenedb == 0 &
    pooled_prediction == 0,
  na.rm = TRUE
)


pooled_summary$Shared_zero_count <-
  pooled_shared_zero_count


pooled_summary$Shared_zero_percent <-
  ifelse(
    pooled_total > 0,
    100 * pooled_shared_zero_count /
      pooled_total,
    NA
  )


write.csv(
  pooled_summary,
  "copy_number_comparison_pooled.csv",
  row.names = FALSE
)


# ============================================================
# 18. Pooled recovery statistics
# ============================================================

pooled_recovery <- calculate_recovery(
  mirgenedb_counts =
    pooled_mirgenedb,
  
  prediction_counts =
    pooled_prediction
)


write.csv(
  pooled_recovery,
  "curated_locus_recovery_pooled.csv",
  row.names = FALSE
)


# ============================================================
# 19. Prepare plotting factors
# ============================================================

recovery_summary$Species <- factor(
  recovery_summary$Species,
  levels = rev(species)
)


# ============================================================
# 20. Plot copy-level recovery by species
# ============================================================

copy_recovery_plot <- ggplot(
  recovery_summary,
  aes(
    x = Species,
    y = Copy_recovery_percent
  )
) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Curated copies recovered (%)"
  ) +
  theme_classic() +
  theme(
    axis.text =
      element_text(size = 10)
  )


print(
  copy_recovery_plot
)


ggsave(
  filename =
    "curated_copy_recovery_by_species.pdf",
  
  plot =
    copy_recovery_plot,
  
  width =
    9,
  
  height =
    8
)


ggsave(
  filename =
    "curated_copy_recovery_by_species.png",
  
  plot =
    copy_recovery_plot,
  
  width =
    9,
  
  height =
    8,
  
  dpi =
    300
)


# ============================================================
# 21. Plot presence recovery by species
# ============================================================

presence_recovery_plot <- ggplot(
  recovery_summary,
  aes(
    x = Species,
    y = Presence_recovery_percent
  )
) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "MirGeneDB-positive families recovered (%)"
  ) +
  theme_classic() +
  theme(
    axis.text =
      element_text(size = 10)
  )


print(
  presence_recovery_plot
)


ggsave(
  filename =
    "presence_recovery_by_species.pdf",
  
  plot =
    presence_recovery_plot,
  
  width =
    9,
  
  height =
    8
)


ggsave(
  filename =
    "presence_recovery_by_species.png",
  
  plot =
    presence_recovery_plot,
  
  width =
    9,
  
  height =
    8,
  
  dpi =
    300
)


# ============================================================
# 22. Plot presence and full-copy recovery together
# ============================================================

recovery_plot_data <- recovery_summary[
  ,
  c(
    "Species",
    "Presence_recovery_percent",
    "Full_copy_recovery_percent"
  )
]


recovery_plot_long <- melt(
  recovery_plot_data,
  id.vars = "Species",
  variable.name = "Recovery_type",
  value.name = "Recovery_percent"
)


recovery_plot_long$Recovery_type <- factor(
  recovery_plot_long$Recovery_type,
  levels = c(
    "Presence_recovery_percent",
    "Full_copy_recovery_percent"
  ),
  labels = c(
    "At least one locus recovered",
    "Full curated copy number recovered"
  )
)


combined_recovery_plot <- ggplot(
  recovery_plot_long,
  aes(
    x = Species,
    y = Recovery_percent,
    fill = Recovery_type
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Recovery (%)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(
      size = 10
    )
  )


print(
  combined_recovery_plot
)


ggsave(
  filename =
    "presence_and_copy_number_recovery.pdf",
  
  plot =
    combined_recovery_plot,
  
  width =
    10,
  
  height =
    8
)


ggsave(
  filename =
    "presence_and_copy_number_recovery.png",
  
  plot =
    combined_recovery_plot,
  
  width =
    10,
  
  height =
    8,
  
  dpi =
    300
)


# ============================================================
# 23. Print key results
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "SPECIES-LEVEL RECOVERY\n"
)

cat(
  "============================================================\n\n"
)


print(
  recovery_summary
)


cat(
  "\n============================================================\n"
)

cat(
  "POOLED RECOVERY\n"
)

cat(
  "============================================================\n\n"
)


print(
  pooled_recovery
)


cat(
  "\nOverall copy-level recovery:\n"
)

cat(
  pooled_recovery$Recovered_curated_loci,
  "of",
  pooled_recovery$Total_curated_loci,
  "curated copies recovered (",
  round(
    pooled_recovery$Copy_recovery_percent,
    2
  ),
  "%)\n",
  sep = " "
)


cat(
  "\nOverall presence recovery:\n"
)

cat(
  pooled_recovery$Presence_recovered_comparisons,
  "of",
  pooled_recovery$MirGeneDB_positive_comparisons,
  "MirGeneDB-positive family-species comparisons recovered (",
  round(
    pooled_recovery$Presence_recovery_percent,
    2
  ),
  "%)\n",
  sep = " "
)


cat(
  "\nFull curated copy-number recovery:\n"
)

cat(
  pooled_recovery$Fully_recovered_comparisons,
  "of",
  pooled_recovery$MirGeneDB_positive_comparisons,
  "MirGeneDB-positive family-species comparisons had",
  "prediction >= curated copy number (",
  round(
    pooled_recovery$Full_copy_recovery_percent,
    2
  ),
  "%)\n",
  sep = " "
)


cat(
  "\nAnalysis complete.\n"
)

cat(
  "Output files were written to:\n",
  getwd(),
  "\n"
)
