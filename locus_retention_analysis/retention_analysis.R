library(readxl)
library(dplyr)
library(ggplot2)
library(broom)


# ------------------------------------------------------------
# 1. Read the Excel file
# ------------------------------------------------------------

dat_raw <- read_excel(
  "eutheria_intronic.xlsx",
  sheet = "Sheet1"
)

# Inspect the imported data
head(dat_raw)
names(dat_raw)
dim(dat_raw)
str(dat_raw)


# ------------------------------------------------------------
# 2. Examine genomic-location categories
# ------------------------------------------------------------

table(dat_raw$location, useNA = "ifany")


# ------------------------------------------------------------
# 3. Retain only unambiguous intronic and intergenic loci
# ------------------------------------------------------------

dat <- dat_raw %>%
  filter(location %in% c("intronic", "intergenic")) %>%
  mutate(
    # Intergenic is the reference category in the GLM
    location = factor(
      location,
      levels = c("intergenic", "intronic")
    ),
    
    # Number of species in which the locus was not retained
    lost_locus = all_synteny_species - retention_locus,
    
    # Recalculate the retention rate to verify the input column
    retention_rate_check =
      retention_locus / all_synteny_species
  )


# Check the retained dataset
dim(dat)
table(dat$location)
head(dat)


# Verify that all families have 58 assessable species
unique(dat$all_synteny_species)

stopifnot(
  all(dat$all_synteny_species == 58)
)


# Check that the supplied retention rate is correct
all.equal(
  dat$retention_rate,
  dat$retention_rate_check
)


# Check for impossible values
stopifnot(
  all(dat$retention_locus >= 0),
  all(dat$lost_locus >= 0),
  all(dat$retention_locus <= dat$all_synteny_species)
)


# ------------------------------------------------------------
# 4. Descriptive statistics
# ------------------------------------------------------------

descriptive_stats <- dat %>%
  group_by(location) %>%
  summarise(
    n_families = n(),
    
    mean_retained_loci =
      mean(retention_locus),
    
    sd_retained_loci =
      sd(retention_locus),
    
    median_retained_loci =
      median(retention_locus),
    
    min_retained_loci =
      min(retention_locus),
    
    max_retained_loci =
      max(retention_locus),
    
    mean_retention_rate =
      mean(retention_rate_check),
    
    sd_retention_rate =
      sd(retention_rate_check),
    
    median_retention_rate =
      median(retention_rate_check),
    
    IQR_retention_rate =
      IQR(retention_rate_check),
    
    .groups = "drop"
  )

print(descriptive_stats)


# Optional: display proportions as percentages
descriptive_stats_percent <- descriptive_stats %>%
  mutate(
    mean_retention_percent =
      100 * mean_retention_rate,
    
    median_retention_percent =
      100 * median_retention_rate
  )

print(descriptive_stats_percent)


# ------------------------------------------------------------
# 5. Fit an ordinary binomial GLM
# ------------------------------------------------------------
# This model is fitted mainly to evaluate overdispersion.
#
# Response:
#   retained loci versus lost loci
#
# Predictor:
#   location, with intergenic as the reference category
# ------------------------------------------------------------

fit_binomial <- glm(
  cbind(retention_locus, lost_locus) ~ location,
  family = binomial(link = "logit"),
  data = dat
)

summary(fit_binomial)


# ------------------------------------------------------------
# 6. Evaluate overdispersion
# ------------------------------------------------------------

pearson_chisq <- sum(
  residuals(fit_binomial, type = "pearson")^2
)

residual_df <- df.residual(fit_binomial)

dispersion_ratio <- pearson_chisq / residual_df

overdispersion_results <- data.frame(
  Pearson_chisquare = pearson_chisq,
  residual_df = residual_df,
  dispersion_ratio = dispersion_ratio
)

print(overdispersion_results)


# A dispersion ratio substantially greater than 1 indicates
# that the ordinary binomial variance is too small.
#
# In that case, use a quasibinomial GLM.


# ------------------------------------------------------------
# 7. Fit the primary quasibinomial GLM
# ------------------------------------------------------------

fit_quasi <- glm(
  cbind(retention_locus, lost_locus) ~ location,
  family = quasibinomial(link = "logit"),
  data = dat
)

summary(fit_quasi)


# ------------------------------------------------------------
# 8. Extract the coefficient table
# ------------------------------------------------------------

coefficient_table <- coef(summary(fit_quasi))

print(coefficient_table)


# Extract the intronic-location coefficient
beta_intronic <- coefficient_table[
  "locationintronic",
  "Estimate"
]

se_intronic <- coefficient_table[
  "locationintronic",
  "Std. Error"
]

t_intronic <- coefficient_table[
  "locationintronic",
  "t value"
]

p_intronic <- coefficient_table[
  "locationintronic",
  "Pr(>|t|)"
]


# ------------------------------------------------------------
# 9. Convert the coefficient to an odds ratio
# ------------------------------------------------------------

odds_ratio <- exp(beta_intronic)

odds_ratio


# Interpretation:
#
# OR > 1:
#   higher estimated retention odds for intronic loci
#
# OR < 1:
#   lower estimated retention odds for intronic loci
#
# OR = 1:
#   no difference


# ------------------------------------------------------------
# 10. Calculate a 95% confidence interval
# ------------------------------------------------------------
# Because this is a quasibinomial model, use the t distribution
# and the residual degrees of freedom.
# ------------------------------------------------------------

df_residual <- df.residual(fit_quasi)

critical_t <- qt(
  0.975,
  df = df_residual
)

beta_CI_lower <- beta_intronic -
  critical_t * se_intronic

beta_CI_upper <- beta_intronic +
  critical_t * se_intronic

OR_CI_lower <- exp(beta_CI_lower)
OR_CI_upper <- exp(beta_CI_upper)


# Combined model-results table
model_results <- data.frame(
  term = "Intronic versus intergenic",
  beta = beta_intronic,
  standard_error = se_intronic,
  t_value = t_intronic,
  residual_df = df_residual,
  p_value = p_intronic,
  odds_ratio = odds_ratio,
  OR_95CI_lower = OR_CI_lower,
  OR_95CI_upper = OR_CI_upper,
  dispersion_parameter = summary(fit_quasi)$dispersion
)

print(model_results)


# ------------------------------------------------------------
# 11. Predicted retention probabilities
# ------------------------------------------------------------
# Predicted probabilities for an average locus in each category
# according to the fitted model.
# ------------------------------------------------------------

new_dat <- data.frame(
  location = factor(
    c("intergenic", "intronic"),
    levels = levels(dat$location)
  )
)

predicted_link <- predict(
  fit_quasi,
  newdata = new_dat,
  type = "link",
  se.fit = TRUE
)

new_dat <- new_dat %>%
  mutate(
    predicted_log_odds =
      as.numeric(predicted_link$fit),
    
    standard_error =
      as.numeric(predicted_link$se.fit),
    
    predicted_probability =
      plogis(predicted_log_odds),
    
    lower_probability =
      plogis(
        predicted_log_odds -
          critical_t * standard_error
      ),
    
    upper_probability =
      plogis(
        predicted_log_odds +
          critical_t * standard_error
      )
  )

print(new_dat)


# ------------------------------------------------------------
# 12. Wilcoxon rank-sum sensitivity analysis
# ------------------------------------------------------------
# This compares the family-level retention-rate distributions
# without assuming a binomial variance structure.
# ------------------------------------------------------------

wilcox_result <- wilcox.test(
  retention_rate_check ~ location,
  data = dat,
  exact = FALSE,
  conf.int = TRUE
)

print(wilcox_result)


# ------------------------------------------------------------
# 13. Optional effect-size difference in descriptive means
# ------------------------------------------------------------

mean_rates <- dat %>%
  group_by(location) %>%
  summarise(
    mean_retention =
      mean(retention_rate_check),
    .groups = "drop"
  )

mean_rates

mean_intronic <- mean_rates$mean_retention[
  mean_rates$location == "intronic"
]

mean_intergenic <- mean_rates$mean_retention[
  mean_rates$location == "intergenic"
]

mean_difference <- mean_intronic - mean_intergenic

mean_difference


# ------------------------------------------------------------
# 14. Plot the observed retention distributions
# ------------------------------------------------------------

retention_plot <- ggplot(
  dat,
  aes(
    x = location,
    y = retention_rate_check
  )
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.8,
    alpha = 0.65
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    x = "Genomic location",
    y = "Genomic locus retention",
    title = "Retention of Eutherian-origin miRNA-like loci"
  ) +
  theme_classic(
    base_size = 12
  )

print(retention_plot)


# Save the plot
ggsave(
  filename = "eutheria_retention_intronic_intergenic.pdf",
  plot = retention_plot,
  width = 5,
  height = 4
)


# ------------------------------------------------------------
# 15. Export results
# ------------------------------------------------------------

write.csv(
  descriptive_stats,
  "eutheria_retention_descriptive_statistics.csv",
  row.names = FALSE
)

write.csv(
  model_results,
  "eutheria_retention_quasibinomial_results.csv",
  row.names = FALSE
)

write.csv(
  new_dat,
  "eutheria_retention_predicted_probabilities.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Compact reporting output
# ------------------------------------------------------------

cat(
  "\nNumber of families:",
  nrow(dat),
  "\nIntronic families:",
  sum(dat$location == "intronic"),
  "\nIntergenic families:",
  sum(dat$location == "intergenic"),
  "\n\nQuasibinomial GLM:",
  "\nbeta =",
  round(beta_intronic, 4),
  "\nSE =",
  round(se_intronic, 4),
  "\nt =",
  round(t_intronic, 3),
  "\ndf =",
  df_residual,
  "\nP =",
  signif(p_intronic, 4),
  "\nOR =",
  round(odds_ratio, 3),
  "\n95% CI for OR =",
  paste0(
    round(OR_CI_lower, 3),
    "–",
    round(OR_CI_upper, 3)
  ),
  "\nDispersion =",
  round(summary(fit_quasi)$dispersion, 3),
  "\n\nWilcoxon P =",
  signif(wilcox_result$p.value, 4),
  "\n"
)