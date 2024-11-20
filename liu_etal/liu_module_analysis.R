# Load required packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(stats)
library(multcomp) # For FDR adjustment

# Define the project folder and file paths
project_folder <- file.path("/fh/fast/gilbert_p/fg_data", "SEATRAC", "TB_hackday_2024")

meta_fn <- file.path(project_folder, "processed_data", "liu_etal_metadata.csv")
module_fn <- file.path(project_folder, "processed_data", "liu_etal_modules.csv")
cts_fn <- file.path(project_folder, "processed_data", "liu_etal_counts_long.csv")
immune_fn <- file.path(project_folder, "processed_data", "darrah_dose_immune_markers.csv")

# Load the data
meta <- read_csv(meta_fn)
modules <- read_csv(module_fn)
ncts <- read_csv(cts_fn)

# Load immune data and rename 'subjid' to 'animalid'
imm <- read_csv(immune_fn) %>%
  rename(animalid = subjid)

# Merge module names with the gene count data
# This allows us to compute module scores for each module
ncts <- ncts %>%
  left_join(modules, by = "gene")

# Compute module scores (average expression of genes in each module using log2-counts)
mod_scores <- ncts %>%
  group_by(sampleid, module) %>%
  summarise(count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  left_join(meta, by = "sampleid")

# Quick plot of module scores by visit pre/post-BCG
ggplot(mod_scores, aes(x = module, y = count, fill = visit)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pre" = "blue", "d2" = "red", "wk2" = "green", "wk4" = "purple", "wk12" = "orange")) +
  labs(title = "Module Scores by Visit", x = "Module", y = "Average Expression") +
  theme_minimal()

# Estimate how M1 module scores at day 2 correlate with various immune responses
mod <- "M1"
mod_visit <- "d2" # Options ['pre', 'd2', 'wk2', 'wk4', 'wk12']

# Filter immune response data to match module scores and calculate correlations
results <- imm %>%
  group_by(key) %>%
  group_map(~ {
    test_df <- mod_scores %>%
      filter(module == mod, visit == mod_visit) %>%
      inner_join(.x, by = "animalid") %>%
      drop_na(value, count)
    
    if (nrow(test_df) > 0) {
      rho <- cor(test_df$count, test_df$value, method = "spearman")
      pvalue <- cor.test(test_df$count, test_df$value, method = "spearman")$p.value
      tibble(
        module = mod,
        mod_visit = mod_visit,
        imm_key = unique(.x$key),
        n = nrow(test_df),
        rho = rho,
        pvalue = pvalue
      )
    } else {
      tibble(
        module = mod,
        mod_visit = mod_visit,
        imm_key = unique(.x$key),
        n = 0,
        rho = NA,
        pvalue = NA
      )
    }
  }) %>%
  bind_rows()

# Adjust p-values for multiple testing (FDR correction)
results <- results %>%
  mutate(fdrq = p.adjust(pvalue, method = "fdr"))

# Sort results by p-value and display the top 20
results <- results %>%
  arrange(pvalue) %>%
  slice_head(n = 20)

# View the results
print(results)