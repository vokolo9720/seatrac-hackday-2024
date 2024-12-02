# Load required libraries
library(dplyr)       # For data manipulation
library(ggplot2)     # For data visualization
library(readr)       # For reading CSV/TSV files
library(tidyr)       # For reshaping data
library(stats)       # For basic statistical tests
library(multcomp)    # For multiple comparison corrections
library(lme4)        # For mixed-effects modeling
library(broom)       # For tidying model outputs

# Define file paths
git_repo <- "C:/Andrew/gitrepo"
data_folder <- file.path(git_repo, "seatrac-hackday-2024", "foreman_etal")

# Section 1: Load and preprocess data
# -----------------------------------
# Read TPM and metadata files
tpm <- read_tsv(file.path(data_folder, "ge_tpm.tsv"))
edds <- read_tsv(file.path(data_folder, "edds.tsv"))

# Select and rename important columns
keep_cols <- c("ColumnName", "biosample", "biosample_name", "subject", "compartment", 
               "sort", "gran", "Granuloma__", "CD4_T_cells", "CD8_T_cells", "CFU",
               "CD11b_Mtb_", "CD11b_Mtb_0", "condition")

edds <- edds %>%
  select(all_of(keep_cols)) %>%
  rename(
    Granuloma_num = Granuloma__,
    sampleid = ColumnName
  ) %>%
  mutate(logCFU = log10(CFU + 1))

# Save preprocessed metadata
write_csv(edds, file.path(data_folder, "foreman_etal_meta.csv"))

# Section 2: Identify sex from Y-chromosome genes
# -----------------------------------------------
ygenes <- c("DDX3Y", "UTY", "ZFY")

# Extract Y-chromosome gene expression and calculate sex
ytpm <- tpm %>%
  filter(gene_id %in% ygenes) %>%
  select(gene_id, all_of(edds$sampleid)) %>%
  column_to_rownames("gene_id")

sex <- colSums(ytpm) %>%
  log10() %>%
  as.data.frame() %>%
  setNames("value") %>%
  rownames_to_column("sampleid") %>%
  mutate(sex = if_else(value < 0.5, "M", "F")) %>%
  select(sampleid, sex)

# Merge sex information with metadata
edds <- edds %>%
  left_join(sex, by = "sampleid")

# Save updated metadata
write_csv(edds, file.path(data_folder, "foreman_etal_meta.csv"))

# Section 3: Log-transform TPM counts
# ------------------------------------
tpm <- tpm %>%
  select(gene_id, all_of(edds$sampleid)) %>%
  mutate(across(-gene_id, ~ log2(. + 1)))

write_csv(tpm, file.path(data_folder, "foreman_etal_counts.csv"))

# Section 4: Filter genes for analysis
# ------------------------------------
ltpm <- tpm %>%
  column_to_rownames("gene_id") %>%
  filter(rowMeans(. > 1) > 0.6)

genes <- rownames(ltpm)

# Merge metadata and log-transformed counts
meta <- edds %>%
  filter(condition %in% c("CD4_gran", "CD4_PBMC")) %>%
  column_to_rownames("sampleid")

mdf <- meta %>%
  bind_cols(t(ltpm)) %>%
  rownames_to_column("sampleid")

# Section 5: Mann-Whitney U Test for Gene Associations
# ----------------------------------------------------
results <- lapply(genes, function(gene) {
  gran_expr <- mdf %>%
    filter(condition == "CD4_gran") %>%
    pull(gene)
  pbmc_expr <- mdf %>%
    filter(condition == "CD4_PBMC") %>%
    pull(gene)
  
  test <- wilcox.test(gran_expr, pbmc_expr)
  mean_diff <- mean(gran_expr) - mean(pbmc_expr)
  
  list(
    gene = gene,
    pvalue = test$p.value,
    stat = test$statistic,
    assoc = ifelse(mean_diff > 0, "GRAN", "PBMC")
  )
})

results_df <- bind_rows(results) %>%
  arrange(pvalue) %>%
  mutate(FDRq = p.adjust(pvalue, method = "fdr"))

print(head(results_df, 20))

# Section 6: Correlation with CFU Using OLS Models
# ------------------------------------------------
sig_genes <- results_df %>%
  filter(FDRq < 0.1) %>%
  pull(gene)

cfu_res <- lapply(sig_genes, function(gene) {
  formula <- as.formula(paste("logCFU ~", gene))
  model <- lm(formula, data = mdf)
  summary <- summary(model)
  
  conf_int <- confint(model)
  rho <- cor(mdf$logCFU, mdf[[gene]], method = "spearman")
  pvalue_corr <- cor.test(mdf$logCFU, mdf[[gene]], method = "spearman")$p.value
  
  list(
    gene = gene,
    est = summary$coefficients[gene, "Estimate"],
    est_lb = conf_int[gene, 1],
    est_ub = conf_int[gene, 2],
    pvalue = summary$coefficients[gene, "Pr(>|t|)"],
    rho = rho,
    corr_pvalue = pvalue_corr
  )
})

cfu_res_df <- bind_rows(cfu_res) %>%
  arrange(pvalue) %>%
  mutate(FDRq = p.adjust(pvalue, method = "fdr"))

print(head(cfu_res_df, 20))

# Section 7: Mixed-Effects Models
# --------------------------------
gene <- sig_genes[1] # Example with one gene
formula <- as.formula(paste("logCFU ~", gene, "+ condition"))

mixed_model <- lmer(
  formula,
  data = mdf,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"),
  random = ~ 1 | subject
)

summary(mixed_model)

# Section 8: Visualization
# ------------------------
# Scatterplot of Gene Expression vs. CFU
ggplot(mdf, aes_string(x = "logCFU", y = gene)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  labs(
    title = paste("Correlation of", gene, "with logCFU"),
    x = "logCFU",
    y = paste("Expression of", gene)
  ) +
  theme_minimal()

# Boxplot of Gene Expression by Condition
ggplot(mdf, aes_string(x = "condition", y = gene)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(
    title = paste("Expression of", gene, "by Condition"),
    x = "Condition",
    y = paste("Expression of", gene)
  ) +
  theme_minimal()
