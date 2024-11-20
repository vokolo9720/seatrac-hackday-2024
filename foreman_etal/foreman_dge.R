library(readr)
library(tidyverse)
library(edgeR)
library(kimma)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load count data and meta-data
# Select the data for modeling (groups and visits)
# Create DGE object and set up the model
# Run the model and volcano plot
# Does model need adjustment for covariates? e.g., sex, route, dose


# Extensions
# Identify genes that are associated with low-birden granulomas in CD8 and/or CD4 T cells
# Do the replicate in Bromley et al. pseudo bulk datasets for CD4 and CD8 populations?

data_folder <- "/fh/fast/gilbert_p/fg_data/SEATRAC/TB_hackday_2024/data/foreman_etal/"
data_folder <- "foreman_etal/"

# These are already log2-CPM
raw <- readr::read_csv(paste(data_folder, "foreman_etal_counts.csv", sep=""))
meta <- readr::read_csv(paste(data_folder, "foreman_etal_meta.csv", sep=""))
# rownames(meta) = meta$sampleid


meta_ss = meta %>% filter(condition == "CD8_gran")
keep_ids = meta_ss %>% pull(sampleid)
keep_ids = c('gene_id', keep_ids)

ncts_ss = raw %>% select(any_of(keep_ids))

# Discard genes that have low counts/prevalence
filter = rowSums(ncts_ss > 1) >= (0.5 * ncol(ncts_ss))
ncts_ss = ncts_ss[filter, ]

# Create the object for differential expression testing
dge_o = DGEList(counts=ncts_ss,
                genes=ncts_ss[, 1],
                samples=meta_ss,
                group=meta_ss[['CFU']])

# Compute weights
ncts_ss <- calcNormFactors(dge_o, method = "TMM")

# Specify the model/design matrix
design_temp = model.matrix(~CFU, data=meta_ss)

# Create the voom object and fit the model
v <- voom(dge_o, design=design_temp, plot=TRUE)

# vvwts <- voomWithQualityWeights(dge_o, design=design_temp, normalize.method="none", plot=TRUE)
fit = lmFit(v, design_temp)

# Estimate contrasts and p-values
fit = eBayes(fit, robust=TRUE)

summary(decideTests(fit, adjust.method="fdr", p.value = 0.05))

topTable(fit, adjust="BH", resort.by="P", coef="CFU")

results <- topTable(fit, adjust="BH", coef="CFU", p.value=1, number=Inf, resort.by="P")

# Create a volcano plot for single-gene association with protection

# Add a column for significance based on FDR
results <- results %>%
  mutate(Significance = ifelse(adj.P.Val < 0.05, "Significant", "Not Significant"))

# Select the top 10 genes based on adjusted p-value for labeling
top_genes <- results %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 10)

max_logFC <- max(abs(results$logFC), na.rm = TRUE)

# Create the volcano plot
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Significance), alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_text_repel(data = top_genes,
                  aes(label = gene_id),
                  max.overlaps = 10,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  size = 3) +
  xlim(c(-max_logFC, max_logFC)) +
  theme_minimal() +
  labs(
    x = "log2 Fold-change",
    y = "-log10 P-value",
    color = "FDR < 0.05") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
print(volcano_plot)


# Redo the analysis using a mixed-effects model to account for the longitudinal design

dge_o = DGEList(counts=ncts_ss,
                genes=ncts_ss[, 1],
                samples=meta_ss,
                group=meta_ss[['CFU']])

ncts_ss <- calcNormFactors(dge_o, method = "TMM")

design_temp=model.matrix(~CFU, data=meta_ss)

v <- voom(dge_o, design=design_temp, plot=FALSE)

# Can't figure out why I get this error here
# lme/lmerel model: expression~protect_outcome+sex+visit+(1|animalid)
# Error in `[.data.frame`(weights.format, order(rownames(weights.format)),  : 
#                           undefined columns selected

klm <- kmFit(dat = v,
             model = "~CFU + (1|subject)",
             run_lme = TRUE,
             libraryID="sampleid",
             patientID="subject",
             use_weights = TRUE,
             metrics = TRUE,
             run_contrast = FALSE,
             processors=1)

summarise_kmFit(fdr = klm$lm)

plot_volcano(model_result = klm, 
             model = "lme", variables = "CFU",
             y_cutoff = 0.05)