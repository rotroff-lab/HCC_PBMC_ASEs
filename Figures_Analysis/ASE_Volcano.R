library(dplyr)
library(ggplot2)
library(readxl)

# Load results
results <- readRDS("meta_analysis_ASE.rds") %>% filter(annotated == "annotated")

# Define significance categories
results <- results %>% mutate(
  sig = case_when(
    abs(meanincleveldiff) >= 0.1 & fdr <= 0.05 ~ "dPSI > 0.1",
    abs(meanincleveldiff) >= 0.05 & fdr <= 0.05 ~ "dPSI > 0.05",
    TRUE ~ "Not Sig"
  )
)

# Define colors
sig_colors <- c("Not Sig" = "grey", "dPSI > 0.05" = "#FF2400", "dPSI > 0.1" = "#880808")

# Volcano plot function
plot_volcano <- function(data, output, legend = TRUE) {
  p <- ggplot(data, aes(x = meanincleveldiff, y = -log10(fdr), color = sig)) +
    geom_point(size = 1) +
    xlim(-0.35, 0.35) +
    scale_color_manual(values = sig_colors) +
    geom_hline(yintercept = 1.30103, linetype = "dashed", color = "#FF2400") +
    theme_bw()
  
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  
  ggsave(output, plot = p, device = "pdf", units = "in", height = 3, width = 6)
}

# Generate volcano plots
plot_volcano(results, "ASE_Volcano.pdf")
plot_volcano(results, "ASE_Volcano_no_legend.pdf", legend = FALSE)

# Load immune gene list
gene_list <- read_excel("data/GO0002376_gene_names.xlsx", col_names = FALSE) %>%
  rename(gene_name = ...1) %>%
  mutate(gene_type = "immune")

# Merge immune annotations
results_immune <- left_join(results, gene_list, by = "gene_name") %>%
  mutate(gene_type = ifelse(is.na(gene_type), "other", "immune"))

# Define immune-specific significance categories
results_immune <- results_immune %>% mutate(
  sig = case_when(
    gene_type == "immune" & abs(meanincleveldiff) >= 0.1 & fdr <= 0.05 ~ "Immune dPSI > 0.1",
    gene_type == "immune" & abs(meanincleveldiff) >= 0.05 & fdr <= 0.05 ~ "Immune dPSI > 0.05",
    gene_type == "immune" ~ "Immune Not Sig",
    TRUE ~ "Other"
  )
)

# Define immune colors
immune_colors <- c("Other" = "#eeeeee", "Immune Not Sig" = "lightblue", "Immune dPSI > 0.05" = "steelblue", "Immune dPSI > 0.1" = "darkblue")

# Generate immune volcano plot
plot_volcano_immune <- function(data, output, legend = TRUE) {
  p <- ggplot(data, aes(x = meanincleveldiff, y = -log10(fdr), color = sig)) +
    geom_point(size = 1) +
    xlim(-0.35, 0.35) +
    scale_color_manual(values = immune_colors) +
    geom_hline(yintercept = 1.30103, linetype = "dashed", color = "steelblue") +
    theme_bw()
  
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  
  ggsave(output, plot = p, device = "pdf", units = "in", height = 3, width = 6)
}

plot_volcano_immune(results_immune, "ASE_Volcano_Immune.pdf")
plot_volcano_immune(results_immune, "ASE_Volcano_Immune_no_legend.pdf", legend = FALSE)

# Summary table for ASE significance
summary_all <- results %>%
  filter(abs(meanincleveldiff) >= 0.05 & fdr <= 0.05) %>%
  mutate(Direction = ifelse(meanincleveldiff > 0, "up", "down")) %>%
  count(type, Direction)

summary_immune <- results_immune %>%
  filter(gene_type == "immune", abs(meanincleveldiff) >= 0.05 & fdr <= 0.05) %>%
  mutate(Direction = ifelse(meanincleveldiff > 0, "up", "down")) %>%
  count(type, Direction)
