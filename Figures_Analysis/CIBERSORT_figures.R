library(tidyverse)
library(umap)
library(ggplot2)
library(pheatmap)
library(openxlsx)


metadata <- readRDS("metadata.rds")

# Join CIBERSORT data with metadata and filter out NA groups
CIBERSORT  <- read.csv("CIBERSORT.csv")
CIBERSORT_annotated <- left_join(CIBERSORT, metadata[, c("RUN", "SRA.STUDY", "GROUP")], by = "RUN")

# Filter CIBERSORT_UMAP based on annotated RUNs
CIBERSORT_UMAP <- filter(CIBERSORT, RUN %in% CIBERSORT_annotated$RUN) %>%
  select(-Correlation, -P.value, -RMSE) %>%
  column_to_rownames(var = "RUN")

# Perform UMAP
PSI_UMAP <- umap(CIBERSORT_UMAP, n_components = 2, random_state = 23, n_neighbors = 50) 
layout <- data.frame(PSI_UMAP[["layout"]])
layout$RUN <- rownames(layout)
df_plot <- left_join(layout,CIBERSORT_annotated, by = "RUN")
df_plot$GROUP <- factor(df_plot$GROUP, levels = c("Healthy", "HBV", "HCC"))

# Generate UMAP plots
plot_umap <- ggplot(df_plot, aes(x = X1, y = X2, color = GROUP)) +
  geom_point(size = 2, alpha = 1) +
  scale_color_manual(values = c("lightblue", "#097969", "#880808")) +
  facet_wrap(~SRA.STUDY, ncol = 3) +
  theme_classic()
ggsave("CIBERSORT_SRA_startified.pdf", plot_umap, width = 7, height = 5, units = "in", dpi = 1000)

plot_umap2 <- ggplot(df_plot, aes(x=X1, y=X2,color=df_plot$GROUP)) +
  geom_point(size=2, alpha = 1) +
  scale_color_manual(values = c("lightblue", "#097969", "#880808"))+
  theme_classic()
ggsave("CIBERSORT_SRA_combined.pdf", plot_umap2, width = 7, height = 5, units = "in", dpi = 1000)

# Create custom color palettes
my_colour <- list(
  GROUP = c(Healthy = "lightblue", HCC = "#880808", HBV = "#097969"),
  SRA.STUDY = c(SRP094502 = "#127734", SRP162958 = "gold", SRP312121 = "grey", SRP332585 = "darkred", SRP412160 = "#2a9df4", SRP330954 = "#ff5f00", SRP201023 = "black")
)

custom_palette <- c("#f2f3f2", colorRampPalette(c("#ADD8E6", "darkred"))(10000))

# Create pheatmap
pheatmap(CIBERSORT_UMAP, 
         annotation_row = as.data.frame(CIBERSORT_annotated[, c("GROUP", "SRA.STUDY")]),
         cluster_cols = TRUE,
         annotation_colors = my_colour,
         show_rownames = FALSE,
         color = custom_palette,
         border_color = NA,
         clustering_method = "ward.D2",
         filename = "figures/CIBERSORT_heatmap.pdf")

# Prepare data for stacked bar plot
stacked_bar <- left_join(CIBERSORT_UMAP, CIBERSORT_annotated[, c("RUN", "GROUP")], by = "rowname") %>%
  mutate(test = paste0(GROUP, rownames(.))) %>%
  pivot_longer(cols = -c(rowname, GROUP, SRA.STUDY, test), names_to = "CellType", values_to = "Proportion") %>%
  filter(CellType != "Correlation" & CellType != "P.value" & CellType != "RMSE")

# Generate stacked bar plot
plot_stacked_bar <- ggplot(stacked_bar, aes(x = test, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y = "Proportion", x = "Sample", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors_22) +
  theme_classic()

ggsave("figures/cibersort_stackedbar.pdf", plot_stacked_bar, width = 10, height = 6)

# Prepare table data
CIBERSORT_table <- left_join(CIBERSORT, CIBERSORT_annotated[, c("RUN", "GROUP", "SRA.STUDY")], by = "RUN") %>%
  select(RUN, SRA.STUDY, GROUP, everything()) %>%
  mutate(across(.cols = -c(1:3), .fns = ~round(.x, 3)))

write.xlsx(CIBERSORT_table, "figures/Table_S4.xlsx")

# Logistic regression function
glmfit <- function(df, cell_types) {
  result_list <- list()
  for (cell_type in cell_types) {
    df_subset <- df[, c(cell_type, "GROUP")]
    for (i in 0:1) {
      for (j in (i + 1):2) {
        subset_ij <- df_subset[df_subset$GROUP %in% c(i, j), ]
        fit <- glm(GROUP ~ ., data = subset_ij, family = "binomial")
        model_summary <- summary(fit)
        p_value <- round(model_summary$coefficients[2, 4], digits = 10)
        beta <- round(model_summary$coefficients[2, 1], digits = 3)
        std_err <- round(model_summary$coefficients[2, 2], digits = 3)
        result <- data.frame(
          Group1 = paste0("Group ", i),
          Group2 = paste0("Group ", j),
          p_value = p_value,
          coefficient = beta,
          std_err = std_err,
          Cell_type = cell_type,
          max = max(subset_ij[, cell_type])
        )
        result_list[[paste0(cell_type, "_Group", i, "_vs_Group", j)]] <- result
      }
    }
  }
  final_result <- do.call(rbind, result_list)
  return(final_result)
}

# Perform logistic regression and save results
result <- glmfit(CIBERSORT_annotated, cell_types)
result$fdr <- p.adjust(result$p_value, method = "fdr")
result$Group1 <- factor(ifelse(result$Group1 == "Group 0", "HCC", ifelse(result$Group1 == "Group 1", "Healthy", "HBV")))
result$Group2 <- factor(ifelse(result$Group2 == "Group 0", "HCC", ifelse(result$Group2 == "Group 1", "Healthy", "HBV")))
result <- result %>% filter(fdr <= 0.05)
result$fdr <- format(result$fdr, scientific = TRUE, digits = 2)

# Boxplot of WBC counts per group
p <- ggplot(CIBERSORT_annotated_SRA, aes(x = GROUP, y = Proportion, color = GROUP)) +
  geom_boxplot(position = position_dodge(), outlier.size = 0.5) +
  geom_dotplot(aes(fill = GROUP), binaxis = 'y', stackdir = 'center', position = position_dodge(width = 0.75), dotsize = 0.5) +
  facet_wrap(~Cell_type, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("lightblue", "#097969", "#880808")) +
  scale_color_manual(values = c("lightblue", "#097969", "#880808")) +
  theme_classic()

p + stat_pvalue_manual(result, label = "fdr", y.position = "max", tip.length = 0.015, hide.ns = TRUE, bracket.size = 0.5, label.size = 3)
ggsave("CIBERSORT_glm.pdf", width = 8, height = 11, units = "in", dpi = 600)
