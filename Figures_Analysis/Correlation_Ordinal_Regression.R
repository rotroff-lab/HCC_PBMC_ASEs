library(dplyr)
library(pheatmap)
library(tidyr)
library(MASS)
library(caret)

# Load metadata
metadata <- readRDS("metadata.txt")

# Load and annotate CIBERSORT data
CIBERSORT <- read.csv("CIBERSORT.csv")
cell_types <- names(CIBERSORT)[-1]

# Load PSI data
psi <- readRDS("PSI.rds")
psi <- as.data.frame(t(psi[-1,]))


features <- readRDS("immune_coding_ASEs.rds")
psi <- psi %>% select(RUN, all_of(features))

all_data <- left_join(CIBERSORT, psi)
all_data[, -1] <- lapply(all_data[, -1], as.numeric)

# Remove low variance ASEs and Cell types
low_var_cells <- nearZeroVar(all_data[cell_types], saveMetrics = TRUE)
cell_types <- rownames(low_var_cells[low_var_cells$nzv == FALSE,])
low_var_ASEs <- nearZeroVar(all_data[features], saveMetrics = TRUE)
ASEs <- rownames(low_var_ASEs[low_var_ASEs$nzv == FALSE,])

# Compute correlations
cor_results <- list()

# Loop through each ASE
for (ase in features) {
  for (cell in cell_types) {
    
    # Compute correlation
    cor_test <- cor.test(all_data[[ase]], all_data[[cell]], method = "spearman")
    
    # Store results
    cor_results <- append(cor_results, list(data.frame(
      ASE = ase,
      Cell_Type = cell,
      Correlation = cor_test$estimate,
      P_Value = cor_test$p.value
    )))
  }
}

cor_results_df <- do.call(rbind, cor_results)
cor_results_df$FDR <- p.adjust(cor_results_df$P_Value, method = "fdr")
cor_sig <- cor_results_df  %>% filter(FDR < 0.05)
write.csv(cor_results_df, "Table_S5.csv", row.names = FALSE)



##Ordinal Regressions
# Recode Staging for Ordinal Regression

metadata <- metadata %>%
  mutate(BCLC_STAGE = case_when(
    BCLC_STAGE %in% c("A0", "I", "Ib") ~ "A",
    BCLC_STAGE %in% c("II", "IIb") ~ "B",
    BCLC_STAGE == "III" ~ "C",
    BCLC_STAGE %in% c("IV", "IVbb") ~ "D",
    TRUE ~ BCLC_STAGE
  ))

# Merge ASE PSI with metadata
psi <- left_join(psi, metadata[, c("RUN", "BCLC_STAGE")], by = "RUN") %>%
  filter(!is.na(BCLC_STAGE) & BCLC_STAGE != "U")

psi$BCLC_STAGE <- factor(psi$BCLC_STAGE, levels = c("A", "B", "C", "D"), ordered = TRUE)

# Merge cell type data with metadata
ciber <- left_join(CIBERSORT, metadata[, c("RUN", "BCLC_STAGE")], by = "RUN") %>%
  filter(!is.na(BCLC_STAGE) & BCLC_STAGE != "U")

ciber$BCLC_STAGE <- factor(ciber$BCLC_STAGE, levels = c("A", "B", "C", "D"), ordered = TRUE)


# Ordinal logistic regression for ASEs
result <- data.frame()
for (ase in features) {
  tryCatch({
    model <- polr(as.formula(paste("BCLC_STAGE ~", ase)), data = psi2, method = "logistic")
    t_value <- coef(model)[1] / sqrt(diag(vcov(model)))[1]
    p_value <- 2 * (1 - pnorm(abs(t_value)))
    result <- rbind(result, data.frame(ASE = ase, Coefficients = coef(model), T_values = t_value, P_values = p_value))
  }, error = function(e) {
    message("Error occurred for ", ase, ", skipping.")
  })
}
result$FDR <- p.adjust(result$P_values, method = "fdr")
write.csv(result, "Table_S7.csv", row.names = FALSE)


# Ordinal logistic regression for cell types
result_cell <- data.frame()
for (cell in cell_types) {
  tryCatch({
    model <- polr(as.formula(paste("BCLC_STAGE ~", cell)), data = ciber, method = "logistic")
    t_value <- coef(model)[1] / sqrt(diag(vcov(model)))[1]
    p_value <- 2 * (1 - pnorm(abs(t_value)))
    result_cell <- rbind(result_cell, data.frame(Cell_Type = cell, Coefficients = coef(model), T_values = t_value, P_values = p_value))
  }, error = function(e) {
    message("Error occurred for ", cell, ", skipping.")
  })
}
result_cell$FDR <- p.adjust(result_cell$P_values, method = "fdr")
write.csv(result_cell, "Table_S8.csv", row.names = FALSE)
