library(ggplot2)
library(dplyr)
library(pROC)
library(tidyr)

# Load Gini Data
gini <- readRDS("gini.rds")
rownames(gini) <- gini$feature
gini$feature <- NULL
gini <- as.matrix(gini)
gini_plot <- data.frame(
  feature = rownames(gini), 
  min = rowMins(gini), 
  max = rowMaxs(gini), 
  median = rowMedians(gini)
) %>% arrange(desc(median))

gini_plot$feature <- factor(gini_plot$feature, levels = rev(gini_plot$feature))

# Plot Gini Scores
ggplot(gini_plot, aes(x = feature, y = median, ymin = min, ymax = max)) +
  geom_pointrange() +
  coord_flip() +
  xlab("ASE") + ylab("Median Gini Score") +
  theme_classic()
ggsave("gini_three_group.pdf", device = "pdf", height = 6, width = 3.5)

# Load Predictions
results_train <- readRDS("training_results.rds")
results_test <- readRDS("test_results.rds")
results_validation <- readRDS("validation_results.rds")

# Function to create two-group labels
convert_two_group_labels <- function(df, col_name, target_group, target_label, non_target_label) {
  df %>%
    mutate({{ col_name }} := ifelse({{ col_name }} == target_group, target_label, non_target_label))
}

# Generate datasets for each condition in train, test, and validation
generate_two_group_datasets <- function(df, col_name) {
  list(
    HCC = convert_two_group_labels(df, {{ col_name }}, "group1", "HCC", "NonHCC"),
    Healthy = convert_two_group_labels(df, {{ col_name }}, "group2", "Healthy", "NonHealthy"),
    HBV = convert_two_group_labels(df, {{ col_name }}, "group3", "HBV", "NonHBV")
  )
}

# Apply transformation for train, test, and validation sets
train_datasets <- generate_two_group_datasets(results_train, diagnosis)
test_datasets <- generate_two_group_datasets(results_test, truth)
validation_datasets <- generate_two_group_datasets(results_validation, truth)

# Function to compute ROC
compute_roc <- function(df, label_col, prob_col) {
  roc(df[[label_col]], df[[prob_col]])
}

# Compute ROC Curves for each class
roc_train <- list(
  HCC = compute_roc(train_datasets$HCC, "diagnosis", "prob.group1"),
  Healthy = compute_roc(train_datasets$Healthy, "diagnosis", "prob.group2"),
  HBV = compute_roc(train_datasets$HBV, "diagnosis", "prob.group3")
)

roc_test <- list(
  HCC = compute_roc(test_datasets$HCC, "truth", "prob.group1"),
  Healthy = compute_roc(test_datasets$Healthy, "truth", "prob.group2"),
  HBV = compute_roc(test_datasets$HBV, "truth", "prob.group3")
)

roc_validation <- list(
  HCC = compute_roc(validation_datasets$HCC, "truth", "prob.group1"),
  HBV = compute_roc(validation_datasets$HBV, "truth", "prob.group3")
)

# Print AUC Values
lapply(roc_test, auc)
lapply(roc_train, auc)
lapply(roc_validation, auc)

# Plot ROC Curves
plot_roc_curves <- function(roc_list, colors, filename) {
  pdf(filename, height = 4, width = 4)
  plot(roc_list[[1]], print.auc = TRUE, col = colors[1], print.auc.y = .45, print.auc.x = 0.3)
  for (i in 2:length(roc_list)) {
    plot(roc_list[[i]], print.auc = TRUE, col = colors[i], print.auc.y = .35 - (i - 2) * 0.1, print.auc.x = 0.3, add = TRUE)
  }
  dev.off()
}

plot_roc_curves(roc_train, c("#FF2400", "#AFE1AF", "lightblue"), "ROC_three_group_training.pdf")
plot_roc_curves(roc_test, c("#880808", "#097969", "darkblue"), "ROC_three_group_test.pdf")
plot_roc_curves(roc_validation, c("#340002", "#003400"), "ROC_three_group_validation.pdf")







##Metrics Bar Plot
results_train <- readRDS("training_results.rds")
results_test <- readRDS("test_results.rds")
results_validation <- readRDS("validation_results.rds")

# Rename Columns for Consistency
rename_columns <- function(df) {
  colnames(df)[1:2] <- c("pred", "diagnosis")
  return(df)
}

results_test <- rename_columns(results_test)
results_validation <- rename_columns(results_validation)

# Function to Compute Accuracy Metrics
compute_metrics <- function(group, df) {
  TP <- sum(df$pred == group & df$diagnosis == group)
  TN <- sum(df$pred != group & df$diagnosis != group)
  FP <- sum(df$pred == group & df$diagnosis != group)
  FN <- sum(df$pred != group & df$diagnosis == group)
  total <- nrow(df)
  
  metrics <- data.frame(
    accuracy = (TP + TN) / total,
    Balanced_Accuracy = ((TP / (TP + FN)) + (TN / (FP + TN))) / 2,
    Misclassification = (FP + FN) / total,
    Sensitivity = TP / (TP + FN),
    Specificity = TN / (TN + FP),
    PPV = TP / (TP + FP),
    NPV = TN / (TN + FP),
    stringsAsFactors = FALSE
  )
  
  # Transpose and format
  metrics <- as.data.frame(t(metrics))
  colnames(metrics) <- group
  return(round(metrics * 100, 2)) 
}

# Compute Metrics for Each Dataset
groups <- c("group1", "group2", "group3")
results_train_sum <- do.call(cbind, lapply(groups, compute_metrics, df = results_train))
results_test_sum <- do.call(cbind, lapply(groups, compute_metrics, df = results_test))
results_validation_sum <- do.call(cbind, lapply(groups, compute_metrics, df = results_validation))

# Convert to Data Frame and Add Group Labels
format_results <- function(results, group_name) {
  results <- as.data.frame(results)
  results$group <- group_name
  results$Metric <- rownames(results)
  return(results)
}

results_train_sum <- format_results(results_train_sum, "train")
results_test_sum <- format_results(results_test_sum, "test")
results_validation_sum <- format_results(results_validation_sum, "validation")

# Combine All Data
test_train_class <- bind_rows(results_train_sum, results_test_sum, results_validation_sum)
colnames(test_train_class) <- c("HCC", "Healthy", "HBV", "group", "Metric")

# Reshape Data for Plotting
test_train_class <- test_train_class %>%
  pivot_longer(cols = c(HCC, Healthy, HBV), names_to = "source", values_to = "combined") %>%
  mutate(source = factor(source, levels = c("HCC", "HBV", "Healthy"))) %>%
  filter(!(group == "validation" & source == "Healthy"))

# Define Group Labels for Faceting
test_train_class$group <- case_when(
  test_train_class$group == "test" & test_train_class$source == "HCC" ~ "test_HCC",
  test_train_class$group == "train" & test_train_class$source == "HCC" ~ "train_HCC",
  test_train_class$group == "test" & test_train_class$source == "Healthy" ~ "test_Healthy",
  test_train_class$group == "train" & test_train_class$source == "Healthy" ~ "train_Healthy",
  test_train_class$group == "test" & test_train_class$source == "HBV" ~ "test_HBV",
  test_train_class$group == "validation" & test_train_class$source == "HCC" ~ "validation_HCC",
  test_train_class$group == "validation" & test_train_class$source == "HBV" ~ "validation_HBV",
  test_train_class$group == "train" & test_train_class$source == "HBV" ~ "train_HBV",
  TRUE ~ test_train_class$group
)

test_train_class$group <- factor(
  test_train_class$group,
  levels = c("train_HCC", "test_HCC", "validation_HCC", "train_HBV", "test_HBV", "validation_HBV", "train_Healthy", "test_Healthy")
)

# Remove Accuracy Metric from Final Plot
test_train_class <- test_train_class %>% filter(Metric != "accuracy")

# Format Metric Names for Readability
test_train_class$Metric <- recode(test_train_class$Metric, "Balanced_Accuracy" = "Balanced Accuracy")
test_train_class$Metric <- factor(test_train_class$Metric, levels = c("Balanced Accuracy", "Sensitivity", "Specificity", "Misclassification", "PPV", "NPV"))

# Plot Results
ggplot(test_train_class, aes(source, combined, fill = group)) +
  geom_col(position = position_dodge(), width = 0.5) +
  ylim(0, 100) +
  facet_wrap(~Metric) +
  scale_fill_manual(values = c("#FF2400", "#880808", "#340002", "#AFE1AF", "#097969", "#003400", "lightblue", "darkblue")) +
  theme_classic()
ggsave("all_metrics_bar_plot_three_group.pdf", device = "pdf", height = 4, width = 8)


# Waffle Plots
source("waffle_plot_three_group.R")
pdf("waffle_plot_two_group.pdf", height = 8, width = 18)
print(p)
dev.off()


# Function to Compute Prediction Proportions
compute_proportions <- function(df, group_label) {
  df %>%
    count(diagnosis, pred) %>%
    mutate(
      pred = recode(pred, "group1" = "HCC", "group2" = "Healthy", "group3" = "HBV"),
      diagnosis = recode(diagnosis, "group1" = "HCC", "group2" = "Healthy", "group3" = "HBV")
    ) %>%
    group_by(pred) %>%
    mutate(prop_pred = (n / sum(n)) * 100) %>%
    ungroup() %>%
    mutate(group = group_label)
}

# Compute Proportions for All Datasets
prop_all <- bind_rows(
  compute_proportions(results_train, "train"),
  compute_proportions(results_test, "test"),
  compute_proportions(results_validation, "validation")
) %>%
  mutate(group = factor(group, levels = c("train", "test", "validation")))

# Plot
ggplot(prop_all, aes(x = diagnosis, y = prop_pred, fill = diagnosis)) +
  geom_bar(stat = "identity") +
  facet_wrap(group ~ pred, nrow = 1) +
  ylim(0, 100) +
  theme_classic() +
  scale_fill_manual(name = NULL, values = diagnosis_colors) +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    panel.spacing = unit(2, "lines")
  )

ggsave("waffle_bar_plot_three_group.pdf", 
       device = "pdf", units = "in", height = 2, width = 14)
