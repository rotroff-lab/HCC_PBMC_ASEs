# Load required libraries
library(ggplot2)
library(pROC)
library(dplyr)

# Load Gini importance scores
gini <- readRDS("gini.rds")
rownames(gini) <- gini$feature
gini$feature <- NULL
gini <- as.matrix(gini)

# Compute summary statistics for Gini scores
gini_plot <- data.frame(
  feature = rownames(gini),
  min = rowMins(gini),
  max = rowMaxs(gini),
  median = rowMedians(gini)
)

# Reverse factor level ordering for labels after coord_flip()
gini_plot <- gini_plot[order(gini_plot$median, decreasing = TRUE), ]
gini_plot$feature <- factor(gini_plot$feature, levels = rev(gini_plot$feature))

# Plot Gini scores
ggplot(gini_plot, aes(x = feature, y = median, ymin = min, ymax = max)) +
  geom_pointrange() +
  coord_flip() +
  xlab("ASE") +
  ylab("Median Gini Score") +
  theme_classic()
ggsave("gini_two_group.pdf", 
       device = "pdf", units = "in", height = 6, width = 3.5)

# Load classification results
results_train <- readRDS("training_results.rds")
results_test <- readRDS("test_results.rds")
results_validation <- readRDS("validation_results.rds")

# Compute ROC curves
roc_hcc_train <- roc(results_train$diagnosis, results_train$prob.group1)
roc_hcc_test <- roc(results_test$truth, results_test$prob.group1)
roc_hcc_validation <- roc(results_validation$truth, results_validation$prob.group1)

# Print AUC values
cat("AUC values:\n")
cat("Train AUC:", auc(roc_hcc_train), "\n")
cat("Test AUC:", auc(roc_hcc_test), "\n")
cat("Validation AUC:", auc(roc_hcc_validation), "\n")

# Plot ROC curves
pdf("ROC_two_group.pdf", height = 4, width = 4)
plot(roc_hcc_train, print.auc = TRUE, col = "lightblue")
plot(roc_hcc_test, print.auc = TRUE, col = "darkblue", print.auc.y = 0.4, add = TRUE)
plot(roc_hcc_validation, print.auc = TRUE, col = "red", print.auc.y = 0.3, add = TRUE)
dev.off()

# Assign splicing biomarker classifications
assign_classifications <- function(df) {
  df$splicingbiomarker <- NA
  df$splicingbiomarker[df$prediction == "group1" & df$truth == "group1"] <- "TP"
  df$splicingbiomarker[df$prediction == "group2" & df$truth == "group2"] <- "TN"
  df$splicingbiomarker[df$prediction == "group2" & df$truth == "group1"] <- "FN"
  df$splicingbiomarker[df$prediction == "group1" & df$truth == "group2"] <- "FP"
  return(df)
}

results_train <- assign_classifications(results_train)
results_test <- assign_classifications(results_test)
results_validation <- assign_classifications(results_validation)

# Function to calculate classification metrics
calculate_metrics <- function(df) {
  TP <- sum(df$splicingbiomarker == "TP")
  TN <- sum(df$splicingbiomarker == "TN")
  FP <- sum(df$splicingbiomarker == "FP")
  FN <- sum(df$splicingbiomarker == "FN")
  
  balanced_accuracy <- (TP / (TP + FN) + TN / (TN + FP)) / 2
  misclassification <- (FP + FN) / nrow(df)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  
  data.frame(
    Metric = c("Balanced Accuracy", "Misclassification", "Sensitivity", 
               "Specificity", "PPV", "NPV"),
    Value = c(balanced_accuracy, misclassification, sensitivity, 
              specificity, PPV, NPV)
  )
}

# Compute and aggregate metrics
results_train_sum <- calculate_metrics(results_train)
results_test_sum <- calculate_metrics(results_test)
results_validation_sum <- calculate_metrics(results_validation)

results_train_sum$group <- "train"
results_test_sum$group <- "test"
results_validation_sum$group <- "validation"

# Combine all results
test_train_class <- rbind(results_train_sum, results_test_sum, results_validation_sum)
test_train_class$group <- factor(test_train_class$group, levels = c("train", "test", "validation"))
test_train_class$Metric <- factor(test_train_class$Metric, levels = c("Balanced Accuracy", "Sensitivity", 
                                                                      "Specificity", "Misclassification", "PPV", "NPV"))
test_train_class$Value <- test_train_class$Value * 100

# Plot classification metrics
ggplot(test_train_class, aes(Metric, Value, fill = group)) +
  geom_col(position = position_dodge(), width = 0.5) +
  scale_fill_manual(values = c("lightblue", "darkblue", "red")) +
  ylim(0, 100) +
  theme_classic()
ggsave("all_metrics_bar_plot.pdf", 
       device = "pdf", units = "in", height = 4, width = 8)

# Generate waffle plots
source("waffle_plot_two_group.R")

pdf("waffle_plot_two_group.pdf", height = 6, width = 16)
print(p)
dev.off()


process_results <- function(data, group_label) {
  data %>%
    count(model_6, pred) %>%
    mutate(
      pred = recode(pred, "group1" = "HCC", "group2" = "Healthy"),
      diagnosis = recode(model_6, "group1" = "HCC", "group2" = "Healthy", "group3" = "HBV"),
      group = group_label
    ) %>%
    group_by(pred) %>%
    mutate(prop_pred = (n / sum(n)) * 100) %>%
    ungroup()
}

# Process training, test, and validation sets
prop_train <- process_results(results_train, "train")
prop_test <- process_results(results_test, "test")
prop_validation <- process_results(results_validation, "validation")

# Combine all data
prop_all <- bind_rows(prop_train, prop_test, prop_validation) %>%
  mutate(
    group = factor(group, levels = c("train", "test", "validation")),
    color = diagnosis,
    diagnosis = ifelse(diagnosis == "HBV", "Healthy", diagnosis)
  )

# Define colors
diagnosis_colors <- c("HCC" = "#880808", "Healthy" = "lightblue", "HBV" = "#097969")

# Plot
ggplot(prop_all, aes(x = diagnosis, y = prop_pred, fill = color)) +
  geom_bar(stat = "identity", width = 0.75) +
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

# Save plot
ggsave("waffle_bar_plot_two_group.pdf", 
       device = "pdf", units = "in", height = 2, width = 12)

