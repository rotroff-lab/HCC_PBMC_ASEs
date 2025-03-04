library(randomForest)
library(readxl)
library(dplyr)
library(matrixStats)
library(pROC)



################################################################################
#ASEs
################################################################################

# Import PSI table. Scale and Center data
PSI <- readRDS("PSI.rds")
PSI[, grep("^SRR", names(PSI))] <- as.data.frame(scale(PSI[, grep("^SRR", names(PSI))], center = TRUE))

# Subset for Immune ASEs
all_immune_ASEs_validation <- readRDS("all_immune_ASEs_validation.rds")
psi <- PSI %>% filter(AS_ID_short %in% all_immune_ASEs_validation) %>% dplyr::select(-ID,-gene_name, -type, -annotated)
psi <- t(psi)
psi <- psi[-1,]

# Import sample partitions
metadata_sample_partitions <- readRDS("partition.rds")
metadata_sample_partitions <- metadata_sample_partitions %>% select(RUN, three_group, model)
metadata_sample_partitions <- metadata_sample_partitions %>%
  dplyr::rename(diagnosis = three_group)
metadata_training <- metadata_sample_partitions %>% filter(model == "train")
metadata_test <- metadata_sample_partitions %>% filter(model == "test")
metadata_validation <- metadata_sample_partitions %>% filter(model == "validation")


# Add diagnosis data to psi table
diagnosis <- metadata_sample_partitions[,c("RUN","diagnosis")]
df_ASEs <- merge(diagnosis, psi, by.x="RUN", by.y="row.names")
df_ASEs$diagnosis <- as.factor(df_ASEs$diagnosis)
rownames(df_ASEs) <- df_ASEs$RUN
df_ASEs$RUN <- NULL

# Extract features
rf_classifier <- readRDS("RF_model.rds")
features_ASE <- rownames(rf_classifier$importance)
df_ASEs <- df[, c("diagnosis", features_ASE)]



################################################################################
#Cell_Types
################################################################################

# Cell Proportions
CIBERSORT  <- read.csv("CIBERSORT.csv")
names(CIBERSORT)[1] <- "RUN"
CIBERSORT <- CIBERSORT %>% dplyr::select(-Correlation, -P.value, -RMSE, -Dendritic.cells.resting)


#prepare df for RF
df_cell <- merge(diagnosis, CIBERSORT, by.x="RUN", by.y="RUN")
df_cell$diagnosis <- as.factor(df_cell$diagnosis)
rownames(df_cell) <- df_cell$RUN
df_cell$RUN <- NULL


df_cell[, -1] <- apply(df_cell[, -1], 2, as.numeric)
features_cell <- names(df_cell[-1])


################################################################################
#ASEs and Cell Types
################################################################################

df_all <- merge(df_ASEs, df_cell, by = c("row.names", "diagnosis"))
rownames(df_all) <- df_all$Row.names
df_all <- df_all[,-1]
df_all[, -1] <- apply(df_all[, -1], 2, as.numeric)
features_all <- names(df_all[-1])


## Functions for RF
run_random_forest <- function(train, test, features) {
  set.seed(10)
  train <- train[, c("diagnosis", features)]
  train[, -1] <- apply(train[, -1], 2, as.numeric)
  test[, -1] <- apply(test[, -1], 2, as.numeric)
  
  rf_classifier <- randomForest(diagnosis ~ ., data = train, importance = TRUE, ntree = 200)
  predictions <- predict(rf_classifier, test)
  prob_predictions <- predict(rf_classifier, test, type = "prob")
  
  return(list(model = rf_classifier, predictions = predictions, prob_predictions = prob_predictions))
}

run_analysis <- function(metadata_train, metadata_test, df, features, output_prefix) {
  train <- df[metadata_train$RUN, ]
  test <- df[metadata_test$RUN, ]
  
  rf_results <- run_random_forest(train, test, features)
  
  results <- data.frame(
    prediction = rf_results$predictions,
    truth = test$diagnosis,
    prob = rf_results$prob_predictions,
    stringsAsFactors = FALSE
  )

  #saveRDS(results, paste0(output_prefix, "_probabilities.rds"))
}

training_accuracy <- function(metadata, df, features, output_prefix) {  
  results <- data.frame(sample=character(), pred=character(), prob=numeric(), stringsAsFactors=FALSE)
  gini <- data.frame(feature=features)
  
  for (i in 1:nrow(metadata)) {
    print(i)
    
    train <- df[metadata[-i, "RUN"], c("diagnosis", features)]
    test <- df[metadata[i, "RUN"], c("diagnosis", features)]
    
    train[, -1] <- apply(train[, -1], 2, as.numeric)
    test[, -1] <- apply(test[, -1], 2, as.numeric)
    
    set.seed(10)
    rf_classifier <- randomForest(diagnosis ~ ., data=train, importance=TRUE, ntree=200)
    
    gini_df <- data.frame(rf_classifier$importance[,"MeanDecreaseGini",drop=F])
    colnames(gini_df) <- paste0("model",i)
    gini <- merge(gini, gini_df, by.x="feature", by.y="row.names")
    
    pred <- predict(rf_classifier, test)
    pred_prob <- predict(rf_classifier, test, type="prob")
    
    results <- rbind(results, data.frame(RUN=rownames(test), pred=pred, prob=pred_prob))
  }
  
  results <- merge(results, metadata[, c("RUN", "diagnosis")], by="RUN")
  results <- results %>%
    dplyr::rename(truth = diagnosis, prediction = pred)
  
  #saveRDS(gini, paste0(output_prefix, "_gini.rds"))
  #saveRDS(results, paste0(output_prefix, "_probabilities.rds"))
  
}

# Run RF for ASEs
training_results_ASEs <- training_accuracy(metadata_training, df_ASEs, features_ASE, "training_ASEs_")
test_metrics_ASEs <- run_analysis(metadata_training, metadata_test, df_ASEs, features_ASE, "test_ASEs_")
validation_metrics_ASEs <- run_analysis(metadata_training, metadata_validation, df_ASEs, features_ASE, "validation_ASEs_")

# Run RF for Cell Type Prop
training_results_Cell <- training_accuracy(metadata_training, df_cell, features_cell, "training_Cell_")
test_metrics_Cell <- run_analysis(metadata_training, metadata_test, df_cell, features_cell, "test_Cell_")
validation_metrics_Cell <- run_analysis(metadata_training, metadata_validation, df_cell, features_cell, "validation_Cell_")

# Run RF for Combined
training_results_Combined <- training_accuracy(metadata_training, df_all, features_all, "training_Combined_")
test_metrics_Combined <- run_analysis(metadata_training, metadata_test, df_all, features_all, "test_Combined_")
validation_metrics_Combined <- run_analysis(metadata_training, metadata_validation, df_all, features_all, "validation_Combined_")
