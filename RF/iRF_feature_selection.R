library(randomForest)
library(dplyr)
library(matrixStats)
library(foreach)
library(doParallel)
library(pROC)

# Import PSI table. Scale and Center data
PSI <- readRDS("PSI.rds")
PSI[, grep("^SRR", names(PSI))] <- as.data.frame(scale(PSI[, grep("^SRR", names(PSI))], center = TRUE))

# Subset for Immune ASEs
all_immune_ASEs_validation <- readRDS("all_immune_ASEs_validation.rds")
psi <- PSI %>% filter(AS_ID_short %in% all_immune_ASEs_validation) %>% dplyr::select(-ID,-gene_name, -type, -annotated)
psi <- t(psi)
psi <- psi[-1,]

# Import sample partitions
metadata_sample_partitions <- metadata_sample_partitions %>% select(RUN, two_group, model)
metadata_sample_partitions <- metadata_sample_partitions %>%
  dplyr::rename(diagnosis = two_group)
metadata_training <- metadata_sample_partitions %>% filter(model == "train")
metadata_test <- metadata_sample_partitions %>% filter(model == "test")
metadata_validation <- metadata_sample_partitions %>% filter(model == "validation")


# Add diagnosis data to psi table
diagnosis <- metadata_sample_partitions[,c("RUN","diagnosis")]
df_ASEs <- merge(diagnosis, psi, by.x="RUN", by.y="row.names")
df_ASEs$diagnosis <- as.factor(df_ASEs$diagnosis)
rownames(df_ASEs) <- df_ASEs$RUN
df_ASEs$RUN <- NULL

# Subset Trainging Samples
train <- metadata_training$RUN
df2 <- df[train,]
df2[, -1] <- apply(df2[, -1], 2, as.numeric)
diagnosis <- df2$diagnosis


# Extract Features
features <- names(df2)[-1]

truthtable <- function(group, pred.truth){
  TP <- nrow(pred.truth[which(pred.truth$pred %in% group & pred.truth$diagnosis %in% group),])
  TN <- nrow(pred.truth[which(!pred.truth$pred %in% group & !pred.truth$diagnosis %in% group),])
  
  FP <- nrow(pred.truth[which(pred.truth$pred %in% group & !pred.truth$diagnosis %in% group),])
  FN <- nrow(pred.truth[which(!pred.truth$pred %in% group & pred.truth$diagnosis %in% group),])
  total <- nrow(pred.truth)
  
  myrow <- data.frame(
    accuracy= (TP+TN)/total,
    balanced_accuracy=((TP/(TP+FN))+(TN/(FP+TN)))/2,
    misclassification = (FP+FN)/total,
    sensitivity= TP/(TP+FN),
    specificity =TN/(TN+FP),
    PPV = TP/(TP+FP),
    NPV = TN/(TN+FP),
    stringsAsFactors = F
  )
  rownames(myrow) <- group
  myrow <- data.frame(t(myrow),stringsAsFactors = F)
  myrow <- apply(myrow, 2, function(x) round(x*100,digits = 2))
  return(myrow)
}


# Initialize df
iRF_list <- data.frame(
  num_features=as.numeric(),
  features=as.character(),
  balanced_accuracy_HCC=as.numeric(),
  balanced_accuracy_Healthy=as.numeric(),
  balanced_accuracy_HBV=as.numeric(),
  misclassification_HCC=as.numeric(),
  misclassification_Healthy=as.numeric(),
  misclassification_HBV=as.numeric(),
  sensitivity_HCC=as.numeric(),
  sensitivity_Healthy=as.numeric(),
  sensitivity_HBV=as.numeric(),
  specificity_HCC=as.numeric(),
  specificity_Healthy=as.numeric(),
  specificity_HBV=as.numeric(),
  PPV_HCC=as.numeric(),
  PPV_Healthy=as.numeric(),
  PPV_HBV=as.numeric(),
  NPV_HCC=as.numeric(),
  NPV_Healthy=as.numeric(),
  NPV_HBV=as.numeric(),
  OOB_min = as.numeric(),
  OOB_mean =as.numeric(),
  OOB_max = as.numeric()
)


# Select ntree
ntrees <- 200

# Run iterative RF LOOCV for feature selection
while (length(features) > 0) {
  results <- data.frame(sample = as.character(), pred = as.character())
  gini <- data.frame(feature = features)
  oob_error_rates <- data.frame(ntree = seq(ntrees))
  
  # Parallel LOOCV
  parallel_results <- foreach(i = 1:nrow(metadata_training), .packages = c('randomForest', 'foreach', 'doParallel')) %dopar% {
    train <- df2[metadata_training[-i, "RUN"], c("diagnosis", features)]  
    test <- df2[metadata_training[i, "RUN"], c("diagnosis", features)]  
    
    # Random forest on training
    train[, -1] <- lapply(train[, -1], as.numeric)
    test[, -1] <- lapply(test[, -1], as.numeric)
    rf_classifier = randomForest(diagnosis ~ ., data=train, importance=TRUE, ntree=ntrees)
    
    # Gini
    gini_df <- data.frame(rf_classifier$importance[, "MeanDecreaseGini", drop = F])
    colnames(gini_df) <- paste0("model", i)
    
    # OOB error rates
    oob_error_rates_df <- data.frame(rf_classifier$err.rate[, "OOB"])
    colnames(oob_error_rates_df) <- paste0("model", i)
    
    # Random forest on LOO sample
    rfpred <- predict(rf_classifier, test)
    rfpred_prob <- predict(rf_classifier, test, type = "prob")
    results2 <- data.frame(sample = rownames(test), pred = rfpred, prob = rfpred_prob)
    
    list(gini_df = gini_df, oob_error_rates_df = oob_error_rates_df, results2 = results2)
  }
  
  # Combine parallel results
  for (res in parallel_results) {
    gini <- merge(gini, res$gini_df, by.x = "feature", by.y = "row.names")
    oob_error_rates <- cbind(oob_error_rates, res$oob_error_rates_df)
    results <- rbind(results, res$results2)
  }
  
  metadata_training2 <- metadata_training[, c("RUN", "model_6")] %>% dplyr::rename(diagnosis = model_6)
  results <- merge(results, metadata_training2, by.x = "sample", by.y = "RUN")
  results_sum <- as.data.frame(do.call(cbind, lapply(c("group1","group2","group3"), truthtable, pred.truth=results)))
  
  ##remove lowest median gini feature. If > 500 remove 100 features. below remove 1 feature
  rownames(gini) <- gini$feature
  gini$feature <- NULL
  gini <- as.matrix(gini)
  gini_plot <- data.frame(feature=rownames(gini), min=rowMins(gini), max= rowMaxs(gini), mean = rowMeans(gini))
  
  features <- gini_plot %>%
    arrange(mean) %>%
    slice_tail(n = max(length(features) - ifelse(length(features) > 2000, 500, ifelse(length(features) > 500, 100, 1)), 1)) %>%
    pull(feature)
  
  ##OOB errors
  oob_error_rates$ntree <- NULL
  oob_error_rates <- as.matrix(oob_error_rates)
  
  oob_error_min <- min(oob_error_rates)
  oob_error_mean <- mean(unlist(oob_error_rates))
  oob_error_max <- max(oob_error_rates)
  
  
  myrow <- data.frame(
    num_features=length(features),
    features=paste0(features, collapse=","),
    balanced_accuracy_HCC= results_sum$group1[2],
    balanced_accuracy_Healthy=results_sum$group2[2],
    balanced_accuracy_HBV=results_sum$group3[2],
    misclassification_HCC =results_sum$group1[3],
    misclassification_Healthy =results_sum$group2[3],
    misclassification_HBV =results_sum$group3[3],
    sensitivity_HCC=results_sum$group1[4],
    sensitivity_Healthy=results_sum$group2[4],
    sensitivity_HBV=results_sum$group3[4],
    specificity_HCC =results_sum$group1[5],
    specificity_Healthy =results_sum$group2[5],
    specificity_HBV =results_sum$group3[5],
    PPV_HCC =results_sum$group1[6],
    PPV_Healthy =results_sum$group2[6],
    PPV_HBV =results_sum$group3[6],
    NPV_HCC =results_sum$group1[7],
    NPV_Healthy =results_sum$group2[7],
    NPV_HBV =results_sum$group3[7],
    OOB_min = oob_error_min,
    OOB_mean = oob_error_mean,
    OOB_max = oob_error_max,
    stringsAsFactors = F
  )
  
  iRF_list[length(features), ] <- myrow
}

# Stop the parallel backend
stopCluster(cl)


saveRDS(iRF_list, "feature_selection.rds")




