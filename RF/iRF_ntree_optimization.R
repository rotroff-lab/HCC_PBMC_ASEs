library(randomForest)
library(readxl)
library(dplyr)
library(matrixStats)

# Import PSI table. Scale and Center data
PSI <- readRDS("PSI.rds")
PSI[, grep("^SRR", names(PSI))] <- as.data.frame(scale(PSI[, grep("^SRR", names(PSI))], center = TRUE))

# Subset for Immune ASEs
all_immune_ASEs_validation <- readRDS("all_immune_ASEs_validation.rds")
psi <- PSI %>% filter(AS_ID_short %in% all_immune_ASEs_validation) %>% dplyr::select(-ID,-gene_name, -type, -annotated)
psi <- t(psi)
psi <- psi[-1,]

# import sample partitions and subset for training set
metadata_sample_partitions <- readRDS("partition.rds")
metadata_training <- metadata_sample_partitions %>% filter(model == "train")


# Add diagnosis data to psi table
diagnosis <- metadata_sample_partitions[,c(1,4)]
colnames(diagnosis) <- c("RUN","diagnosis")
df <- merge(diagnosis, psi, by.x="RUN", by.y="row.names")
df$diagnosis <- as.factor(df$diagnosis)
rownames(df) <- df$RUN
df$RUN <- NULL

# Extract ASEs
features <- names(df)
features <- features[-1]


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
  OOB_max = as.numeric(),
  ntreeRF = as.numeric()
)

# Run RF
ntrees <- c(100,200,300,400)
row_num <- 1 

for (ntree in ntrees) {
  results <- data.frame(sample=as.character(),
                        pred=as.character())
  gini <- data.frame(feature=features)
  oob_error_rates <- data.frame(ntree = seq(ntree))

  
  #LOOCV
  
  for(i in 1:nrow(metadata_training)){
    print(i)
    train <- metadata_training[-i,"RUN"]
    test <- metadata_training[i,"RUN"]
    train <- df[train,]
    test <- df[test,]
    x <- train
    x$diagnosis <- NULL
    x <- as.matrix(x)
    x <- apply(x, c(1, 2), as.numeric)
    y <- train$diagnosis
    
    #random forest on training
    train <- train[,c("diagnosis",features)]
    
    if (length(features) > 1) {
      
      train[, -1] <- apply(train[, -1], 2, as.numeric)
      test[, -1] <- apply(test[, -1], 2, as.numeric)
      set.seed(10)
      
    } else {
      train[,2] <- as.numeric(train[,2])
      test[,2] <- as.numeric(test[,2])
      
    }
    rf_classifier = randomForest(diagnosis ~ ., data=train, importance=TRUE, ntree=ntree)
    
    #gini
    gini_df <- data.frame(rf_classifier$importance[,"MeanDecreaseGini",drop=F])
    colnames(gini_df) <- paste0("model",i)
    gini <- merge(gini, gini_df, by.x="feature", by.y="row.names")
    
    #oob
    oob_error_rates_df <- data.frame(rf_classifier$err.rate[, "OOB"])
    colnames(oob_error_rates_df) <- paste0("model",i)
    oob_error_rates <- merge(oob_error_rates, oob_error_rates_df, by.x="ntree", by.y="row.names")
    
    #random forest on test
    rfpred <- predict(rf_classifier, test)
    rfpred_prob <- predict(rf_classifier, test, type = "prob")
    results2 <- data.frame(sample=rownames(test),pred=rfpred, prob=rfpred_prob)
    results <- rbind(results, results2)
  }
  
  results <- merge(results, train[,"diagnosis",drop=F],by.x="sample", by.y="row.names")
  results_sum <- as.data.frame(do.call(cbind, lapply(c("group1","group2","group3"), truthtable, pred.truth=results)))
  
  rownames(gini) <- gini$feature
  gini$feature <- NULL
  gini <- as.matrix(gini)
  gini_plot <- data.frame(feature=rownames(gini), min=rowMins(gini), max= rowMaxs(gini), mean = rowMeans(gini))
  
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
    ntreeRF = ntree,
    stringsAsFactors = F
  )
  
  iRF_list[row_num, ] <- myrow
  
  row_num <- row_num + 1 
  
}

saveRDS(iRF_list, "ntree_optimization.rds")




