
setwd("/proj/RotroffDLab/noah/HCC_splicing/2024_HCC_PBMCs/")


library(poolr)
library(matrixStats)
library(tidyr)
library(dplyr)
library(metafor)
library(foreach)
library(doParallel)

cl <- makePSOCKcluster(20)
registerDoParallel(cl)

comparison <- "hcc_v_nohcc"
splice_IDs <- readRDS("splice_ID_filtered.rds")

# Import meta data
metadata <- readRDS("metadata.rds")

##join validation, test and training
psi <- readRDS("PSI.rds")
psi <- psi[complete.cases(psi),]

rownames(psi) <- psi$ID
metadata <- metadata[which(metadata$RUN %in% colnames(psi)),]
PBMC_studies <- unique(metadata$SRA.STUDY)

findmissing <- function(study){
  samples <- metadata[which(metadata$SRA.STUDY==study),"RUN"]
  psi <- psi[,which(colnames(psi) %in% c("ID",samples))]
  na_count_per_row <- apply(psi, 1, function(x) sum(is.na(x)))/(ncol(psi)-1)
  df <- data.frame(ID=psi$ID, sum_missing=na_count_per_row)
  return(df[which(df$sum_missing<0.1),"ID"])
}

keep <- lapply(unique(metadata$SRA.STUDY),findmissing)
features <- unlist(keep)
feature_counts <- data.frame(table(features))

# ASEs measured in all six studies
best <- feature_counts[which(feature_counts$Freq>=6),"features"]
psi <- psi[which(psi$ID %in% best),]

#split by chr for faster processing
events <- data.frame(ID=psi$ID)
events$chr <- sub("RI,|SE,|A3SS,|A5SS,","",events$ID)
events$chr <- sub(",.*","",events$chr)

packages=c("metafor")


# Run Meta Analysis
foreach(chr = unique(events$chr), .packages = packages) %dopar% {
  events_sub <- events[which(events$chr==chr),]

  df1 <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(df1) <- c("yi","vi","sei","zi","pval","ci.lb","ci.up","w","incleveldiff","study","AS_ID")
  saveRDS(df1,paste0("meta_analysis/",comparison,"_forest_numbers_",chr,".rds"))
  
  df2 <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(df2) <- c("AS_ID", "pvalue", "meanincleveldiff", "sample_num","ci.lb","ci.up")
  saveRDS(df2, paste0("meta_analysis/",comparison,"_meta_numbers_",chr,".rds"))
  
  for(k in 1:length(events_sub$ID)){
    AS_ID <- events_sub$ID[[k]]
    variance_list <- data.frame()
    
    for(i in 1:length(PBMC_studies)){
      study <- PBMC_studies[[i]]
      met <- metadata[which(metadata$SRA.STUDY %in% study),c("RUN", comparison)]
      colnames(met) <- c("RUN","group")
      sample1 <- met[which(met$group=="group1"),"RUN"]
      sample2 <- met[which(met$group=="group2"),"RUN"]
      n1i <- length(sample1)
      n2i <- length(sample2)
      as1 <- as.numeric(psi[AS_ID,sample1])
      as2 <- as.numeric(psi[AS_ID,sample2])
      prop_nas <- length(as1[which(is.na(as1))]) + length(as1[which(is.na(as2))])/(n1i + n2i)
      
      if(prop_nas<0.1){
        m1i <- mean(as1, na.rm=T)
        sd1i <- sd(as1, na.rm=T)
        m2i <- mean(as2, na.rm=T)
        sd2i <- sd(as2, na.rm=T)
        variance <- summary(escalc(measure="MD", m1i=m1i, m2i=m2i, sd1i=sd1i, sd2i=sd2i, n1i=n1i, n2i=n2i))
        variance$w <- n1i+n2i
        variance$incleveldiff <- m1i-m2i
        variance$study <- study
        variance_list <- rbind(variance_list, variance)
      }
    }
    variance_list$AS_ID <- AS_ID
    
    forest_all <- readRDS(paste0("meta_analysis/",comparison,"_forest_numbers_",chr,".rds"))
    forest_all <- rbind(forest_all, variance_list)
    saveRDS(forest_all,paste0("meta_analysis/",comparison,"_forest_numbers_",chr,".rds"))
    
    meanincleveldiff <- mean(variance_list$incleveldiff)
    sample_num <- nrow(variance_list)
    
    meta_model <- tryCatch({rma.mv(yi=variance_list$yi, V=variance_list$vi, W=variance_list$w)}, error=function(msg){return(NA)})

    if(all(is.na(meta_model))){
      meta_results <- data.frame(AS_ID=AS_ID, pvalue=NA, meanincleveldiff=meanincleveldiff, sample_num=sample_num, ci.lb=NA, ci.ub=NA)
    }else{
      res1 <- summary(meta_model)
      meta_results <- data.frame(AS_ID=AS_ID, pvalue=res1$pval, meanincleveldiff=meanincleveldiff, 
                                 sample_num=sample_num, ci.lb=res1$ci.lb, ci.ub=res1$ci.ub)
    }
    
    meta_all <- readRDS(paste0("meta_analysis/",comparison,"_meta_numbers_",chr,".rds"))
    meta_all <- rbind(meta_all, meta_results)
    saveRDS(meta_all,paste0("meta_analysis/",comparison,"_meta_numbers_",chr,".rds"))
    
  }
  
}


##Combine Meta Analysis and adjust p value
meta_list <- list.files("meta_analysis", pattern=paste0(comparison,"_meta_numbers_"), full.names = T)
make_list <- function(meta){
  meta_results <- readRDS(meta)
  return(meta_results)
}
meta_all <- lapply(meta_list, make_list)
meta_results <- do.call("rbind",meta_all)

meta_results <- meta_results[!duplicated(meta_results),]
meta_results <- merge(meta_results, splice_IDs, by.x ="AS_ID")
meta_results <- meta_results %>% filter(annotated == "annotated")
meta_results$fdr <- p.adjust(meta_results$pvalue, method="fdr")
saveRDS(meta_results, "meta_analysis_ASE.rds")

