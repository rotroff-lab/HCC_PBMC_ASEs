library(stringr)
SRP <- as.character(commandArgs(TRUE)[1])
comparison <- as.character(commandArgs(TRUE)[2])
path <- as.character(commandArgs(TRUE)[3])
output <- as.character(commandArgs(TRUE)[4])

numeric_df <- function(df){
  df <- data.frame(df)
  saverownames <- rownames(df)
  savecolnames <- colnames(df)
  df <- data.frame(sapply(df, function(x) as.numeric(x)))
  rownames(df) <- saverownames
  colnames(df) <- savecolnames
  return(df)
}


loop_files <- function(file1, file2){
  #read in table with inclusion levels
  df_post <- read.table(paste0(path,"/rmats/",comparison,"/", SRP,"/", output, "/",file1), header=T)
  df_pval <- read.table(paste0(path,"/rmats/",comparison,"/", SRP,"/", output, "/",file2), header=T)
  
  df <- df_post[which(df_post$ID %in% df_pval$ID),]
  df_post <- df_post[which(!df_post$ID %in% df_pval$ID),]
  
  df <- df[order(df_pval$ID),]
  
  df$IJC_SAMPLE_1 <- df_pval$IJC_SAMPLE_1
  df$SJC_SAMPLE_1 <- df_pval$SJC_SAMPLE_1
  df$IJC_SAMPLE_2 <- df_pval$IJC_SAMPLE_2
  df$SJC_SAMPLE_2 <- df_pval$SJC_SAMPLE_2
  df$PValue <- df_pval$PValue
  
  #Calculate IncLevel1 and IncLevel2
  #IncLevel1
  IJC_1 <- str_split_fixed(df$IJC_SAMPLE_1, ",", unique(str_count(df$IJC_SAMPLE_1, ',')+1))
  SJC_1 <- str_split_fixed(df$SJC_SAMPLE_1, ",", unique(str_count(df$SJC_SAMPLE_1, ',')+1))
  I_1 <- as.matrix(numeric_df(IJC_1))/df$IncFormLen
  S_1 <- as.matrix(numeric_df(SJC_1))/df$SkipFormLen
  IncLevel1 <- round(I_1/(I_1+S_1),digits=3)
  IncLevel1_vec <- apply(IncLevel1,1,function(x) paste0(x, collapse = ","))
  IncLevel1_mean <- apply(IncLevel1,1,function(x) mean(x, na.rm=T))
  df$IncLevel1 <- IncLevel1_vec
  
  IJC_2 <- str_split_fixed(df$IJC_SAMPLE_2, ",", unique(str_count(df$IJC_SAMPLE_2, ',')+1))
  SJC_2 <- str_split_fixed(df$SJC_SAMPLE_2, ",", unique(str_count(df$SJC_SAMPLE_2, ',')+1))
  I_2 <- as.matrix(numeric_df(IJC_2))/df$IncFormLen
  S_2 <- as.matrix(numeric_df(SJC_2))/df$SkipFormLen
  IncLevel2 <- round(I_2/(I_2+S_2),digits=3)
  IncLevel2_vec <- apply(IncLevel2,1,function(x) paste0(x, collapse = ","))
  IncLevel2_mean <- apply(IncLevel2,1,function(x) mean(x, na.rm=T))
  df$IncLevel2 <- IncLevel2_vec

  #Calculate IncLevelDifference
  #mean of Inclevel1 - mean of Inclevel2
  IncLevelDifference_vec <- IncLevel1_mean - IncLevel2_mean
  df$IncLevelDifference <- round(IncLevelDifference_vec, digits=3)
  
  df_imputed <- rbind(df, df_post)
  
  df_imputed$FDR <- NA
  
  #if IncLevel1 or IncLevel2 contain any NAs change pvalue to NA
  removeP <- function(x){
    if(grepl("NA|NaN",x[["IncLevel1"]])){
      x[["PValue"]] <- NA
    }
    if(grepl("NA|NaN",x[["IncLevel2"]])){
      x[["PValue"]] <- NA
    }
    return(x)
  }
  df_post <- data.frame(t(apply(df_imputed, 1, removeP)))
  
  write.table(df_imputed, paste0(path,"/rmats/",comparison,"/", SRP,"/", output, "/adjusted_",file1), sep="\t", col.names = T, row.names = F, quote = F)
}  

rmats_files1 <- list.files(paste0("rmats/",comparison,"/",SRP,"/", output, "/"), pattern = "MATS.JCEC")
rmats_files2 <- list.files(paste0("rmats/",comparison,"/",SRP,"/", output, "/"), pattern = "Result_P")
mapply(loop_files, rmats_files1, rmats_files2)

