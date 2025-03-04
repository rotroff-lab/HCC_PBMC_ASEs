library(parallel)
library(stringr)

#file to impute missing read values if less than 10% of junctions have no coverage

comparison <- as.character(commandArgs(TRUE)[2])
SRP <- as.character(commandArgs(TRUE)[1])
path <- as.character(commandArgs(TRUE)[3])
output <- as.character(commandArgs(TRUE)[4])

loop_files <- function(file){
  df <- read.table(paste0(path,"/rmats/",comparison,"/", SRP,"/",output,"/",file), header=T)
  myrownames <- paste0("ID",df$ID)
  rownames(df) <- myrownames
  
  numeric_df <- function(df){
    df <- data.frame(df)
    saverownames <- rownames(df)
    savecolnames <- colnames(df)
    df <- data.frame(sapply(df, function(x) as.numeric(x)))
    rownames(df) <- saverownames
    colnames(df) <- savecolnames
    return(df)
  }

  #for each individual AS event, determine the number of NAs
  
  #proportion of NAs
  SJC_1 <- str_split_fixed(df$SJC_SAMPLE_1, ",", unique(str_count(df$SJC_SAMPLE_1, ',')+1))
  rownames(SJC_1) <- myrownames
  SJC_2 <- str_split_fixed(df$SJC_SAMPLE_2, ",", unique(str_count(df$SJC_SAMPLE_2, ',')+1))
  rownames(SJC_2) <- myrownames
  
  IJC_1 <- str_split_fixed(df$IJC_SAMPLE_1, ",", unique(str_count(df$IJC_SAMPLE_1, ',')+1))
  rownames(IJC_1) <- myrownames
  IJC_2 <- str_split_fixed(df$IJC_SAMPLE_2, ",", unique(str_count(df$IJC_SAMPLE_2, ',')+1))
  rownames(IJC_2) <- myrownames
  
  SJC <- as.matrix(numeric_df(cbind(SJC_1, SJC_2)))
  rownames(SJC) <- myrownames
  IJC <- as.matrix(numeric_df(cbind(IJC_1, IJC_2)))
  rownames(IJC) <- myrownames
  
  total <- data.frame(SJC + IJC)
  total$num_0 <- apply(total, 1, function(x) sum(x == 0))
  total$percent <- total$num_0/ncol(total)
  total$impute <- ifelse(total$percent<0.1&total$percent>0,"yes","no")
  
  if(any(total$impute=="yes")){
    
    #which to impute
    total_1 <- data.frame(as.matrix(numeric_df(SJC_1))+as.matrix(numeric_df(IJC_1)))
    rownames(total_1) <- myrownames
    total_1$impute <- total$impute
    total_1 <- total_1[which(total_1$impute=="yes"),]
  
    total_2 <- data.frame(as.matrix(numeric_df(SJC_2))+as.matrix(numeric_df(IJC_2)))
    rownames(total_2) <- myrownames
    total_2$impute <- total$impute
    total_2 <- total_2[which(total_2$impute=="yes"),]
  
    
    #get medians
    SJC_medians <- apply(SJC, 1, function(x) round(median(x)))
    names(SJC_medians) <- myrownames
    IJC_medians <- apply(IJC, 1, function(x) round(median(x)))
    names(IJC_medians) <- myrownames
    
    impute_loop <- function(i, df, impute, median){
      impute_pos <- which(impute[i,]==0)
  
      if(any(!is.na(impute_pos))){
        df[i,impute_pos] <- median[[i]]
      }
        return(df[i,]) 
    }
    
    IJC_1_imputed <- do.call("rbind",lapply(rownames(total_1), impute_loop, impute=total_1, df=IJC_1, median=IJC_medians)) #works
    SJC_1_imputed <- do.call("rbind",lapply(rownames(total_1), impute_loop, impute=total_1, df=SJC_1, median=SJC_medians))
    SJC_2_imputed <- do.call("rbind",lapply(rownames(total_2), impute_loop, impute=total_2, df=SJC_2, median=SJC_medians))
    IJC_2_imputed <- do.call("rbind",lapply(rownames(total_2), impute_loop, impute=total_2, df=IJC_2, median=IJC_medians))
    
    df2 <- df[rownames(total_1),]
    df2$IJC_SAMPLE_1 <- apply(IJC_1_imputed,1, function(x) paste(x,collapse=","))
    df2$SJC_SAMPLE_1 <- apply(SJC_1_imputed,1, function(x) paste(x,collapse=","))
    df2$IJC_SAMPLE_2 <- apply(IJC_2_imputed,1, function(x) paste(x,collapse=","))
    df2$SJC_SAMPLE_2 <- apply(SJC_2_imputed,1, function(x) paste(x,collapse=","))
    
    write.table(df2, paste0(path,"/rmats/",comparison,"/", SRP,"/", output, "/imputed_",file), sep="\t", col.names = T, row.names = F, quote=F)
  }else{
  write.table(df, paste0(path,"/rmats/",comparison,"/", SRP,"/", output, "/imputed_",file), sep="\t", col.names = T, row.names = F, quote=F)
  }
}  


rmats_files <- list.files(paste0("rmats/",comparison,"/",SRP,"/", output, "/"), pattern = "JCEC.raw")
mclapply(rmats_files, loop_files)

