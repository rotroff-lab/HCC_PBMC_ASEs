library(stringr)

comparison <- as.character(commandArgs(TRUE[1]))
outputname <- "output"


numeric_df <- function(df){
  saverownames <- rownames(df)
  savecolnames <- colnames(df)
  df <- data.frame(sapply(df, function(x) as.numeric(x)))
  rownames(df) <- saverownames
  colnames(df) <- savecolnames
  return(df)
}

rmats_files <- list.files(paste0("rmats/",comparison), pattern = "SRP")
dir.create(paste0("rmats_results"))
dir.create(paste0("rmats_results/",comparison))
dir.create(paste0("rmats_results/",comparison,"/tables"))

rmats <- function(SRP_file, AStype, outputname){
  print(SRP_file)
  path <- paste0("rmats/",comparison,"/",SRP_file)
  df <- read.csv(paste0(path,"/",outputname,"/",AStype,".MATS.JC.txt"), header=T, sep="\t")
  print("rmats file read")
  
  samplenames_1 <- read.table(paste0(path,"/input/condition.txt"), sep=",")
  samplenames_1 <- str_extract(samplenames_1, "SRR\\d*")
  samplenames_2 <- read.table(paste0(path,"/input/control.txt"), sep=",")
  samplenames_2 <- str_extract(samplenames_2, "SRR\\d*")
  
  print("sample data read")
  
  SJC_1 <- str_split_fixed(df$SJC_SAMPLE_1, ",", length(samplenames_1))
  colnames(SJC_1) <- samplenames_1
  
  SJC_2 <- str_split_fixed(df$SJC_SAMPLE_2, ",", length(samplenames_2))
  colnames(SJC_2) <- samplenames_2
  
  SJC <- data.frame(cbind(SJC_1, SJC_2))
  
  if(AStype=="SE"){
    rownames(SJC) <- paste(AStype, df$chr, df$strand, df$exonStart_0base, df$exonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, sep=",")
  }
  if(AStype %in% c("A3SS","A5SS")){
    rownames(SJC) <- paste(AStype, df$chr, df$strand, df$longExonStart_0base, df$longExonEnd, df$shortES, df$shortEE, df$flankingES, df$flankingEE, sep=",")
  }
  if(AStype %in% c("RI")){
    rownames(SJC) <- paste(AStype, df$chr, df$strand, df$riExonStart_0base, df$riExonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, sep=",")
  }
  
  SJC <- numeric_df(SJC)
  SJC$ID <-rownames(SJC)
  
  #IJC
  IJC_1 <- str_split_fixed(df$IJC_SAMPLE_1, ",", length(samplenames_1))
  colnames(IJC_1) <- samplenames_1
  
  IJC_2 <- str_split_fixed(df$IJC_SAMPLE_2, ",", length(samplenames_2))
  colnames(IJC_2) <- samplenames_2
  
  IJC <- data.frame(cbind(IJC_1, IJC_2))
  
  if(AStype=="SE"){
    rownames(IJC) <- paste(AStype, df$chr, df$strand, df$exonStart_0base, df$exonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, sep=",")
  }
  if(AStype %in% c("A3SS","A5SS")){
    rownames(IJC) <- paste(AStype, df$chr, df$strand, df$longExonStart_0base, df$longExonEnd, df$shortES, df$shortEE, df$flankingES, df$flankingEE, sep=",")
  }
  if(AStype %in% c("RI")){
    rownames(IJC) <- paste(AStype, df$chr, df$strand, df$riExonStart_0base, df$riExonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, sep=",")
  }
  
  IJC <- numeric_df(IJC)
  IJC$ID <-rownames(IJC)
  
  mylist <- list(SJC, IJC)
  return(mylist)
}


#SE____________________________________________________________________________
removemissing <- function(SJC, IJC){
  #for each ASE event, get read counts
  SJC$ID <- NULL
  IJC$ID <- NULL 
  loopASE <- function(ASE_ID, SJC, IJC){
    counts <- as.numeric(SJC[ASE_ID,]) + as.numeric(IJC[ASE_ID,])
    total <- length(counts)
    sumzero <- length(counts[which(counts==0)])
    prop <- sumzero/total
    ##print(prop)
    if(prop<=0.1){
      return(ASE_ID)
    }
  }
  lapply(rownames(SJC), loopASE, SJC=SJC, IJC=IJC)
}


removeme <- function(SRP, keep){
  ##print(nrow(SRP))
  SRP <- SRP[which(SRP$ID %in% keep),]
  ##print(nrow(SRP))
  return(SRP)
}


JC_list <- lapply(rmats_files, rmats, AStype="SE", outputname=outputname)
SJC_list <- lapply(JC_list, function(x) x[[1]])
IJC_list <- lapply(JC_list, function(x) x[[2]])

keep <- unlist(mapply(removemissing, SJC_list, IJC_list))

IJC_list <- lapply(IJC_list, removeme, keep=keep)
SJC_list <- lapply(SJC_list, removeme, keep=keep)

SE_SJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), SJC_list)
SE_IJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), IJC_list)

#A3____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="A3SS", outputname=outputname)
SJC_list <- lapply(JC_list, function(x) x[[1]])
IJC_list <- lapply(JC_list, function(x) x[[2]])

keep <- unlist(mapply(removemissing, SJC_list, IJC_list))
  
IJC_list <- lapply(IJC_list, removeme, keep=keep)
SJC_list <- lapply(SJC_list, removeme, keep=keep)

A3SS_SJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), SJC_list)
A3SS_IJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), IJC_list)

#A5____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="A5SS", outputname=outputname)
SJC_list <- lapply(JC_list, function(x) x[[1]])
IJC_list <- lapply(JC_list, function(x) x[[2]])

keep <- unlist(mapply(removemissing, SJC_list, IJC_list))

IJC_list <- lapply(IJC_list, removeme, keep=keep)
SJC_list <- lapply(SJC_list, removeme, keep=keep)

A5SS_SJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), SJC_list)
A5SS_IJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), IJC_list)

#RI____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="RI", outputname=outputname)
SJC_list <- lapply(JC_list, function(x) x[[1]])
IJC_list <- lapply(JC_list, function(x) x[[2]])

keep <- unlist(mapply(removemissing, SJC_list, IJC_list))

IJC_list <- lapply(IJC_list, removeme, keep=keep)
SJC_list <- lapply(SJC_list, removeme, keep=keep)

RI_SJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), SJC_list)
RI_IJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), IJC_list)

#final tables__________________________________________________________________

SJC <- do.call("rbind", list(SE_SJC, RI_SJC, A3SS_SJC, A5SS_SJC))
saveRDS(SJC, paste0("rmats_results/",comparison,"/tables/SJC_JC.rds"))

IJC <- do.call("rbind", list(SE_IJC, RI_IJC, A3SS_IJC, A5SS_IJC))
saveRDS(IJC, paste0("rmats_results/",comparison,"/tables/IJC_JC.rds"))

ID_table <- data.frame(ID=SJC$ID)
saveRDS(ID_table, paste0("rmats_results/",comparison,"/tables/splice_ID_JC.rds"))

