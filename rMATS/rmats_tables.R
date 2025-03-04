library(stringr)

comparison <- as.character(commandArgs(TRUE[1]))
outputname <- "output"


ID_table <- readRDS(paste0("rmats_results/",comparison,"/tables/splice_ID_filtered.rds"))
ID_keep <- ID_table$ID

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
  df <- read.csv(paste0(path,"/",outputname,"/adjusted_",AStype,".MATS.JCEC.txt"), header=T, sep="\t")

    if(nrow(df[is.na(df$chr),])>0){
    print("Weird NA Problem:")
    print(AStype)
  }
  df <- df[!is.na(df$chr),]
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
  
  pvalues <- data.frame(ID=IJC$ID, pvalues=df$PValue)
  colnames(pvalues) <- c("ID",SRP_file)
  incleveldiff <- data.frame(ID=IJC$ID, incleveldiff=df$IncLevelDifference)
  colnames(incleveldiff) <- c("ID",SRP_file)
  
  #filter for AS_events with junction counts in at least 10 studies
  SJC <- SJC[which(SJC$ID %in% ID_keep),]
  IJC <- IJC[which(IJC$ID %in% ID_keep),]
  pvalues <- pvalues[which(pvalues$ID %in% ID_keep),]
  incleveldiff <- incleveldiff[which(incleveldiff$ID %in% ID_keep),]
  
  mylist <- list(SJC, IJC)
  mylist[[3]] <- pvalues
  mylist[[4]] <- incleveldiff
  return(mylist)
}


#SE____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="SE", outputname=outputname)
SJC_list <- lapply(JC_list, function(x) x[[1]])
IJC_list <- lapply(JC_list, function(x) x[[2]])
pval_list <- lapply(JC_list, function(x) x[[3]])
incleveldiff_list <- lapply(JC_list, function(x) x[[4]])


SE_SJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), SJC_list)
SE_IJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), IJC_list)
SE_pval <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), pval_list)
SE_incleveldiff <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), incleveldiff_list)

#A3____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="A3SS", outputname=outputname)
SJC_list <- lapply(JC_list, function(x) x[[1]])
IJC_list <- lapply(JC_list, function(x) x[[2]])
pval_list <- lapply(JC_list, function(x) x[[3]])
incleveldiff_list <- lapply(JC_list, function(x) x[[4]])

A3SS_SJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), SJC_list)
A3SS_IJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), IJC_list)
A3SS_pval <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), pval_list)
A3SS_incleveldiff <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), incleveldiff_list)

#A5____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="A5SS", outputname=outputname)
SJC_list <- lapply(JC_list, function(x) x[[1]])
IJC_list <- lapply(JC_list, function(x) x[[2]])
pval_list <- lapply(JC_list, function(x) x[[3]])
incleveldiff_list <- lapply(JC_list, function(x) x[[4]])

A5SS_SJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), SJC_list)
A5SS_IJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), IJC_list)
A5SS_pval <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), pval_list)
A5SS_incleveldiff <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), incleveldiff_list)

#RI____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="RI", outputname=outputname)
SJC_list <- lapply(JC_list, function(x) x[[1]])
IJC_list <- lapply(JC_list, function(x) x[[2]])
pval_list <- lapply(JC_list, function(x) x[[3]])
incleveldiff_list <- lapply(JC_list, function(x) x[[4]])

RI_SJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), SJC_list)
RI_IJC <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), IJC_list)
RI_pval <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), pval_list)
RI_incleveldiff <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), incleveldiff_list)

#final tables__________________________________________________________________

SJC <- do.call("rbind", list(SE_SJC, RI_SJC, A3SS_SJC, A5SS_SJC))
saveRDS(SJC, paste0("rmats_results/",comparison,"/tables/SJC.rds"))

IJC <- do.call("rbind", list(SE_IJC, RI_IJC, A3SS_IJC, A5SS_IJC))
saveRDS(IJC, paste0("rmats_results/",comparison,"/tables/IJC.rds"))

ID_table <- data.frame(ID=SJC$ID)
saveRDS(ID_table, paste0("rmats_results/",comparison,"/tables/splice_ID.rds"))

pvalues <- do.call("rbind", list(SE_pval, RI_pval, A3SS_pval, A5SS_pval))
saveRDS(pvalues, paste0("rmats_results/",comparison,"/tables/pvals.rds"))

incleveldiff <- do.call("rbind", list(SE_incleveldiff, RI_incleveldiff, A3SS_incleveldiff, A5SS_incleveldiff))
saveRDS(incleveldiff, paste0("rmats_results/",comparison,"/tables/incleveldiff.rds"))



# Generate PSI tables
rmats <- function(SRP_file, AStype, outputname){
  print(SRP_file)
  path <- paste0("rmats/",comparison,"/",SRP_file)
  df <- read.csv(paste0(path,"/",outputname,"/adjusted_",AStype,".MATS.JCEC.txt"), header=T, sep="\t")
  if(nrow(df[is.na(df$chr),])>0){
    print("Weird NA Problem:")
    print(AStype)
  }
  df <- df[!is.na(df$chr),]
  print("rmats file read")
  
  samplenames_1 <- read.table(paste0(path,"/input/condition.txt"), sep=",")
  samplenames_1 <- str_extract(samplenames_1, "SRR\\d*")
  samplenames_2 <- read.table(paste0(path,"/input/control.txt"), sep=",")
  samplenames_2 <- str_extract(samplenames_2, "SRR\\d*")
  
  print("sample data read")
  
  IncLevel1 <- str_split_fixed(df$IncLevel1, ",", length(samplenames_1))
  colnames(IncLevel1) <- samplenames_1
  
  IncLevel2 <- str_split_fixed(df$IncLevel2, ",", length(samplenames_2))
  colnames(IncLevel2) <- samplenames_2
  
  IncLevel <- data.frame(cbind(IncLevel1, IncLevel2))
  
  if(AStype=="SE"){
    rownames(IncLevel) <- paste(AStype, df$chr, df$strand, df$exonStart_0base, df$exonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, sep=",")
  }
  if(AStype %in% c("A3SS","A5SS")){
    rownames(IncLevel) <- paste(AStype, df$chr, df$strand, df$longExonStart_0base, df$longExonEnd, df$shortES, df$shortEE, df$flankingES, df$flankingEE, sep=",")
  }
  if(AStype %in% c("RI")){
    rownames(IncLevel) <- paste(AStype, df$chr, df$strand, df$riExonStart_0base, df$riExonEnd, df$upstreamES, df$upstreamEE, df$downstreamES, df$downstreamEE, sep=",")
  }
  
  IncLevel <- numeric_df(IncLevel)
  IncLevel$ID <-rownames(IncLevel)
  
  
  
  #filter for AS_events with junction counts in at least 10 studies
  IncLevel <- IncLevel[which(IncLevel$ID %in% ID_keep),]
  
  return(IncLevel)
}


#SE____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="SE", outputname=outputname)
SE_IncLevel <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), JC_list)

#A3____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="A3SS", outputname=outputname)
A3SS_IncLevel <-  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), JC_list)

#A5____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="A5SS", outputname=outputname)
A5SS_IncLevel <-  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), JC_list)

#RI____________________________________________________________________________
JC_list <- lapply(rmats_files, rmats, AStype="RI", outputname=outputname)
RI_IncLevel <-  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ID", all.x = TRUE, all.y = TRUE), JC_list)

#final tables__________________________________________________________________

IncLevel <- do.call("rbind", list(SE_IncLevel, RI_IncLevel, A3SS_IncLevel, A5SS_IncLevel))
saveRDS(IncLevel, paste0("rmats_results/",comparison,"/tables/PSI.rds"))
