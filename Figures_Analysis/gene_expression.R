library(dplyr)
library(edgeR)
library(limma)
library(MetaVolcanoR)
library(org.Hs.eg.db)

# Load metadata
metadata <- readRDS("metadata.rds")

# Get unique SRP studies
SRP_list <- unique(metadata$SRA.STUDY)

# Standardize group labels
metadata$GROUP <- ifelse(metadata$GROUP == "HCC", "group1", "group2")

# Function to process each SRP dataset
eachSRP <- function(SRP, metadata) {
  countdata_path <- paste0("gene_counts/", SRP, "_counts.txt")
  annot_path <- paste0("gene_counts/", SRP, "_annotation.txt")
  
  # Read data
  countdata <- readRDS(countdata_path)
  annot <- readRDS(annot_path)
  
  # Convert ENTREZ IDs to SYMBOLs
  countdata_names <- bitr(rownames(countdata), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop=FALSE)
  countdata_names$SYMBOL <- ifelse(is.na(countdata_names$SYMBOL), countdata_names$ENTREZID, countdata_names$SYMBOL)
  rownames(countdata) <- countdata_names$SYMBOL
  
  # Clean column names
  colnames(countdata) <- sub("\\.bam", "", colnames(countdata))
  
  # Subset metadata
  coldata <- metadata %>%
    filter(SRA.STUDY == SRP) %>%
    select(RUN, GROUP) %>%
    unique() %>%
    complete.cases() %>%
    as.data.frame()
  
  colnames(coldata) <- c("SampleName", "group")
  rownames(coldata) <- coldata$SampleName
  countdata <- countdata[, rownames(coldata)]
  
  # Differential expression analysis using limma
  tryCatch({
    expfactor <- factor(coldata$group, levels = c("group2", "group1"))
    design <- model.matrix(~expfactor)
    colnames(design) <- c("nonHCC", "HCC_vs_nonHCC")
    
    dge <- DGEList(counts=countdata)
    keep <- filterByExpr(dge, design)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend=TRUE)
    diffexp <- topTable(fit, coef="HCC_vs_nonHCC", confint=TRUE, number=100000)
    
    diffexp$Symbol <- rownames(diffexp)
    write.table(diffexp, paste0("gene_counts/", SRP, "_limma_res.txt"))
    return(diffexp)
  }, error=function(e) {
    message("Error processing ", SRP, ": ", e$message)
    return(NULL)
  })
}

# Apply function to all SRP datasets
results_list <- lapply(SRP_list, eachSRP, metadata=metadata)
names(results_list) <- SRP_list

# Save results
saveRDS(results_list, "gene_counts/all_PBMC_results.rds")

# Meta-analysis using MetaVolcanoR
diffexp_list <- readRDS("gene_counts/all_PBMC_results.rds")

meta_degs_rem <- rem_mv(
  diffexp=diffexp_list,
  pcriteria="P.Value",
  foldchangecol='logFC', 
  genenamecol="Symbol",
  geneidcol=NULL,
  collaps=FALSE,
  llcol='CI.L',
  rlcol='CI.R',
  vcol=NULL, 
  cvar=TRUE,
  metathr=0.01,
  draw='HTML',
  ncores=4
)

# Extract meta-analysis results and adjust FDR
meta_analysis <- meta_degs_rem@metaresult
meta_analysis$FDR <- p.adjust(meta_analysis$randomP, method="fdr")

# Save meta-analysis results
saveRDS(meta_analysis, "gene_counts/pbmc_gene_meta.rds")
saveRDS(meta_degs_rem, "gene_counts/pbmc_gene_meta_object.rds")

