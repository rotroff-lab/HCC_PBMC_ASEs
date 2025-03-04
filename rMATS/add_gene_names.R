#Adding Names to rMATS Gene IDs
library(dplyr)
library(tidyr)
library(GenomicRanges)

comparison <- as.character(commandArgs(TRUE[1]))


rmats <- read.table(paste0("rmats_results/", comparison, "/tables/all_rmats_inc.bed"))
colnames(rmats) <- c("chr", "start", "end", "ID", "score", "strand")
rmats$type <- sub(",.*", "", rmats$ID)
rmatsbed <- rmats


##refseq instead of ucsc
refseq_genes <- read.delim("data/refseq_hg38_exon_count.txt", header=FALSE)

colnames(refseq_genes) <- c("ID","chr","strand","start","end", "exonCount", "gene_name")

##remove genes with only 1 exon
refseq_genes <- refseq_genes %>% filter(exonCount > 1) %>%
  select(-exonCount)



#Assigning column types
rmatsbed$start <- as.integer(as.character(rmatsbed$start))
rmatsbed$end <- as.integer(as.character(rmatsbed$end))




#Make GRanges object to be used with UCSC GRanges objects
rMATScoord <- makeGRangesFromDataFrame(rmats,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=FALSE,
                                       seqinfo=NULL,
                                       seqnames.field="chr", 
                                       start.field="start",
                                       end.field="end",
                                       strand.field="strand",
                                       starts.in.df.are.0based=TRUE)

UCSCcoord <- makeGRangesFromDataFrame(refseq_genes,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=FALSE,
                                      seqinfo=NULL,
                                      seqnames.field="chr", 
                                      start.field="start",
                                      end.field="end",
                                      strand.field="strand",
                                      starts.in.df.are.0based=TRUE)


ranges <- subsetByOverlaps(rMATScoord, UCSCcoord, type = "within")

hits <- findOverlaps(rMATScoord, UCSCcoord,
                     maxgap=0L, minoverlap=0L,
                     type="within",
                     select=c("all"),
                     ignore.strand=FALSE)

rmats_hits <- rmatsbed[queryHits(hits),]
genenames <- refseq_genes[subjectHits(hits),]

ID_names <- cbind(rmats_hits,genenames)
colnames(ID_names) <- c("chr_AS","start_AS","end_AS","AS_ID","score_AS","strand_AS","type","gene_ID","chr_gene","strand_gene","start_gene","end_gene","gene_name")
ID_names <- ID_names[which(ID_names$strand_AS==ID_names$strand_gene),]

#####add brackets around genes with a dash
ID_names <- ID_names %>%
  mutate(gene_name = if_else(grepl("-", gene_name), paste("[", gene_name, "]", sep = ""), gene_name))


ID_names <-ID_names %>%
  group_by(AS_ID) %>%
  dplyr::summarize(gene_name = paste(unique(gene_name), collapse = "."))

ID_names <- ID_names[,c("AS_ID", "gene_name")]


#Merge list of gene names with rMATS table
rmats_named <- merge(ID_names, rmats, by.x="AS_ID", by.y="ID", all.y=T)
rmats_named$ID <- rmats_named$AS_ID


##remove concatanted read-through genes
rmats_named <- rmats_named %>%
  mutate(gene_name = case_when(
    grepl("\\.", gene_name) ~ if_else(grepl("\\[.*?\\]", gene_name), gsub("\\.|\\[.*?\\]", "", gene_name), gene_name),
    TRUE ~ gene_name))

##Remove the brackets from dashed genes. Originally used for subsetting
rmats_named <- rmats_named %>%
  mutate(gene_name = gsub("\\[|\\]", "", gene_name))

rmats_named <- unique(rmats_named)
###Write tables!

#remove unnamed
rmats_named <- rmats_named[!is.na(rmats_named$gene_name),]

saveRDS(rmats_named,paste0("rmats_results/", comparison, "/tables/splice_ID_named_refSeq.rds"))


