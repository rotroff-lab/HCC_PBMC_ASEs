library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(purrr)
library(org.Hs.eg.db)
library(pathview)
library(dplyr)
library(ggplot2)

# Load results
results <- readRDS("meta_analysis_ASE.rds")

# Add direction column
results <- results %>%
  mutate(direction = ifelse(meanincleveldiff > 0, "up", "down"))

# Filter significant results
results_sig <- results %>%
  filter(abs(meanincleveldiff) >= 0.05 & fdr <= 0.05)

# Filter for specific exon types
exon_types <- c("SE", "A3SS", "A5SS", "RI")
results_all_exon <- results_sig %>%
  filter(type %in% exon_types & annotated == "annotated")

# Get unique gene names and their frequency
all_genes_dup <- as.data.frame(table(results_all_exon$gene_name))
names(all_genes_dup)[1] <- "SYMBOL"
all_genes <- unique(results_all_exon$gene_name)

# Convert gene symbols to ENTREZ IDs
all_genes_entrez <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") %>%
  left_join(all_genes_dup, by = "SYMBOL")

# Perform KEGG enrichment analysis
mkk <- enrichKEGG(
  gene = all_genes_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 1,
  qvalueCutoff = 0.05
)

# Extract KEGG results
mkk_df <- mkk@result
mkk_df_top20 <- mkk_df %>%
  filter(qvalue <= 0.05) %>%
  pull(ID)

# Named list of ENTREZ IDs with ASE frequency
geneList2 <- setNames(all_genes_entrez$Freq, all_genes_entrez$ENTREZID)

# Define color gradient for visualization
gradient <- colorRampPalette(c("lightgrey", "red"))
new_colors <- gradient(3)

# Set output directory
output_dir <- "pathway/kegg/"
dir.create(output_dir, showWarnings = FALSE)
setwd(output_dir)

# Generate pathway visualizations
for (ID in mkk_df_top20) {
  print(paste("Processing pathway:", ID))
  pathview(
    gene.data  = geneList2,
    gene.idtype = "entrez",
    pathway.id = ID,
    species    = "hsa",
    limit      = c(0,3),
    low        = "#D3D3D3",
    mid        = "#E96969",
    high       = "#FF0000"
  )
}

# Create KEGG dotplot
my_palette <- c("#000B18", "#B8D4FF")
my_color_scale <- scale_color_gradient(low = my_palette[1], high = my_palette[2])

dotplot(mkk, showCategory = 30) + my_color_scale
ggsave(filename = "kegg_dotplot.pdf", height = 8, width = 7, dpi = 1000)
