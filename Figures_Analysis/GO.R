##Go enrichment of significant negatively correlated intron/gene expression
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(purrr)
library(org.Hs.eg.db)


results <- readRDS("meta_analysis_ASE.rds")
results$direction <- ifelse(results$meanincleveldiff > 0, "up", "down")
results_sig <- results %>% filter(abs(meanincleveldiff) >= 0.05 & fdr <= 0.05)

results_up <- results_sig %>% filter(direction =="up")
results_down <- results_sig %>% filter(direction =="down")

genes_up <- results_up %>% pull(gene_name) %>% unique()
genes_down <- results_down %>% pull(gene_name) %>% unique()

##Go terms
GO_ont <- c("BP", "MF", "CC")

##dot color pallete for treeplot
my_palette <- c( "#000B18", "#e0e0e0")
my_color_scale <- scale_color_gradient(low = my_palette[1], high = my_palette[length(my_palette)])

exon_type <- c("SE", "A3SS", "A5SS", "RI")
results_all_exon <- results_sig %>% filter(type %in% exon_type & annotated == "annotated")
all_genes <- results_all_exon %>% pull(gene_name) %>% unique()

gene_list <- list(all_genes)

GO_obj <- enrichGO(gene          = unlist(gene_list[[1]]), 
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

print(treeplot(pairwise_termsim(GO_obj), cluster.params = list(color = c("#000B18", "#00246d", "#32527b","#B8D4FF", "#e0e0e0")), showCategory = 60) + my_color_scale + ggtitle("BP"))
ggsave("GO_sig_ASEs.pdf", width = 10, height = 10, units = "in", dpi = 1000)
write.xlsx(GO_obj@result, "figures/Table_S2.xlsx")
