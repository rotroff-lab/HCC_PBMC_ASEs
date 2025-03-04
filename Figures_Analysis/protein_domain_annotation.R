library(dplyr)
library(ggplot2)
library(ggforce)

## Read annotated splicing events
ann_exons <- read.delim("annotated_exons_PBMC.txt")

## Define gene types
genetype <- c("other", "immune")

## Loop through each gene type
for (gene in genetype) {
  
  ## Initialize data frames
  combined_df_tier1 <- data.frame()
  combined_df_tier2 <- data.frame()
  
  ## Import significant splicing dataset
  all_splicing <- readRDS("meta_analysis_ASE.rds")
  all_splicing <- all_splicing %>% filter(annotated == "annotated")
  
  ## Read gene names associated with immune category
  GO0002376_gene_names <- read_excel("GO0002376_gene_names.xlsx", col_names = FALSE)
  names(GO0002376_gene_names) <- "gene_name"   
  GO0002376_gene_names$gene_type <- "immune"
  
  ## Merge with all splicing data
  results_table_immune <- left_join(all_splicing, GO0002376_gene_names, by = "gene_name")
  results_table_immune$gene_type <- ifelse(is.na(results_table_immune$gene_type), "other", "immune")
  results_table_immune1 <- results_table_immune %>% filter(gene_type == "immune")
  results_table_immune2 <- results_table_immune %>% filter(gene_type == "other")
  all_splicing <- rbind(results_table_immune1, results_table_immune2)
  
  ## Filter for significant events
  sig_splicing <- all_splicing %>% filter(abs(meanincleveldiff) >= 0.05 & fdr <= 0.05)
  
  ## Filter for all measured AS events
  all_measured_AS_events <- all_splicing
  
  if (gene == "immune") {
    sig_splicing <- sig_splicing %>% filter(gene_type == "immune")
    all_measured_AS_events <- all_measured_AS_events %>% filter(gene_type == "immune")
  }
  
  ## Remove RI events
  sig_splicing_no_RI <- sig_splicing[!grepl("^RI", sig_splicing$AS_ID_short), ]
  all_measured_AS_events_no_RI <- all_measured_AS_events[!grepl("^RI", all_measured_AS_events$AS_ID_short), ]
  
  ## Anti-join to remove significant AS events from all measured AS events
  all_measured_AS_events_no_RI <- anti_join(all_measured_AS_events_no_RI, sig_splicing_no_RI, by = "AS_ID") 
  
  ## Merge with annotated exons
  sig_annotated_splicing <- left_join(sig_splicing_no_RI, ann_exons, by = c("ID" = "AS_ID"))
  all_measured_AS_events_no_RI_annotated <- left_join(all_measured_AS_events_no_RI, ann_exons, by = c("ID" = "AS_ID"))
  
  ## Get unique domain IDs
  sig_annotated_splicing_unique <- sig_annotated_splicing %>% distinct(AS_ID, x.name, .keep_all = TRUE)
  all_measured_AS_events_no_RI_annotated_unique <- all_measured_AS_events_no_RI_annotated %>% distinct(AS_ID, x.name, .keep_all = TRUE)
  
  ## Get list of domain IDs to be tested
  domain_ids <- sig_annotated_splicing_unique %>% distinct(x.name) %>% pull(x.name)
  domain_ids <- na.omit(domain_ids)
  
  result_df <- data.frame()
  contingency_table_list <- list()
  
  ## Loop through each domain ID
  for (i in domain_ids) {
    
    ## Calculate contingency table
    total_SIG_AS_in_DOMAIN <- nrow(sig_annotated_splicing_unique %>% filter(x.name == i))
    total_SIG_AS_not_in_DOMAIN <- nrow(sig_annotated_splicing_unique %>% filter(x.name != i | is.na(x.name)))
    total_MEASURED_AS_in_DOMAIN <- nrow(all_measured_AS_events_no_RI_annotated_unique %>% filter(x.name == i))
    total_MEASURED_AS_not_in_DOMAIN <- nrow(all_measured_AS_events_no_RI_annotated_unique %>% filter(x.name != i | is.na(x.name)))
    
    contingency_table <- matrix(c(
      total_SIG_AS_in_DOMAIN,
      total_MEASURED_AS_in_DOMAIN,
      total_SIG_AS_not_in_DOMAIN,
      total_MEASURED_AS_not_in_DOMAIN
    ), nrow = 2, byrow = TRUE)
    
    ## Perform Fisher's Exact Test
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    
    ## Calculate ExonRatio
    ExonRatio <- total_SIG_AS_in_DOMAIN / nrow(sig_annotated_splicing_unique)
    
    ## Get major domain category name
    temp_df <- sig_annotated_splicing_unique %>% filter(x.name == i)
    ancestor_name <- unique(temp_df$ancestor_name)
    
    ## Create new row with results
    new_result <- data.frame(
      domain = i,
      SIG_AS_in_domian = total_SIG_AS_in_DOMAIN,
      SIG_AS_not_in_DOMAIN = total_SIG_AS_not_in_DOMAIN,
      NON_SIG_AS_in_DOMAIN = total_MEASURED_AS_in_DOMAIN,
      NON_SIG_AS_not_in_DOMAIN = total_MEASURED_AS_not_in_DOMAIN,
      major_domain = ancestor_name,
      odds_ratio = fisher_result[[3]],
      p_val = fisher_result[[1]]
    )
    
    ## Add adjusted p-value using BH procedure
    new_result$adj_p_val <- p.adjust(new_result$p_val, method = "BH")
    new_result$log10_P <- -log10(new_result$adj_p_val)
    
    ## Append to result dataframe
    result_df <- rbind(result_df, new_result)
    contingency_table_list[[i]] <- contingency_table
  }
  
  ## Set column names for result dataframe
  colnames(result_df) <- c("domain", "SIG_AS_in_domian", "SIG_AS_not_in_DOMAIN", "NON_SIG_AS_in_DOMAIN", "NON_SIG_AS_not_in_DOMAIN", "major_domain", "odds_ratio", "p-val", "adj_p-val", "log10_P")
  
  ## Count significantly over-represented domains
  n_significant_domains <- nrow(result_df %>% filter(adj_p_val <= 0.05))
  
  ## Bar plot depicting number of AS events per protein domain
  sig_annotated_splicing_unique$direction <- ifelse(sig_annotated_splicing_unique$meanincleveldiff > 0, "up", "down")
  sig_annotated_splicing_unique <- sig_annotated_splicing_unique %>% distinct(AS_ID, ancestor_name, .keep_all = TRUE)
  df_domain_direction <- as.data.frame(table(sig_annotated_splicing_unique[, c("ancestor_name", "direction")]))
  
  ## Plotting
  p <- ggplot(df_domain_direction, aes(x = reorder(ancestor_name, -Freq), y = Freq, fill = direction)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = c("#b7c9e2", "#0D043C")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), text = element_text(size = 16)) +
    xlab(label = "Protein Major Domains") +
    ylab(label = "Number of AS Events Encompassing Protein Domain")
  
  ## Save plot
  ggsave(paste0(gene, "_protein_domain_ASEs.pdf"), plot = p, width = 6, height = 8, units = "in", dpi = 600)
  
  
}
