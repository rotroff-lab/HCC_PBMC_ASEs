library(tidyr)
library(dplyr)
library(ggplot2)

# Ntree optimization
ntree_opt <- readRDS("ntree_optimization.rds")

ntree_opt <- ntree_opt %>%
  select(
    contains("ntreeRF"),
    contains("balanced_accuracy"),
    contains("sensitivity"),
    contains("specificity"),
    contains("misclassification")
  )

ntree_opt <- ntree_opt %>%
  gather(key = "classification", value = "value", -ntreeRF)
ntree_opt$group <- sub(".*_", "", ntree_opt$classification)
ntree_opt$classification2 <- sub("_[^_]*$", "", ntree_opt$classification)


ggplot(ntree_opt, aes(x = ntreeRF, y = value, group = group)) +
  geom_line(aes(color = group)) +
  geom_point(aes(color = group)) +
  scale_x_continuous(breaks = c(100, 200, 300, 400)) +
  facet_wrap(~classification2, scales = "free") +
  scale_color_manual(values = c("#097969", "#880808", "lightblue")) +
  theme_classic()
ggsave("ntree_optimization.pdf", width = 7, height = 4)



# Highlight features
features <- readRDS("ml_model_features.rds")
features <- rownames(features$importance)
gini <- readRDS("gini_all_features.rds")
rownames(gini) <- gini$feature
gini$feature <- NULL
gini <- as.matrix(gini)
gini_plot <- data.frame(feature=rownames(gini), min=rowMins(gini), max= rowMaxs(gini), median= rowMedians(gini))
gini_plot <- gini_plot[order(gini_plot$median, decreasing = T),]
gini_plot$feature <- factor(gini_plot$feature, levels=rev(gini_plot$feature))
gini_highlight <- gini_plot %>% filter(feature %in% features)

ggplot(gini_plot, aes(x = feature, y = median)) +
  geom_errorbar(aes(ymin = min, ymax = max), color = "darkgrey", width = 0.01, size = 0.05) +
  geom_errorbar(aes(ymin = min, ymax = max), color = "darkgrey", width = 0.01, size = 0.05) +
  geom_point(size = 0.01) +
  theme_classic() +
  geom_pointrange(data = gini_highlight, aes(x = feature, y = median, ymin = min, ymax = max), color = "red", size=0.01) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
ggsave("all_features_gini_score.pdf", width = 10, height = 4)  





# Highlight ml model feature OOB mean
results_RF <- readRDS("iRF_LOOCV_ntree200.rds")

optimal_features <- results_RF %>% slice_min(OOB_mean, n = 1) %>% pull(num_features)
results_RF_highlight <- results_RF %>% filter(num_features == optimal_features)


ggplot(results_RF, aes(x = num_features, y = OOB_mean)) +
  geom_errorbar(aes(ymin = OOB_min, ymax = OOB_max), color = "darkgrey", width = 0.01, size = 0.05) +
  geom_errorbar(aes(ymin = OOB_min, ymax = OOB_mean), color = "darkgrey", width = 0.01, size = 0.05) +
  geom_point(size = 0.5) +
  theme_classic() +
  geom_point(data = subset(results_RF_highlight), color = "red", size = 1) +
  geom_errorbar(data = subset(results_RF_highlight), aes(ymin = OOB_min, ymax = OOB_max), color = "red", width = 0.01, size = 0.5) +
  geom_errorbar(data = subset(results_RF_highlight), aes(ymin = OOB_min, ymax = OOB_mean), color = "red", width = 0.01, size = 0.5) +
  scale_x_reverse()
ggsave("8_features.pdf", width = 10, height = 4)  
