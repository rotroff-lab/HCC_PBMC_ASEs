library(ggplot2)
library(dplyr)
library(ggthemes)
library(hrbrthemes)
library(waffle)

names(results_test)[1:2] <- c("pred", "diagnosis")
names(results_validation)[1:2] <- c("pred", "diagnosis")


results_validation$group <- "validation"
results_test$group <- "test"
results_train$group <- "train"

results_test$RUN <- rownames(results_test)
results_validation$RUN <- rownames(results_validation)
names(results_train)[1] <- "RUN"

all_results <- rbind(results_train[,c("RUN","pred","diagnosis","group")], results_test[,c("RUN","pred","diagnosis","group")]) %>% rbind(., results_validation[,c("RUN","pred","diagnosis","group")])


##add staging to predictions
metadata_all <- readRDS("metadata.rds") 
metadata_all$BCLC_STAGE <- ifelse(metadata_all$GROUP == "NoHCC", "Healthy", metadata_all$BCLC_STAGE)
metadata_all$BCLC_STAGE <- ifelse(is.na(metadata_all$TNM_STAGE), metadata_all$BCLC_STAGE, metadata_all$TNM_STAGE)
metadata_all$BCLC_STAGE <- ifelse(metadata_all$GROUP == "HCC" & is.na(metadata_all$BCLC_STAGE), "U", metadata_all$BCLC_STAGE)

all_results <- left_join(all_results, metadata_all, by="RUN")

#summarize predicitons and staging for waffle plot
count_train <- all_results %>%
  dplyr::count(diagnosis, pred, group, BCLC_STAGE) %>%
  mutate(
    pred = recode(pred, "group1" = "HCC", "group2" = "Healthy", "group3" = "HBV"),
    diagnosis = recode(diagnosis, "group1" = "HCC", "group2" = "Healthy", "group3" = "HBV")
  )
count_train$group <- factor(count_train$group, levels = c("train", "test", "validation"))



##Order staging by prediciton and truth for waffle plot annotations
count_train2 <- count_train %>% arrange(pred)
count_train2$BCLC_STAGE <- gsub("b", "", count_train2$BCLC_STAGE)
count_train2$BCLC_STAGE <- ifelse(count_train2$BCLC_STAGE == "Healthy" | count_train2$BCLC_STAGE == "HBV" , "", count_train2$BCLC_STAGE)
count_train2$BCLC_STAGE <- factor(count_train2$BCLC_STAGE, levels = c("A0", "A", "I", "B", "II", "C", "III", "D", "IV", "U", ""))
count_train2 <- count_train2 %>% arrange(group, pred, BCLC_STAGE)
count_train2$BCLC_STAGE <- as.character(count_train2$BCLC_STAGE)
count_train2$group <- as.character(count_train2$group)
count_train2$diagnosis <- as.character(count_train2$diagnosis)
count_train2$pred <- as.character(count_train2$pred)

new_df <- data.frame(
  diagnosis = character(),
  pred = character(),
  group = character(),
  BCLC_STAGE = character(),
  x = numeric(),
  y = numeric(),
  stringsAsFactors = FALSE
)



for (i in 1:nrow(count_train2)) {
  for (j in 1:count_train2$n[i]) {
    # Add a new row to new_df
    new_df[nrow(new_df) + 1, ] <- count_train2[i, ]
  }
}

new_df$group <- factor(new_df$group, levels = c("train", "test", "validation"))
new_df$pred <- factor(new_df$pred, levels = c("HCC", "Healthy", "HBV"))

numbers <- 1:10
repeated_numbers <- rep(numbers, each = 10)
repeated_numbers_y <- rep(numbers, times = 10)

new_df_1 <- new_df %>% filter(pred == "HCC" & group == "train")
new_df_2 <- new_df %>% filter(pred == "Healthy" & group == "train")
new_df_3 <- new_df %>% filter(pred == "HBV" & group == "train")

new_df_4 <- new_df %>% filter(pred == "HCC" & group == "test")
new_df_5 <- new_df %>% filter(pred == "Healthy" & group == "test")
new_df_6 <- new_df %>% filter(pred == "HBV" & group == "test")

new_df_7 <- new_df %>% filter(pred == "HCC" & group == "validation")
new_df_8 <- new_df %>% filter(pred == "Healthy" & group == "validation")
new_df_9 <- new_df %>% filter(pred == "HBV" & group == "validation")


new_df_1$x <- repeated_numbers[1:nrow(new_df_1)]
new_df_2$x <- repeated_numbers[1:nrow(new_df_2)]
new_df_3$x <- repeated_numbers[1:nrow(new_df_3)]
new_df_4$x <- repeated_numbers[1:nrow(new_df_4)]
new_df_5$x <- repeated_numbers[1:nrow(new_df_5)]
new_df_6$x <- repeated_numbers[1:nrow(new_df_6)]
new_df_7$x <- repeated_numbers[1:nrow(new_df_7)]
new_df_8$x <- repeated_numbers[1:nrow(new_df_8)]

if (nrow(new_df_9)>0) {
new_df_9$x <- repeated_numbers[1:nrow(new_df_9)]
}

new_df_1$y <- repeated_numbers_y[1:nrow(new_df_1)]
new_df_2$y <- repeated_numbers_y[1:nrow(new_df_2)]
new_df_3$y <- repeated_numbers_y[1:nrow(new_df_3)]
new_df_4$y <- repeated_numbers_y[1:nrow(new_df_4)]
new_df_5$y <- repeated_numbers_y[1:nrow(new_df_5)]
new_df_6$y <- repeated_numbers_y[1:nrow(new_df_6)]
new_df_7$y <- repeated_numbers_y[1:nrow(new_df_7)]
new_df_8$y <- repeated_numbers_y[1:nrow(new_df_8)]
if (nrow(new_df_9)>0) {
new_df_9$y <- repeated_numbers_y[1:nrow(new_df_9)]
}



diagnosis_colors <- c(
  "HCC" = "#880808", 
  "Healthy" = "lightblue", 
  "HBV" = "#097969"
)

#generate waffle plot
p <- ggplot(data = count_train, aes(fill = diagnosis)) +
  geom_waffle(aes(values = n), color = "white", size = 1.125, n_rows = 10) +
  facet_wrap(group~pred, nrow = 1) +
  scale_x_discrete(expand = c(0, 0, 0, 0)) +
  scale_y_discrete(expand = c(0, 0, 0, 0)) +
  ggthemes::scale_fill_tableau(name = NULL) +
  coord_equal() +
  theme_ipsum_rc(grid = "", base_family = "Helvetica") +
  scale_fill_manual(
    name = NULL,
    values = diagnosis_colors
  )


##iterate through classificationn rows and add stage label
add_annotations_to_plot <- function(plot, df1, df2) {
  for (i in 1:nrow(df1)) {
    df_subset <- as.data.frame(df1[i,])
    plot <- plot + geom_text(data = df_subset, 
                             aes(x = x, y = y, label = BCLC_STAGE), 
                             size = 4, colour = "white", fontface = "bold")
  }
  
  for (i in 1:nrow(df2)) {
    df_subset <- as.data.frame(df2[i,])
    plot <- plot + geom_text(data = df_subset, 
                             aes(x = x, y = y, label = BCLC_STAGE), 
                             size = 4, colour = "white", fontface = "bold")
  }
  
  plot <- plot + theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank())
  
  return(plot)
}


p <- add_annotations_to_plot(p, new_df_1, new_df_2)
p <- add_annotations_to_plot(p, new_df_3, new_df_4)
p <- add_annotations_to_plot(p, new_df_5, new_df_6)
p <- add_annotations_to_plot(p, new_df_7, new_df_8)
if (nrow(new_df_9)>0) {
p <- add_annotations_to_plot(p, new_df_9, new_df_9)
}
p
