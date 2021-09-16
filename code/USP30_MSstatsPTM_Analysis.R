
library(tidyverse)
library(MSstatsPTM)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(file = "../data/ups30_msstatsptm_model.rda")#model_ptm
load(file = "../data/ups30_msstatsptm_summary.rda")#summarized_ptm


original_ptms <- unique(model_ptm$PTM.Model$Protein)
adj_ptms <- unique(model_ptm$ADJUSTED.Model$Protein)
global_prot <- unique(model_ptm$PROTEIN.Model$Protein)

match_list <- c()
for (i in seq_along(global_prot)){
  for (y in seq_along(original_ptms))
    if (grepl(global_prot[[i]], original_ptms[[y]]) & 
        !global_prot[[i]] %in% match_list) {
      match_list <- c(match_list, global_prot[[i]])
    }
}

data.frame("Number" = c(original_ptms, adj_ptms, global_prot), Label = c("Total PTMs", "Matching PTMs", "Global Proteins")) %>% 
  ggplot() + geom_col(aes(y = Number, x = Label))

prot <- sapply(original_ptms, function(x){substring(str_extract_all(x, "^(.*?)_")[[1]], 
                                                    1, nchar(str_extract_all(x, "^(.*?)_")[[1]]) - 1)})
length(intersect(unique(prot), global_prot))

venn.diagram(
  x = list(prot, global_prot),
  category.names = c("PTMs" , "Global Proteins"),
  filename = "test.png",
  output=FALSE
)


## Profile Plots
combo <- merge(model_ptm$ADJUSTED.Model, model_ptm$PTM.Model, all.x = TRUE, by = c("Protein", "Label"))

combo$abs_diff <- abs(combo$log2FC.x - combo$log2FC.y)
pot <- (combo %>% arrange(desc(abs_diff)) %>% head(100) %>% top_n(100) %>% select(Protein))[[1]]


for (i in seq_along(pot)){
  
  temp <- summarized_ptm$PTM$FeatureLevelData %>% filter(PROTEIN == pot[[i]])
  num <- length(unique(temp$FEATURE))
  if (num >= 2){
    print(num)
    print(pot[[i]])
  }
}


ptm_list1 <- c("Q16186_K034")

# pSTY_features$BioReplicate <- factor(pSTY_features$BioReplicate, 
#                                      levels = c("WT_Uninfect_1", "WT_Uninfect_2", "WT_Uninfect_3" , "WT_Uninfect_4", 
#                                                 "WT_Early_1", "WT_Early_2", "WT_Early_3", "WT_Early_4",
#                                                 "WT_Late_1", "WT_Late_2", "WT_Late_3", "WT_Late_4", 
#                                                 "KO_Uninfect_1", "KO_Uninfect_2" , "KO_Uninfect_3", 
#                                                 "KO_Early_1", "KO_Early_2", "KO_Early_3", "KO_Early_4",
#                                                 "KO_Late_1", "KO_Late_2", "KO_Late_3", "KO_Late_4"))

features <- summarized_ptm$PTM$FeatureLevelData %>% filter(PROTEIN == ptm_list1[[1]])
features$Feature <- "PSM"
features$Protein <- features$PROTEIN
features$Abundance <- features$ABUNDANCE 

summarized <- summarized_ptm$PTM$ProteinLevelData %>% filter(Protein == ptm_list1[[1]])
summarized$Abundance <- summarized$LogIntensities 
summarized$Feature <- "Summary"

test <- rbindlist(list(features, summarized), fill = TRUE)

p1 <- test %>% ggplot() + 
  geom_line(aes(x = originalRUN, y = Abundance , group = FEATURE, color = Feature)) + #,  size = Feature)
  geom_point(aes(x = originalRUN, y = Abundance , group = FEATURE, color = Feature), size = 3) +
  geom_vline(data=data.frame(x = c(4.5, 8.5, 12.5)), 
             aes(xintercept=as.numeric(x)), linetype = "dashed") + 
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) + 
  scale_size_manual(values = c(1, 2)) + 
  labs(title = ptm_list1[[1]], x = "BioReplicate", y = "Abundance") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12), 
        axis.text.y = element_text(size = 12), 
        legend.text=element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        title = element_text(size = 16),
        strip.text = element_text(size = 12),
        legend.position = "None")  + 
  annotate("text", x = 2.5, y = 29.5, label = "CCCP", size = 4.5) +
  annotate("text", x = 6.5, y = 29.5, label = "Combo", size = 4.5) +
  annotate("text", x = 10.5, y = 29.5, label = "Ctrl", size = 4.5) + 
  annotate("text", x = 14.5, y = 29.5, label = "USP30", size = 4.5) + 
  ylim(23, 29.5)

p1

prot <- 'TTP_MOUSE|P22893'

global_sum_msstats_norm$ProteinName <- global_sum_msstats_norm$Protein
global_sum_msstats_norm$log2Intensity <- global_sum_msstats_norm$Abundance
global_sum_msstats_norm$Feature <- "Summary"

prot_plot <- py_sum$PROTEIN$FeatureLevelData %>% filter(ProteinName == prot)
prot_plot$Feature <- "PSM"
test <- rbindlist(list(prot_plot, global_sum_msstats_norm %>% filter(ProteinName == prot)), 
                  fill = TRUE)

test$BioReplicate <- factor(test$BioReplicate, 
                            levels = c("WT_Uninfect_1", "WT_Uninfect_2", "WT_Uninfect_3" , "WT_Uninfect_4", 
                                       "WT_Early_1", "WT_Early_2", "WT_Early_3", "WT_Early_4",
                                       "WT_Late_1", "WT_Late_2", "WT_Late_3", "WT_Late_4", 
                                       "KO_Uninfect_1", "KO_Uninfect_2" , "KO_Uninfect_3", 
                                       "KO_Early_1", "KO_Early_2", "KO_Early_3", "KO_Early_4",
                                       "KO_Late_1", "KO_Late_2", "KO_Late_3", "KO_Late_4"))

p2 <- test %>% ggplot() +
  geom_line(aes(x = BioReplicate, y = log2Intensity, group = PSM, color = Feature, size =Feature)) + 
  geom_point(aes(x = BioReplicate, y = log2Intensity, group = PSM, color = Feature), size = 3) + 
  geom_vline(data=data.frame(x = c(2.5, 4.5, 6.5, 8.5, 10.5)), 
             aes(xintercept=as.numeric(x)), linetype = "dashed") + 
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) + 
  scale_size_manual(values = c(1, 2)) + 
  labs(title = prot, x = "BioReplicate", y = "Abundance") + 
  facet_grid(Mixture~.) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12), 
        axis.text.y = element_text(size = 12), 
        legend.text=element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        title = element_text(size = 16),
        strip.text = element_text(size = 12),
        legend.position = "None")  + 
  annotate("text", x = 1.5, y = 22.75, label = "WT_Uninfect", size = 4.5) +
  annotate("text", x = 3.5, y = 22.75, label = "WT_Early", size = 4.5) +
  annotate("text", x = 5.5, y = 22.75, label = "WT_Late", size = 4.5) + 
  annotate("text", x = 7.5, y = 22.75, label = "KO_Uninfect", size = 4.5) + 
  annotate("text", x = 9.5, y = 22.75, label = "KO_Early", size = 4.5) + 
  annotate("text", x = 11.5, y = 22.75, label = "KO_Late", size = 4.5) + 
  ylim(9, 23)

grid.arrange(p1, p2, nrow = 1)
