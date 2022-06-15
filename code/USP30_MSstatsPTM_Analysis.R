
library(tidyverse)
library(MSstatsPTM)
library(data.table)
library(gridExtra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load(file = "../data/MSstatsPTM_USP30_model.rda")#model_ptm
load(file = "../data/MSstatsPTM_USP30_summary.rda")#summarized_ptm


## Calc # Prot and # pep
summarized_ptm$PTM$FeatureLevelData$global_prot = sapply(summarized_ptm$PTM$FeatureLevelData$PROTEIN, function(x){str_split(x, "_")[[1]][1]})
summarized_ptm$PTM$FeatureLevelData$global_prot %>% n_distinct()
summarized_ptm$PTM$FeatureLevelData$PROTEIN %>% n_distinct()
summarized_ptm$PTM$FeatureLevelData$PEPTIDE %>% n_distinct()

summarized_ptm$PROTEIN$FeatureLevelData$PROTEIN %>% n_distinct()
summarized_ptm$PROTEIN$FeatureLevelData$FEATURE  %>% n_distinct()

## SE analyis
pyst_model

model_ptm$PTM.Model %>% ggplot + geom_boxplot(aes(y = SE))
model_ptm$PROTEIN.Model %>% ggplot + geom_boxplot(aes(y = SE))

model_ptm$PTM.Model %>% summarize(mean(SE, na.rm= TRUE))
model_ptm$PROTEIN.Model %>% summarize(mean(SE, na.rm= TRUE))

##

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
pot <- combo %>% arrange(desc(abs_diff)) %>% head(1000) %>% top_n(1000)


for (i in seq_len(nrow(pot))){

  temp <- summarized_ptm$PTM$FeatureLevelData %>% filter(PROTEIN == pot[i, Protein])
  num <- length(unique(temp$FEATURE))

  if (num >= 2){
    temp_p <- summarized_ptm$PROTEIN$FeatureLevelData %>%
      filter(PROTEIN == pot[i, GlobalProtein])

    num2 <- length(unique(temp_p$FEATURE))

    if (num2 >= 2){
      print(num)
      print(num2)
      print(pot[i, Protein])
    }
  }
}


ptm_list1 <- c("P52209_K059")

features <- summarized_ptm$PTM$FeatureLevelData %>% filter(PROTEIN == ptm_list1[[1]])
features$FeatureType <- "Peptide"
features$Protein <- features$PROTEIN
features$Abundance <- features$ABUNDANCE

summarized <- summarized_ptm$PTM$ProteinLevelData %>% filter(Protein == ptm_list1[[1]])
summarized$Abundance <- summarized$LogIntensities
summarized$FeatureType <- "Model"

test <- rbindlist(list(features, summarized), fill = TRUE)

test[test$FeatureType == 'Model'][['FeatureType']] <- "PTM Summarized"
test[test$FeatureType == 'Peptide'][['FeatureType']] <- "PTM Feature"

p1 <- test %>% ggplot() +
  geom_line(aes(x = originalRUN, y = Abundance , group = FEATURE, color = FeatureType), size =2) +
  geom_point(aes(x = originalRUN, y = Abundance , group = FEATURE, color = FeatureType), size = 5) +
  geom_vline(data=data.frame(x = c(4.5, 8.5, 12.5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed") +
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) +
  #scale_size_manual(values = c(1, 2)) +
  labs(title = "P52209 with ubiquitination at site K059", x = "BioReplicate", y = "Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        title = element_text(size = 22),
        strip.text = element_text(size = 16),
        legend.title =  element_blank(),
        legend.direction = "horizontal",
        legend.position = c(.5, .05)) +
  annotate("text", x = 2.5, y = 28, label = "CCCP", size = 8) +
  annotate("text", x = 6.5, y = 28, label = "Combo", size = 8) +
  annotate("text", x = 10.5, y = 28, label = "Ctrl", size = 8) +
  annotate("text", x = 14.5, y = 28, label = "USP30", size = 8) +
  ylim(19, 28)

p1

prot <- 'P52209'

features <- summarized_ptm$PROTEIN$FeatureLevelData %>% filter(PROTEIN == prot)
features$FeatureType <- "Peptide"
features$Protein <- features$PROTEIN
features$Abundance <- features$ABUNDANCE

temp <- data.frame(originalRUN = c("CCCP-B1T1", "CCCP-B1T2", "CCCP-B2T1" , "CCCP-B2T2",
                                   "Combo-B1T1", "Combo-B1T2", "Combo-B2T1" , "Combo-B2T2",
                                   "Ctrl-B1T1", "Ctrl-B1T2", "Ctrl-B2T1" , "Ctrl-B2T2",
                                   "USP30_OE-B1T1", "USP30_OE-B1T2", "USP30_OE-B2T1" , "USP30_OE-B2T2"),
                   Abundance = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
temp$FeatureType <- "Model"
temp$FEATURE <- "temp2"
features <- features %>% filter(censored == FALSE)

summarized <- summarized_ptm$PROTEIN$ProteinLevelData %>% filter(Protein == prot)
summarized$Abundance <- summarized$LogIntensities
summarized$FeatureType <- "Model"

test <- rbindlist(list(features, summarized,temp), fill = TRUE)

test[test$FeatureType == 'Model'][['FeatureType']] <- "Protein Summarized"
test[test$FeatureType == 'Peptide'][['FeatureType']] <- "Protein Feature"
#test <- test %>% filter(!is.na(FeatureType))

p2 <- test %>% ggplot() +
  geom_line(aes(x = originalRUN, y = Abundance , group = FEATURE, 
                color = FeatureType),  size = 2) +
  geom_point(aes(x = originalRUN, y = Abundance , group = FEATURE, 
                 color = FeatureType), size = 5) +
  geom_vline(data=data.frame(x = c(4.5, 8.5, 12.5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed") +
  scale_colour_manual(values = c("#C3C3C3", "#D55E00", "#D55E00")) +
  scale_size_manual(values = c(1, 2,3)) +
  labs(title = "P52209 Protein", x = "BioReplicate", y = "Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust=1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        title = element_text(size = 22),
        strip.text = element_text(size = 16),
        legend.title =  element_blank(),
        legend.direction = "horizontal",
        legend.position = c(.5, .05)) +
  annotate("text", x = 2.5, y = 28, label = "CCCP", size = 8) +
  annotate("text", x = 6.5, y = 28, label = "Combo", size = 8) +
  annotate("text", x = 10.5, y = 28, label = "Ctrl", size = 8) +
  annotate("text", x = 14.5, y = 28, label = "USP30", size = 8) +
  ylim(19, 28)
p2
grid.arrange(p1, p2, nrow = 1)

## Venn Diagramm
unadj_ptm_model <- model_ptm$PTM.Model
adj_ptm_model <- model_ptm$ADJUSTED.Model

unadj_ptm_model$Protein_Label <- paste(unadj_ptm_model$Protein, unadj_ptm_model$Label)
adj_ptm_model$Protein_Label <- paste(adj_ptm_model$Protein, adj_ptm_model$Label)

sig_unadj_ptm_model <- unadj_ptm_model %>% filter(Protein_Label %in% adj_ptm_model$Protein_Label)

sig_unadj_ptm_model <- sig_unadj_ptm_model %>% filter(adj.pvalue < .05)
sig_adj_ptm_model <- adj_ptm_model %>% filter(adj.pvalue < .05)

colors <- c("#E69F00", "#56B4E9")

venn.diagram(
  x = list(sig_unadj_ptm_model$Protein_Label, sig_adj_ptm_model$Protein_Label),
  category.names = c("Unadjusted" , "Adjusted"),
  filename = "usp30_venn_diagramm.png",
  output=TRUE,
  imagetype="png" ,
  height = 1400,
  width = 1400,
  resolution = 100,
  lwd = 2,
  fill = colors,
  main.fontface = "bold",
  main.fontfamily="sans",
  main.pos = c(.5,1),
  fontface = "bold",
  cat.fontface = "bold",
  cex = 3,
  cat.cex = 3,
  cat.pos = c(-40, 30),
  cat.dist = c(.037, .03),
  main.cex = 2.8,
  main = "Overlap between significant adjusted and unadjusted PTMs"
)

sig_unadj_ptm_model <- sig_unadj_ptm_model %>% filter(adj.pvalue < .05 & is.finite(log2FC))
sig_adj_ptm_model <- adj_ptm_model %>% filter(adj.pvalue < .05 & is.finite(log2FC))

venn.diagram(
  x = list(sig_unadj_ptm_model$Protein_Label, sig_adj_ptm_model$Protein_Label),
  category.names = c("Unadjusted" , "Adjusted"),
  filename = "usp30_venn_diagramm_matching_only.png",
  output=TRUE,
  imagetype="png" ,
  height = 1400,
  width = 1400,
  resolution = 100,
  lwd = 2,
  fill = colors,
  main.fontface = "bold",
  main.fontfamily="sans",
  main.pos = c(.5,1),
  fontface = "bold",
  cat.fontface = "bold",
  cex = 3,
  cat.cex = 3,
  cat.pos = c(-40, 30),
  cat.dist = c(.037, .03),
  main.cex = 2.8,
  main = "Signficant adjusted and unadjusted PTMs (matching only)"
)

combo <- adj_ptm_model %>% merge(unadj_ptm_model, by = "Protein_Label", all.x = TRUE, all.y = TRUE)

combo %>% filter(log2FC.x < log2FC.y*1.09 & log2FC.x > log2FC.y*.91 &
                   adj.pvalue.x >= .05 & adj.pvalue.y < .05)




