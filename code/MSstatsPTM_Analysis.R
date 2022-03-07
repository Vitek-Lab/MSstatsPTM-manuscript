
library(MSstatsPTM)
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(VennDiagram)

## Load msstats input data
load(file = "pST_MSstatsPTM_Format.rda")#shigella_peptide_pST
load(file = "global_MSstatsPTM_Format.rda")#shigella_global
load(file = "pY_MSstatsPTM_Format.rda")#shigella_peptide_pY

## Setup for summarization
shigella_global <- shigella_global %>% group_by(ProteinName, PeptideSequence, 
  Charge, PSM, Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate) %>%
  summarize(Intensity = max(Intensity))

pY_input <- list(PTM = data.frame(shigella_peptide_pY), 
                 PROTEIN = data.frame(shigella_global))
pST_input <- list(PTM = data.frame(shigella_peptide_pST), 
                  PROTEIN = data.frame(shigella_global))

## Data summarization ----------------------------------------------------------
py_summarization <- dataSummarizationPTM_TMT(pY_input)
pST_summarization <- dataSummarizationPTM_TMT(pST_input)

# save(py_summarization, file = "pY_Summarization.rda")
# save(pST_summarization, file = "pST_Summarization.rda")

py_summarized_data <- py_summarization$PTM$ProteinLevelData
py_feature_data <- py_summarization$PTM$FeatureLevelData

pst_summarized_data <- pST_summarization$PTM$ProteinLevelData
pst_feature_data <- pST_summarization$PTM$FeatureLevelData

global_summarized_data <- pST_summarization$PROTEIN$ProteinLevelData
global_feature_data <- pST_summarization$PROTEIN$FeatureLevelData


## Create pYST dataset
## Create pY/pST Dataset
## check the overlapped sites between pST and pY
shared_sites <- intersect(unique(py_summarized_data$Protein), 
                          unique(pst_summarized_data$Protein))
shared_ST <- shared_sites[grepl("_S", shared_sites)|grepl("_T", shared_sites)]
shared_Y <- shared_sites[grepl("_Y", shared_sites)]

## remove the overlapped Y sites from data_ST
data_ST <- pst_summarized_data %>% filter(!Protein %in% shared_Y)
## remove the overlapped ST sites from data_Y
data_Y <- py_summarized_data %>% filter(!Protein %in% shared_ST)

## combine pST and pY
pSTY_sum_msstats <- rbind(data_ST, data_Y)

## remove the overlapped Y sites from data_ST
data_ST_feature <- pst_feature_data %>% filter(!ProteinName %in% shared_Y)
## remove the overlapped ST sites from data_Y
data_Y_feature <- py_feature_data %>% filter(!ProteinName %in% shared_ST)

## combine pST and pY
pSTY_features <- rbind(data_ST_feature, data_Y_feature)

## Plot some interesting PTMs --------------------------------------------------
## TTP_MOUSE
colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ptm_list1 <- c("TTP_MOUSE|P22893_S178")

pSTY_features$BioReplicate <- str_replace(pSTY_features$BioReplicate, "3", "1")
pSTY_features$BioReplicate <- str_replace(pSTY_features$BioReplicate, "4", "2")

pSTY_features$BioReplicate <- factor(pSTY_features$BioReplicate,
                 levels = c("WT_Uninfect_1", "WT_Uninfect_2", "WT_Uninfect_3" , "WT_Uninfect_4",
                            "WT_Early_1", "WT_Early_2", "WT_Early_3", "WT_Early_4",
                            "WT_Late_1", "WT_Late_2", "WT_Late_3", "WT_Late_4",
                            "KO_Uninfect_1", "KO_Uninfect_2" , "KO_Uninfect_3",
                            "KO_Early_1", "KO_Early_2", "KO_Early_3", "KO_Early_4",
                            "KO_Late_1", "KO_Late_2", "KO_Late_3", "KO_Late_4"))


pSTY_sum_msstats$BioReplicate <- str_replace(pSTY_sum_msstats$BioReplicate, "3", "1")
pSTY_sum_msstats$BioReplicate <- str_replace(pSTY_sum_msstats$BioReplicate, "4", "2")

temp_plot1 <- pSTY_features %>% filter(ProteinName == ptm_list1[[1]])
temp_plot1$FeatureType <- "Peptide"

pSTY_sum_msstats$ProteinName <- pSTY_sum_msstats$Protein
pSTY_sum_msstats$log2Intensity <- pSTY_sum_msstats$Abundance
pSTY_sum_msstats$FeatureType <- "Model"

test <- rbindlist(list(temp_plot1, pSTY_sum_msstats %>% filter(
  ProteinName == ptm_list1[[1]])), fill = TRUE)

test[test$FeatureType == 'Model'][['FeatureType']] <- "PTM Summarized"
test[test$FeatureType == 'Peptide'][['FeatureType']] <- "PTM Feature"

p1 <- test %>% ggplot() +
  geom_line(aes(x = BioReplicate, y = log2Intensity, group = PSM, color = FeatureType), size = 2) +
  geom_point(aes(x = BioReplicate, y = log2Intensity, group = PSM, color = FeatureType), size = 5) +
  geom_vline(data=data.frame(x = c(2.5, 4.5, 6.5, 8.5, 10.5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed") +
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) +
  scale_size_manual(values = c(1, 2)) +
  labs(title = "TTP with phosphorylation at site S178", x = "BioReplicate", y = "Abundance") +
  facet_grid(Mixture~.) +
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
  annotate("text", x = 1.75, y = 24.75, label = "WT_Uninfect", size = 6) +
  annotate("text", x = 3.5, y = 23.25, label = "WT_Early", size = 6) +
  annotate("text", x = 5.5, y = 24.75, label = "WT_Late", size = 6) +
  annotate("text", x = 7.5, y = 23.25, label = "KO_Uninfect", size = 6) +
  annotate("text", x = 9.5, y = 24.75, label = "KO_Early", size = 6) +
  annotate("text", x = 11.5, y = 23.25, label = "KO_Late", size = 6) +
  ylim(7, 25)

p1

prot <- 'TTP_MOUSE|P22893'

global_summarized_data$ProteinName <- global_summarized_data$Protein
global_summarized_data$log2Intensity <- global_summarized_data$Abundance
global_summarized_data$FeatureType <- "Model"

prot_plot <- global_feature_data %>% filter(ProteinName == prot)
prot_plot$FeatureType <- "Peptide"

test <- rbindlist(list(prot_plot, global_summarized_data %>% filter(
  ProteinName == prot)), fill = TRUE)

test$BioReplicate <- str_replace(test$BioReplicate, "3", "1")
test$BioReplicate <- str_replace(test$BioReplicate, "4", "2")
test$BioReplicate <- factor(test$BioReplicate,
      levels = c("WT_Uninfect_1", "WT_Uninfect_2", "WT_Uninfect_3" , "WT_Uninfect_4",
                 "WT_Early_1", "WT_Early_2", "WT_Early_3", "WT_Early_4",
                 "WT_Late_1", "WT_Late_2", "WT_Late_3", "WT_Late_4",
                 "KO_Uninfect_1", "KO_Uninfect_2" , "KO_Uninfect_3",
                 "KO_Early_1", "KO_Early_2", "KO_Early_3", "KO_Early_4",
                 "KO_Late_1", "KO_Late_2", "KO_Late_3", "KO_Late_4"))

test[test$FeatureType == 'Model'][['FeatureType']] <- "Protein Summarized"
test[test$FeatureType == 'Peptide'][['FeatureType']] <- "Protein Feature"

p2 <- test %>% ggplot() +
  geom_line(aes(x = BioReplicate, y = log2Intensity, group = PSM, 
                color = FeatureType), size = 2) +
  geom_point(aes(x = BioReplicate, y = log2Intensity, group = PSM, 
                 color = FeatureType), size = 5) +
  geom_vline(data=data.frame(x = c(2.5, 4.5, 6.5, 8.5, 10.5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed") +
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) +
  scale_size_manual(values = c(1, 2)) +
  labs(title = "TTP Protein", x = "BioReplicate", y = "Abundance") +
  facet_grid(Mixture~.) +
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
  annotate("text", x = 1.75, y = 24.75, label = "WT_Uninfect", size = 6) +
  annotate("text", x = 3.5, y = 23.25, label = "WT_Early", size = 6) +
  annotate("text", x = 5.5, y = 24.75, label = "WT_Late", size = 6) +
  annotate("text", x = 7.5, y = 23.25, label = "KO_Uninfect", size = 6) +
  annotate("text", x = 9.5, y = 24.75, label = "KO_Early", size = 6) +
  annotate("text", x = 11.5, y = 23.25, label = "KO_Late", size = 6) +
  ylim(7, 25)
p2

grid.arrange(p1, p2, nrow = 1)

## KI67_MOUSE
ptm_list1 <- c("KI67_MOUSE|E9PVX6_T215")
prot <- 'KI67_MOUSE|E9PVX6'

pSTY_features$BioReplicate <- str_replace(pSTY_features$BioReplicate, "3", "1")
pSTY_features$BioReplicate <- str_replace(pSTY_features$BioReplicate, "4", "2")

pSTY_features$BioReplicate <- factor(pSTY_features$BioReplicate,
                 levels = c("WT_Uninfect_1", "WT_Uninfect_2", "WT_Uninfect_3" , "WT_Uninfect_4",
                            "WT_Early_1", "WT_Early_2", "WT_Early_3", "WT_Early_4",
                            "WT_Late_1", "WT_Late_2", "WT_Late_3", "WT_Late_4",
                            "KO_Uninfect_1", "KO_Uninfect_2" , "KO_Uninfect_3",
                            "KO_Early_1", "KO_Early_2", "KO_Early_3", "KO_Early_4",
                            "KO_Late_1", "KO_Late_2", "KO_Late_3", "KO_Late_4"))


pSTY_sum_msstats$BioReplicate <- str_replace(pSTY_sum_msstats$BioReplicate, "3", "1")
pSTY_sum_msstats$BioReplicate <- str_replace(pSTY_sum_msstats$BioReplicate, "4", "2")

temp_plot1 <- pSTY_features %>% filter(ProteinName == ptm_list1[[1]])
temp_plot1$FeatureType <- "Peptide"

pSTY_sum_msstats$ProteinName <- pSTY_sum_msstats$Protein
pSTY_sum_msstats$log2Intensity <- pSTY_sum_msstats$Abundance
pSTY_sum_msstats$FeatureType <- "Model"

test <- rbindlist(list(temp_plot1, pSTY_sum_msstats %>% filter(
  ProteinName == ptm_list1[[1]])), fill = TRUE)

test[test$FeatureType == 'Model'][['FeatureType']] <- "PTM Summarized"
test[test$FeatureType == 'Peptide'][['FeatureType']] <- "PTM Feature"

p1 <- test %>% ggplot() +
  geom_line(aes(x = BioReplicate, y = log2Intensity, group = PSM, 
                color = FeatureType), size = 2) +
  geom_point(aes(x = BioReplicate, y = log2Intensity, group = PSM, 
                 color = FeatureType), size = 5) +
  geom_vline(data=data.frame(x = c(2.5, 4.5, 6.5, 8.5, 10.5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed") +
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) +
  scale_size_manual(values = c(1, 2)) +
  labs(title = "KI67 with phosphorylation at site T215", x = "BioReplicate", 
       y = "Abundance") +
  facet_grid(Mixture~.) +
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
  annotate("text", x = 1.75, y = 24.75, label = "WT_Uninfect", size = 6) +
  annotate("text", x = 3.5, y = 23.25, label = "WT_Early", size = 6) +
  annotate("text", x = 5.5, y = 24.75, label = "WT_Late", size = 6) +
  annotate("text", x = 7.5, y = 23.25, label = "KO_Uninfect", size = 6) +
  annotate("text", x = 9.5, y = 24.75, label = "KO_Early", size = 6) +
  annotate("text", x = 11.5, y = 23.25, label = "KO_Late", size = 6) +
  ylim(7, 25)

p1

global_summarized_data$ProteinName <- global_summarized_data$Protein
global_summarized_data$log2Intensity <- global_summarized_data$Abundance
global_summarized_data$FeatureType <- "Model"

prot_plot <- global_feature_data %>% filter(ProteinName == prot)
prot_plot$FeatureType <- "Peptide"

test <- rbindlist(list(prot_plot, global_summarized_data %>% filter(
  ProteinName == prot)), fill = TRUE)

test$BioReplicate <- str_replace(test$BioReplicate, "3", "1")
test$BioReplicate <- str_replace(test$BioReplicate, "4", "2")
test$BioReplicate <- factor(test$BioReplicate,
                            levels = c("WT_Uninfect_1", "WT_Uninfect_2", "WT_Uninfect_3" , "WT_Uninfect_4",
                                       "WT_Early_1", "WT_Early_2", "WT_Early_3", "WT_Early_4",
                                       "WT_Late_1", "WT_Late_2", "WT_Late_3", "WT_Late_4",
                                       "KO_Uninfect_1", "KO_Uninfect_2" , "KO_Uninfect_3",
                                       "KO_Early_1", "KO_Early_2", "KO_Early_3", "KO_Early_4",
                                       "KO_Late_1", "KO_Late_2", "KO_Late_3", "KO_Late_4"))

test[test$FeatureType == 'Model'][['FeatureType']] <- "Protein Summarized"
test[test$FeatureType == 'Peptide'][['FeatureType']] <- "Protein Feature"

p2 <- test %>% ggplot() +
  geom_line(aes(x = BioReplicate, y = log2Intensity, group = PSM, 
                color = FeatureType), size = 2) +
  geom_point(aes(x = BioReplicate, y = log2Intensity, group = PSM, 
                 color = FeatureType), size = 5) +
  geom_vline(data=data.frame(x = c(2.5, 4.5, 6.5, 8.5, 10.5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed") +
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) +
  scale_size_manual(values = c(1, 2)) +
  labs(title = "KI67 Protein", x = "BioReplicate", y = "Abundance") +
  facet_grid(Mixture~.) +
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
  annotate("text", x = 1.75, y = 24.75, label = "WT_Uninfect", size = 6) +
  annotate("text", x = 3.5, y = 23.25, label = "WT_Early", size = 6) +
  annotate("text", x = 5.5, y = 24.75, label = "WT_Late", size = 6) +
  annotate("text", x = 7.5, y = 23.25, label = "KO_Uninfect", size = 6) +
  annotate("text", x = 9.5, y = 24.75, label = "KO_Early", size = 6) +
  annotate("text", x = 11.5, y = 23.25, label = "KO_Late", size = 6) +
  ylim(7, 25)
p2

grid.arrange(p1, p2, nrow = 1)

## Run Models ------------------------------------------------------------------
## Contrast Matrix
comparison <- read.delim("comparison_matrix.txt")
rownames(comparison) <- comparison$X
comparison <- comparison %>% dplyr::select(-X)
comparison <- as.matrix(comparison) # this step is necessary for contrast comparison

pSTY_model <- groupComparisonPTM(pyst_sum_input, data.type = "TMT", 
                                 contrast.matrix = comparison)

## Volcano Plot
model_df <- pSTY_model$PTM.Model %>% filter(Label == "WT_Late-Wt_Uninfected")
special_model_df <- model_df %>% filter(Protein == "TTP_MOUSE|P22893_S178")

v1 <- ggplot() + geom_point(data = model_df, mapping = aes(x = log2FC, y = -log10(adj.pvalue)), color = "grey") +
  geom_point(data = special_model_df, mapping = aes(x = log2FC, y = -log10(adj.pvalue)), color = "red", size = 5) +
  geom_vline(data=data.frame(x = c(-.5, .5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed", size = 1.25) +
  geom_hline(data=data.frame(x = c(-log10(.05))),
             aes(yintercept=as.numeric(x)), linetype = "dashed", size = 1.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        title = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.position = "none") +
  xlim(-5, 5) +
  ylim(0, 15.5) +
  labs(title = "WT Late vs Uninfected, by PTM quant only", 
       y = "-Log Adj. Pvalue", x = "Estimated log2 Fold Change") +
  geom_label_repel(
    data = data.frame(x = 2.899, y = -log10(.00087),
                      label = "TTP S178 \n Log2FC: 2.90 \n Adj.pvalue: .0009"),
    aes(x = x, y = y, label = label),
    label.padding = unit(0.55, "lines"),
    size = 8,
    nudge_y = 3,
    nudge_x = 0,
    size = 5,
    color = "black",
    fill="#69b3a2")
v1

adj_model_df <- pSTY_model$ADJUSTED.Model %>% filter(Label == "WT_Late-Wt_Uninfected")
adj_special_model_df <- adj_model_df %>% filter(Protein == "TTP_MOUSE|P22893_S178")

v2 <- adj_model_df %>% ggplot() +
  geom_point(mapping = aes(x = log2FC, y = -log10(adj.pvalue)), color = "grey") +
  geom_point(data = adj_special_model_df,
             mapping = aes(x = log2FC, y = -log10(adj.pvalue)),
             color = "black", size = 5) +
  geom_vline(data=data.frame(x = c(-.5, .5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed", size = 1.25) +
  geom_hline(data=data.frame(x = c(-log10(.05))),
             aes(yintercept=as.numeric(x)), linetype = "dashed", size = 1.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        title = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.position = "none") +
  xlim(-5, 5) +
  ylim(0, 15.5) +
  labs(title = "WT Late vs Uninfected, PTM adjusted by Protein", 
       y = "-Log Adj. Pvalue", x = "Estimated log2 Fold Change") +
  geom_label_repel(
    data = data.frame(x = .886, y = -log10(.248), 
                      label = "TTP S178 \n Log2FC: .886 \n Adj.pvalue: .248"),
    aes(x = x, y = y, label = label),
    label.padding = unit(0.55, "lines"),
    size = 8,
    nudge_y = 3,
    nudge_x = 2,
    size = 5,
    color = "black",
    fill="#69b3a2")

grid.arrange(v1, v2, nrow = 1)

## KI67_MOUSE
model_df <- pSTY_model$PTM.Model %>% filter(Label == "WT_Early-Wt_Uninfected")
special_model_df <- model_df %>% filter(Protein == "KI67_MOUSE|E9PVX6_T215")

v1 <- ggplot() + geom_point(data = model_df, mapping = aes(x = log2FC, y = -log10(adj.pvalue)), color = "grey") +
  geom_point(data = special_model_df, mapping = aes(x = log2FC, y = -log10(adj.pvalue)),
             color = "black", size = 5) +
  geom_vline(data=data.frame(x = c(-.5, .5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed", size = 1.25) +
  geom_hline(data=data.frame(x = c(-log10(.05))),
             aes(yintercept=as.numeric(x)), linetype = "dashed", size = 1.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        title = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.position = "none") +
  xlim(-5, 5) +
  ylim(0, 12) +
  labs(title = "WT Early vs Uninfected, by PTM quant only", 
       y = "-Log Adj. Pvalue", x = "Estimated log2 Fold Change") +
  geom_label_repel(
    data = data.frame(x = 0.185, y = -log10(0.338),
                      label = "KI67 T215 \n Log2FC: .185 \n Adj.pvalue: .338"),
    aes(x = x, y = y, label = label),
    label.padding = unit(0.55, "lines"),
    size = 8,
    nudge_y = 3,
    nudge_x = 2,
    size = 5,
    color = "black",
    fill="#69b3a2")


adj_model_df <- pSTY_model$ADJUSTED.Model %>% filter(Label == "WT_Early-Wt_Uninfected")
adj_special_model_df <- adj_model_df %>% filter(Protein == "KI67_MOUSE|E9PVX6_T215")

v2 <- adj_model_df %>% ggplot() +
  geom_point(mapping = aes(x = log2FC, y = -log10(adj.pvalue)), color = "grey") +
  geom_point(data = adj_special_model_df,
             mapping = aes(x = log2FC, y = -log10(adj.pvalue)),
             color = "red", size = 5) +
  geom_vline(data=data.frame(x = c(-.5, .5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed", size = 1.25) +
  geom_hline(data=data.frame(x = c(-log10(.05))),
             aes(yintercept=as.numeric(x)), linetype = "dashed", size = 1.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        title = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.position = "none") +
  xlim(-5, 5) +
  ylim(0, 12) +
  labs(title = "WT Early vs Uninfected, PTM adjusted by Protein", 
       y = "-Log Adj. Pvalue", x = "Estimated log2 Fold Change") +
  geom_label_repel(
    data = data.frame(x = .804, y = -log10(.01559), 
                      label = "KI67 T215 \n Log2FC: .804 \n Adj.pvalue: .0156"),
    aes(x = x, y = y, label = label),
    label.padding = unit(0.55, "lines"),
    size = 8,
    nudge_y = 3,
    nudge_x = 2,
    size = 5,
    color = "black",
    fill="#69b3a2")

grid.arrange(v1, v2, nrow = 1)

## Venndiagram -----------------------------------------------------------------

colors <- c("#E69F00", "#56B4E9")

unadj_ptm_model <- pSTY_model$PTM.Model
adj_ptm_model <- pSTY_model$ADJUSTED.Model

unadj_ptm_model$Protein_Label <- paste(unadj_ptm_model$Protein, unadj_ptm_model$Label)
adj_ptm_model$Protein_Label <- paste(adj_ptm_model$Protein, adj_ptm_model$Label)

sig_unadj_ptm_model <- unadj_ptm_model %>% filter(Protein_Label %in% adj_ptm_model$Protein_Label)

sig_unadj_ptm_model <- sig_unadj_ptm_model %>% filter(adj.pvalue < .05 & is.finite(log2FC))
sig_adj_ptm_model <- adj_ptm_model %>% filter(adj.pvalue < .05 & is.finite(log2FC))


venn.diagram(
  x = list(sig_unadj_ptm_model$Protein_Label, sig_adj_ptm_model$Protein_Label),
  category.names = c("Unadjusted" , "Adjusted"),
  filename = "shig_venn_diagramm_time_series.png",
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
  main = "Overlap between signficant adjusted and unadjusted PTMs"
)

combo <- adj_ptm_model %>% merge(unadj_ptm_model, by = "Protein_Label")

combo %>% filter(log2FC.x < log2FC.y*1.09 & log2FC.x > log2FC.y*.91 &
                   adj.pvalue.x >= .05 & adj.pvalue.y < .05)


