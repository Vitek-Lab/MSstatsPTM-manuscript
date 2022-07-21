
library(MSstatsPTM)
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(VennDiagram)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ptm_input <- read.csv("../data/IpaH_PTM_Input.txt")
protein_input <- read.csv("../data/Ipah_Global_Input.txt", sep = "\t")
protein_input$FeatureType <- "Peptide"
protein_input$Abundance <- log2(protein_input$Intensity)

ptm_summarization <- read.csv("../data/Ipah_KGG_TMT11_corrected_msstats_quant_combined_results.txt",
                              sep = "\t")
summarized_model <- ptm_summarization %>% filter(PeptideSequence == "Model")
ptm_summarization <- ptm_summarization %>% filter(PeptideSequence != "Model")
protein_summarization <- read.csv("../data/Ipah_Global_TMT11_msstats_quant_results.txt",
                                  sep = "\t")

protein_summarization = read.csv("../data/Ipah_Global_Input.txt",
         sep = "\t")

ptm_summarization %>% head()
ptm_summarization$global_prot = sapply(ptm_summarization$ProteinName, function(x){str_split(x, "\\|")[[1]][1]})

ptm_summarization$global_prot %>% n_distinct()
ptm_input$ProteinName %>% n_distinct()
ptm_input$PSM %>% n_distinct()

protein_summarization$ProteinName %>% n_distinct()
protein_summarization$PeptideSequence %>% n_distinct()


protein_summarization$PeptideSequence <- "Model"
protein_summarization$FeatureType <- "Model"
protein_summarization_combined <- rbindlist(list(protein_input, protein_summarization), fill = TRUE)


unadj_ptm_model <- read.csv("../data/Ipah_PTM_test_results.txt", sep = "\t")
protein_model <- read.csv("../data/Ipah_Global_TMT11_msstats_test_results.txt", sep = "\t")
adj_ptm_model <- read.csv("../data/Ipah_PTM_Adjusted_test_results.txt", sep = "\t")

combined_models <- merge(adj_ptm_model, unadj_ptm_model, all.x=TRUE, all.y=TRUE, by = c("Protein", "Label"))

combined_models$FC_Diff <- abs(combined_models$log2FC.x - combined_models$log2FC.y)
combined_models %>% arrange(desc(FC_Diff))

## SE analysis
unadj_ptm_model %>% ggplot + geom_boxplot(aes(y = SE))
protein_model %>% ggplot + geom_boxplot(aes(y = SE))

## Profile Plots
"GSDMD_HUMAN|P57764_K204"
"GSDMD_HUMAN|P57764_K62"

gsmd_sum <- protein_summarization_combined %>% filter(ProteinName == "GSDMD_HUMAN|P57764")
gsmd_204_sum <- ptm_summarization %>% filter(ProteinName == "GSDMD_HUMAN|P57764_K204")
gsmd_62_sum <- ptm_summarization %>% filter(ProteinName == "GSDMD_HUMAN|P57764_K62")

gsmd_62_sum$BioReplicate <- factor(gsmd_62_sum$BioReplicate,
                                     levels = c("NoDox0hr_1", "NoDox0hr_2", "NoDox6hr_1",
                                                "NoDox6hr_2", "Dox1hr_1", "Dox2hr_1",
                                                "Dox2hr_2", "Dox4hr_1",
                                                "Dox4hr_2", "Dox6hr_1", "Dox6hr_2"))

gsmd_sum$BioReplicate <- factor(gsmd_sum$BioReplicate,
                                levels = c("NoDox0hr_1", "NoDox0hr_2", "NoDox6hr_1",
                                           "NoDox6hr_2", "Dox1hr_1", "Dox2hr_1",
                                           "Dox2hr_2", "Dox4hr_1",
                                           "Dox4hr_2", "Dox6hr_1", "Dox6hr_2"))
gsmd_62_sum <- as.data.table(gsmd_62_sum)
gsmd_62_sum[gsmd_62_sum$FeatureType == 'Model'][['FeatureType']] <- "PTM Summarized"
gsmd_62_sum[gsmd_62_sum$FeatureType == 'Peptide'][['FeatureType']] <- "PTM Feature"

gsmd_sum[gsmd_sum$FeatureType == 'Model'][['FeatureType']] <- "Protein Summarized"
gsmd_sum[gsmd_sum$FeatureType == 'Peptide'][['FeatureType']] <- "Protein Feature"

p1 <- gsmd_62_sum %>% ggplot() + geom_line(aes(x = BioReplicate, y = Abundance,
                                               color = FeatureType,
                                               group = PSM, size = FeatureType)) +
  geom_point(aes(x = BioReplicate, y = Abundance, color = FeatureType), size = 5) +
  geom_vline(data=data.frame(x = c(2.5, 4.5, 5.5, 7.5, 9.5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed") +
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) +
  scale_size_manual(values = c(1.5,1)) +
  labs(title = "GSDMD with ubiquitination at site K62", x = "BioReplicate",
       y = "Abundance") +
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
  ylim(11.5, 23) +
  annotate("text", x = 1.68, y = 22.75, label = "No_Dox0hr", size = 8) +
  annotate("text", x = 3.5, y = 22., label = "No_Dox6hr", size = 8) +
  annotate("text", x = 5, y = 22.75, label = "Dox1hr", size = 8) +
  annotate("text", x = 6.5, y = 22., label = "Dox2hr", size = 8) +
  annotate("text", x = 8.5, y = 22.75, label = "Dox4hr", size = 8) +
  annotate("text", x = 10.5, y = 22., label = "Dox6hr", size = 8)

p2 <- gsmd_sum %>% ggplot() + geom_line(aes(x = BioReplicate, y = Abundance,
                                            color = FeatureType,
                                         group = PSM, size = FeatureType)) +
  geom_point(aes(x = BioReplicate, y = Abundance, color = FeatureType), size = 5) +
  geom_vline(data=data.frame(x = c(2.5, 4.5, 5.5, 7.5, 9.5)),
             aes(xintercept=as.numeric(x)), linetype = "dashed") +
  scale_colour_manual(values = c("#C3C3C3", "#D55E00")) +
  scale_size_manual(values = c(1.5,1)) +
  labs(title = "GSDMD Protein", x = "BioReplicate", y = "Abundance") +
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
  ylim(11.5, 23) +
  annotate("text", x = 1.68, y = 22.75, label = "No_Dox0hr", size = 8) +
  annotate("text", x = 3.5, y = 22., label = "No_Dox6hr", size = 8) +
  annotate("text", x = 5, y = 22.75, label = "Dox1hr", size = 8) +
  annotate("text", x = 6.5, y = 22., label = "Dox2hr", size = 8) +
  annotate("text", x = 8.5, y = 22.75, label = "Dox4hr", size = 8) +
  annotate("text", x = 10.5, y = 22., label = "Dox6hr", size = 8)

grid.arrange(p1, p2, nrow = 1)

## Volcano Plot


unadj_plot_df <- unadj_ptm_model %>% filter(Label == "Dox1hr_v_Dox4hr")
unadj_special_plot_df <- unadj_plot_df %>% filter(Protein == "GSDMD_HUMAN|P57764_K62")

v1 <- unadj_plot_df %>% ggplot() +
  geom_point(mapping = aes(x = log2FC, y = -log10(adj.pvalue)), color = "grey") +
  geom_point(data = unadj_special_plot_df,
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
        title = element_text(size = 22),
        strip.text = element_text(size = 16),
        legend.position = "none") +
  xlim(-5, 5) +
  ylim(0, 8) +
  labs(title = "Dox 4 hr vs Dox 1 hr, by PTM quant only",
       y = "-Log Adj. Pvalue",
       x = "Estimated log2 Fold Change") +
  geom_label_repel(
    data = data.frame(x = -.500, y = -log10(.064),
                      label = "GSDMD K62 \n Log2FC: -.501 \n Adj.pvalue: .0644"),
    aes(x = x, y = y, label = label),
    label.padding = unit(0.55, "lines"),
    size = 8,
    nudge_y = 3,
    nudge_x = 3,
    size = 5,
    color = "black",
    fill="#69b3a2")

adj_plot_df <- adj_ptm_model %>% filter(Label == "Dox1hr_v_Dox4hr")
adj_special_plot_df <- adj_plot_df %>% filter(Protein == "GSDMD_HUMAN|P57764_K62")

v2 <- adj_plot_df %>% ggplot() +
  geom_point(mapping = aes(x = log2FC, y = -log10(adj.pvalue)), color = "grey") +
  geom_point(data = adj_special_plot_df,
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
        title = element_text(size = 22),
        strip.text = element_text(size = 16),
        legend.position = "none") +
  xlim(-5, 5) +
  ylim(0, 8) +
  labs(title = "Dox 4 hr vs Dox 1 hr, PTM adjusted by Protein",
       y = "-Log Adj. Pvalue",
       x = "Estimated log2 Fold Change") +
  geom_label_repel(
    data = data.frame(x = 2.789747, y = -log10(5.249283e-08),
                      label = "GSDMD K62 \n Log2FC: 2.79 \n Adj.pvalue: 5.25e-08"),
    aes(x = x, y = y, label = label),
    label.padding = unit(0.55, "lines"),
    size = 8,
    nudge_y = 0,
    nudge_x = -3,
    size = 5,
    color = "black",
    fill="#69b3a2")

grid.arrange(v1, v2, nrow = 1)


## Venn Diagram
colors <- c("#E69F00", "#56B4E9")

unadj_ptm_model <- read.csv("../data/Ipah_PTM_test_results.txt", sep = "\t")
adj_ptm_model <- read.csv("../data/Ipah_PTM_Adjusted_test_results.txt", sep = "\t")

unadj_ptm_model$Protein_Label <- paste(unadj_ptm_model$Protein, unadj_ptm_model$Label)
adj_ptm_model$Protein_Label <- paste(adj_ptm_model$Protein, adj_ptm_model$Label)

sig_unadj_ptm_model <- unadj_ptm_model %>% filter(Protein_Label %in% adj_ptm_model$Protein_Label)

sig_unadj_ptm_model <- sig_unadj_ptm_model %>% filter(adj.pvalue < .05)
sig_adj_ptm_model <- adj_ptm_model %>% filter(adj.pvalue < .05)


venn.diagram(
  x = list(sig_unadj_ptm_model$Protein_Label, sig_adj_ptm_model$Protein_Label),
  category.names = c("Unadjusted" , "Adjusted"),
  filename = "ipah_venn_diagramm.png",
  output=TRUE,
  imagetype="png" ,
  height = 1500,
  width = 1500,
  resolution = 100,
  lwd = 2,
  fill = colors,
  main.fontface = "bold",
  main.fontfamily="sans",
  main.pos = c(.5,.965),
  fontface = "bold",
  cat.fontface = "bold",
  cex = 3,
  cat.cex = 3,
  cat.pos = c(-40, 30),
  cat.dist = c(.037, .03),
  main.cex = 2.8,
  main = "Dataset 4: Differentially abundant adjusted and unadjusted PTMs"
)

combined_models <- merge(adj_ptm_model, unadj_ptm_model, all.x=TRUE, all.y=TRUE, by = c("Protein_Label"))
combined_models %>% filter(log2FC.x < log2FC.y*1.1 & log2FC.x > log2FC.y*.9 &
                   adj.pvalue.x >= .05 & adj.pvalue.y < .05)
