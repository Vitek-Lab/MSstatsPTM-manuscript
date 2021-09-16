
library(tidyverse)
library(MSstatsPTM)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("../../PTM_Refractor/Label_Free_PTM_Data/USP30/dta_usp30.RData")
head(dta_usp30)

processed_df <- readRDS("../../PTM_Refractor/Label_Free_PTM_Data/USP30/usp30-site.RDS")
head(processed_df)

processed_df$FragmentIon <- NA
processed_df$ProductCharge <- NA
processed_df$IsotopeLabelType <- 'L'
processed_df <- processed_df %>% mutate(PrecursorCharge = stringr::str_split(processed_df$feature, "_", simplify = TRUE)[,2])

colnames(processed_df) <- c('ProteinName', 'PeptideSequence', 'is_mod', 'Feature', 
                            'Condition', "Batch", "Run", 'Intensity', 'site', 'FragmentIon',
                       'ProductCharge', 'IsotopeLabelType', 'PrecursorCharge')

#processed_df <- processed_df %>% mutate(BioReplicate = stringr::str_split(processed_df$Run, "-", simplify = TRUE)[,2])
processed_df$BioReplicate <- processed_df$Batch
processed_df$Intensity <- 2^processed_df$Intensity

## Combine proteins and sites so that we can estimate PTMs (like MSstatsTMTPTM)
PTM_df <- processed_df %>% filter(is_mod == TRUE)
PTM_df$ProteinName <- paste(PTM_df$ProteinName, PTM_df$site, sep = '_')
protein_df <- processed_df %>% filter(is_mod == FALSE)

# mean((PTM_df %>% group_by(ProteinName) %>% summarize(num_feat = n_distinct(Feature)))$num_feat)
# mean((protein_df %>% group_by(ProteinName) %>% summarize(num_feat = n_distinct(Feature)))$num_feat)

drops <- c("is_mod", "site", "Batch", "Feature")
PTM_df <- PTM_df[ , !(names(PTM_df) %in% drops)]
protein_df <- protein_df[ , !(names(protein_df) %in% drops)]

processed_df <- list('PTM' = PTM_df, 'PROTEIN' = protein_df)

## Fit MSstatsPTM Model --------------------------------------------------------
summarized_ptm <- dataSummarizationPTM(processed_df, min_feature_count.PTM = 3)

model_ptm <- groupComparisonPTM(summarized_ptm, data.type = "LabelFree")

## ----------
## Test some proteins that might have been missed in analysis
adj_model <- model_ptm$ADJUSTED.Model
ptm_model <- model_ptm$PTM.Model

adj_model %>% merge(ptm_model, by = c("Protein", "Label")) %>% filter(adj.pvalue.x < .05 & adj.pvalue.y >= .05 & Label == "CCCP vs USP30_OE")

## Fit Limma -------------------------------------------------------------------
library(limma)

## Fit limma model (both not and with adjust) for given dataset
fit_limma <- function(summarized_data, conditions, runs){
  
  ## Convert into format required for limma
  input <- data.frame(summarized_data %>% select(PTM, Run, Condition, Abundance.x) %>% 
                        pivot_wider(names_from = c(Condition, Run), values_from = Abundance.x,
                                    names_sort = TRUE))
  rownames(input) <- input$PTM
  input <- input %>% select(-PTM)
  input_adj <- data.frame(summarized_data %>% select(PTM, Run, Condition, Adj_Abundance) %>% 
                            pivot_wider(names_from = c(Condition, Run), values_from = Adj_Abundance,
                                        names_sort = TRUE))
  rownames(input_adj) <- input_adj$PTM
  input_adj <- input_adj %>% select(-PTM)
  
  ## Create contrast matrix
  class <- c()
  for (x in seq_len(conditions)){
    cond <- rep(paste0("mix_", as.character(x)), runs)
    class <- c(class, cond)
  }
  class <- as.factor(class)
  design <- model.matrix(~0+class)
  
  input.matrix <- as.matrix(input)
  input_adj.matrix <- as.matrix(input_adj)
  
  ## Run models
  fit <- lmFit(input.matrix, design=design)
  fit_adj <- lmFit(input_adj.matrix, design=design)
  
  contrast.matrix <- design.pairs(colnames(design))
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  fit2_adj <- contrasts.fit(fit_adj, contrast.matrix)
  fit2_adj <- eBayes(fit2_adj)
  
  ## Output to data.table
  comparisons <- data.table()
  comparisons_adj <- data.table()
  
  for (g in seq_along(colnames(fit2$coefficients))){
    comparisons <- rbindlist(list(comparisons, 
                                  data.table(PTM = rownames(fit2$coefficients),
                                             Label = colnames(fit2$coefficients)[g],
                                             Log2FC = as.vector(fit2$coefficients[,g]),
                                             pvalue = as.vector(fit2$p.value[,g]),
                                             df = as.vector(fit2$df.residual),
                                             se = as.vector(fit2$sigma))))
    comparisons_adj <- rbindlist(list(comparisons_adj, 
                                      data.table(PTM = rownames(fit2_adj$coefficients),
                                                 Label = colnames(fit2_adj$coefficients)[g],
                                                 Log2FC = as.vector(fit2_adj$coefficients[,g]),
                                                 pvalue = as.vector(fit2_adj$p.value[,g]),
                                                 df = as.vector(fit2_adj$df.residual),
                                                 se = as.vector(fit2_adj$sigma))))
  }
  
  return(list(limma_test = comparisons, limma_adj_test = comparisons_adj))
}
## Limma pairwise function
design.pairs <- function(levels) {
  n <- length(levels)
  design <- matrix(0,n,choose(n,2))
  rownames(design) <- levels
  colnames(design) <- 1:choose(n,2)
  k <- 0
  for (i in 1:(n-1))
    for (j in (i+1):n) {
      k <- k+1
      design[i,k] <- 1
      design[j,k] <- -1
      colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
    }
  design
}

PTM_df_limma <- PTM_df
PTM_df_limma$Intensity <- log2(PTM_df_limma$Intensity)
protein_df_limma <- protein_df
protein_df_limma$Intensity <- log2(protein_df_limma$Intensity)

## Normalize - equalize medians
## KGG
kgg_med <- PTM_df_limma %>% group_by(Run) %>% 
  summarize(med = median(Intensity, na.rm=TRUE)) %>% 
  mutate(inty_expected = median(med), adjIntensity = inty_expected - med) %>% 
  ungroup() %>% select(Run, adjIntensity)

PTM_df_limma <- PTM_df_limma %>% left_join(kgg_med, by = 'Run')
PTM_df_limma$Intensity <- PTM_df_limma$Intensity + PTM_df_limma$adjIntensity
PTM_df_limma$Intensity <- 2^PTM_df_limma$Intensity

## Prot
prot_med <- protein_df_limma %>% group_by(Run) %>% 
  summarize(med = median(Intensity, na.rm=TRUE)) %>% 
  mutate(inty_expected = median(med), adjIntensity = inty_expected - med) %>% 
  ungroup() %>% select(Run, adjIntensity)
protein_df_limma <- protein_df_limma %>% left_join(prot_med, by = 'Run')
protein_df_limma$Intensity <- protein_df_limma$Intensity + protein_df_limma$adjIntensity
protein_df_limma$Intensity <- 2^protein_df_limma$Intensity

## Summarize
summarized_ptm <- data.table()
summarized_proteins <- data.table()
runs <- unique(PTM_df_limma$Run)
for (r in seq_along(runs)){
  sum_runs <- PTM_df_limma %>% filter(Run == runs[[r]]) %>%
    group_by(ProteinName, Condition, Run) %>%
    summarize(Abundance = log2(sum(Intensity, na.rm = TRUE)))
  summarized_ptm <- rbindlist(list(summarized_ptm, sum_runs))
  
  sum_runs_prot <- protein_df_limma %>% filter(Run == runs[[r]]) %>%
    group_by(ProteinName, Condition, Run) %>%
    summarize(Abundance = log2(sum(Intensity, na.rm = TRUE)))
  summarized_proteins <- rbindlist(list(summarized_proteins, sum_runs_prot))
}

summarized_ptm$PTM <- summarized_ptm$ProteinName
summarized_ptm$ProteinName <- sapply(summarized_ptm$PTM, function(x) {str_split(x, "_",2)[[1]][1]})
joined <- merge(summarized_ptm, summarized_proteins, by = c("ProteinName", "Run", "Condition"), all.x = TRUE)
joined <- joined %>% filter(is.finite(Abundance.y))
joined$Adj_Abundance <- joined$Abundance.x - joined$Abundance.y

limma_test_res <- fit_limma(joined, 4, 4)
limma_test_res$limma_test$adj.pvalue <- 0.
limma_test_res$limma_adj_test$adj.pvalue <- 0.

## adjust pvalue
for (l in seq_along(unique(limma_test_res$limma_test$Label))){
  limma_test_res$limma_test[Label == unique(limma_test_res$limma_test$Label)[l]]$adj.pvalue <- p.adjust(
    limma_test_res$limma_test[Label == unique(limma_test_res$limma_test$Label)[l]]$pvalue, method = "BH")
  limma_test_res$limma_adj_test[Label == unique(limma_test_res$limma_adj_test$Label)[l]]$adj.pvalue <- p.adjust(
    limma_test_res$limma_adj_test[Label == unique(limma_test_res$limma_adj_test$Label)[l]]$pvalue, method = "BH")
}

limma_test_res$limma_adj_test <- limma_test_res$limma_adj_test %>% mutate(Label = ifelse(Label == "classmix_1-classmix_2", "CCCP vs Combo", 
                                                                         ifelse(Label == "classmix_1-classmix_3", "CCCP vs Ctrl",
                                                                                ifelse(Label == "classmix_1-classmix_4", "CCCP vs USP30_OE",
                                                                                       ifelse(Label == "classmix_2-classmix_3", "Combo vs Ctrl",
                                                                                              ifelse(Label == "classmix_2-classmix_4", "Combo vs USP30_OE",
                                                                                                     "Ctrl vs USP30_OE"))))))

limma_test_res$limma_test <- limma_test_res$limma_test %>% mutate(Label = ifelse(Label == "classmix_1-classmix_2", "CCCP vs Combo", 
                                                                         ifelse(Label == "classmix_1-classmix_3", "CCCP vs Ctrl",
                                                                                ifelse(Label == "classmix_1-classmix_4", "CCCP vs USP30_OE",
                                                                                       ifelse(Label == "classmix_2-classmix_3", "Combo vs Ctrl",
                                                                                              ifelse(Label == "classmix_2-classmix_4", "Combo vs USP30_OE",
                                                                                                     "Ctrl vs USP30_OE"))))))

## Fit T-test ------------------------------------------------------------------
run_ttest <- function(summarized_data){
  
  ptms <- unique(summarized_data$PTM)
  
  ttest_temp <- data.table()
  ttest_adj_temp <- data.table()
  ## Run ttest for each ptm
  for (p in seq_along(ptms)){
    temp_joined <- summarized_data %>% filter(PTM == ptms[[p]])
    groups <- (temp_joined %>% distinct(Condition))[[1]]
    t <- 2
    
    ## Loop over groups
    for (g in 1:(length(groups)-1)){
      for (g2 in (g+1):length(groups)){
        tryCatch({ttest_ptm <- t.test((temp_joined %>% filter(Condition == groups[g]) %>% select(Abundance.x))[[1]], 
                                      (temp_joined %>% filter(Condition == groups[g2]) %>% select(Abundance.x))[[1]])
        ttest_temp <- rbindlist(list(ttest_temp, data.table(ptm = ptms[[p]], 
                                                            label = paste(groups[g], "vs", groups[g2], sep = " "), 
                                                            pval = ttest_ptm$p.value, 
                                                            tstat = ttest_ptm$statistic[[1]],
                                                            SE = ttest_ptm$stderr,
                                                            df = ttest_ptm$parameter[[1]],
                                                            estimate = ttest_ptm$estimate[[2]] - ttest_ptm$estimate[[1]])))
        },
        error=function(e){cat("ERROR :", as.character(p), "\n")}
        )
        tryCatch({ttest_adj_ptm <- t.test((temp_joined %>% filter(Condition == groups[g]) %>% select(Adj_Abundance))[[1]], 
                                          (temp_joined %>% filter(Condition == groups[g2]) %>% select(Adj_Abundance))[[1]])
        ttest_adj_temp <- rbindlist(list(ttest_adj_temp, data.table(ptm = ptms[[p]], 
                                                                    label = paste(groups[g], "vs", groups[g2], sep = " "), 
                                                                    pval = ttest_adj_ptm$p.value, 
                                                                    tstat = ttest_adj_ptm$statistic[[1]],
                                                                    SE = ttest_adj_ptm$stderr,
                                                                    df = ttest_adj_ptm$parameter[[1]],
                                                                    estimate = ttest_adj_ptm$estimate[[2]] - ttest_adj_ptm$estimate[[1]])))},
        error=function(e){cat("ERROR :", as.character(p), "\n")})
      }
    }
  }
  return(list(ttest_temp = ttest_temp, ttest_adj_temp = ttest_adj_temp))
}

ttest_model <- run_ttest(joined)

ttest_model$ttest_adj_temp$adj.pvalue <- 0.
ttest_model$ttest_temp$adj.pvalue <- 0.

for (l in seq_along(unique(ttest_model$ttest_adj_temp$label))){
  ttest_model$ttest_adj_temp[label == unique(ttest_model$ttest_adj_temp$label)[l]]$adj.pvalue <- p.adjust(
    ttest_model$ttest_adj_temp[label == unique(ttest_model$ttest_adj_temp$label)[l]]$pval, method = "BH")
  ttest_model$ttest_temp[label == unique(ttest_model$ttest_temp$label)[l]]$adj.pvalue <- p.adjust(
    ttest_model$ttest_temp[label == unique(ttest_model$ttest_temp$label)[l]]$pval, method = "BH")
}

## Analyze model differences ---------------------------------------------------
library(VennDiagram)

msstats_sig <- model_ptm$ADJUSTED.Model %>% filter(adj.pvalue < .05) %>% mutate(MSstats_labels = paste0(Protein, "_", Label))
limma_sig <- limma_test_res$limma_adj_test %>% filter(adj.pvalue < .05) %>% mutate(Limma_labels = paste0(PTM, "_", Label))
ttest_sig <- ttest_model$ttest_adj_temp %>% filter(adj.pvalue < .05) %>% mutate(ttest_labels = paste0(ptm, "_", label))

venn.diagram(
  x = list(msstats_sig$MSstats_labels, limma_sig$Limma_labels, ttest_sig$ttest_labels),
  category.names = c("MSstatsPTM" , "Limma" , "T-test"),
  filename = "usp30_ven.png",
  output=TRUE
)

msstats_sig %>% select(MSstats_labels) %>% 
  full_join(limma_sig %>% select(Limma_labels), 
            by = c("MSstats_labels" = "Limma_labels"), keep = TRUE) %>% 
  filter(is.na(MSstats_labels))




