
library(data.table)
library(MSstatsPTM)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Function from MSstats to add in missing intensity peaks
.makeBalancedDesign = function(input, fill_missing) {
  feature = NULL

  is_tmt = is.element("Channel", colnames(input))
  input$feature = paste(input$PSM, input$FragmentIon, input$ProductCharge, sep = '_')
  input$Fraction = 1

  if (fill_missing) {
    cols = intersect(colnames(input),
                     c("ProteinName", "feature", "PeptideSequence", "PSM",
                       "PrecursorCharge", "FragmentIon", "ProductCharge"))
    annotation_cols = intersect(colnames(input),
                                c("Run", "Condition", "BioReplicate", "Channel",
                                  "Mixture", "TechRepMixture", "TechReplicate"))
    intensity_ids = intersect(c("feature", "Run", "Channel", "IsotopeLabelType",
                                "Fraction"), colnames(input))
    if (is_tmt) {
      group_col = "Run"
      measurement_col = "Channel"
      cols = intersect(c(cols, "Fraction"),
                       colnames(input))
    } else {
      group_col = "Fraction"
      measurement_col = "Run"
    }
    all_possibilities = .getFullDesign(input, group_col, "feature",
                                       measurement_col, is_tmt)
    all_possibilities = merge(
      all_possibilities,
      unique(input[, c(cols, group_col), with = FALSE]),
      all.x = TRUE, by = c("feature", group_col), allow.cartesian = TRUE)
    all_possibilities = merge(
      all_possibilities,
      unique(input[, annotation_cols, with = FALSE]),
      all.x = TRUE, by = unique(c("Run", measurement_col)))
    intensities = intersect(c(intensity_ids, "Intensity", "isZero"),
                            colnames(input))
    input = merge(all_possibilities,
                  unique(input[, intensities, with = FALSE]),
                  all.x = TRUE, by = intensity_ids)
  } else {
    if (!is_tmt) {
      any_missing = as.character(unique(.getMissingRunsPerFeature(input)[, feature]))
    }
  }
  input
}
.getFullDesign = function(input, group_col, feature_col,
                          measurement_col, is_tmt
) {
  if (is_tmt) {
    labels = "L"
    groups = unique(input[[group_col]])
    by_group = vector("list", length(groups))
    measurements = unique(input[[measurement_col]])
    for (group_id in seq_along(groups)) {
      group = groups[group_id]
      group_filter = input[[group_col]] == group
      by_group[[group_id]] = data.table::as.data.table(
        expand.grid(labels = labels,
                    features = unique(input[[feature_col]][group_filter]),
                    measurements = measurements))
      by_group[[group_id]]$group = group
    }
    result = data.table::rbindlist(by_group)
    colnames(result) = c("IsotopeLabelType", feature_col,
                         measurement_col, group_col)
    result[, 2:4, with = FALSE]
  } else {
    labels = unique(input[["IsotopeLabelType"]])
    groups = unique(input[[group_col]])
    by_group = vector("list", length(groups))
    measurements = unique(input[[measurement_col]])
    for (group_id in seq_along(groups)) {
      group = groups[group_id]
      group_filter = input[[group_col]] == group
      by_group[[group_id]] = data.table::as.data.table(
        expand.grid(
          labels = labels,
          features = unique(input[[feature_col]][group_filter]),
          measurements = unique(input[[measurement_col]][group_filter])
        ))
      by_group[[group_id]]$group = group
    }
    result = data.table::rbindlist(by_group)
    colnames(result) = c("IsotopeLabelType", feature_col,
                         measurement_col, group_col)
    result
  }
}


# Prepare data set -------------------------------------------------------------

# Design for global profiling data
design_prot1 <- read.table("000547_ptm_benchmark/ptm_tft_global_ip1.txt",
  header = TRUE, sep = "\t")
design_prot1 <- as.data.table(design_prot1)
design_prot2 <- read.table("000547_ptm_benchmark/ptm_tft_global_ip2.txt",
  header = TRUE, sep = "\t")
design_prot2 <- as.data.table(design_prot2)
design_prot2[,id_injectionset := 2]
design_prot <- rbindlist(list(design_prot1, design_prot2))

# Design for KGG enriched data
design_kgg1 <- read.table("000547_ptm_benchmark/ptm_tft_kgg_ip1.txt",
  header = TRUE, sep = "\t")
design_kgg1 <- as.data.table(design_kgg1)
design_kgg1[, batch := "BCH1"]
design_kgg2 <- read.table("000547_ptm_benchmark/ptm_tft_kgg_ip2.txt",
  header = TRUE, sep = "\t")
design_kgg2 <- as.data.table(design_kgg2)
design_kgg2[, batch := "BCH2"]

design_kgg <- rbindlist(list(design_kgg1, design_kgg2))

# Prepare data set with columns of interest for later use
fn_kgg <- paste0("vistaga_", design_kgg$run_id,
                 "_kgg_with_unmodified_peptides.tsv")
fn_prot <- paste0("vistaga_", design_prot$run_id, ".tsv")

sub_cols <- c(
  "Reference", "run_id", "peptide_id", "vista_peak_area_light",
  "vista_confidence_score", "peptide_m_z", "peptide_sequence",
  "peptide_trypticity", "peptide_miscleavages", "peptide_validity",
  "label1", "ms2_rt", "ms2_charge", "is_xq"
)

dta_kgg <- vector("list", length(fn_kgg))
for (i in seq_along(fn_kgg)) {
  path_file <- paste0("000547_ptm_benchmark/", fn_kgg[i])
  dta_kgg[[i]] <- read.table(
    path_file, header = TRUE, sep = "\t", fill = TRUE
  ) %>%
    as_tibble() %>%
    select(all_of(sub_cols))
}
dta_kgg <- bind_rows(dta_kgg)

dta_prot <- vector("list", length(fn_prot))
for (i in seq_along(fn_prot)) {
  path_file <- paste0("000547_ptm_benchmark/", fn_prot[i])
  dta_prot[[i]] <- read.table(
    path_file, header = TRUE, sep = "\t", fill = TRUE
  ) %>%
    as_tibble() %>%
    select(all_of(sub_cols))
}
dta_prot <- bind_rows(dta_prot)

# AQUA peptides enriched by immunoprecipitation (IP)
aqua1_1 <- read_csv("000547_ptm_benchmark/KGG_IP1_SpikeIn_First25_MSPlorer.csv")
aqua1_2 <- read_csv("000547_ptm_benchmark/KGG_IP1_SpikeIn_Second25_MSPlorer.csv")
aqua2_1 <- read_csv("000547_ptm_benchmark/KGG_IP2_SpikeIn_First25_MSPlorer.csv")
aqua2_2 <- read_csv("000547_ptm_benchmark/KGG_IP2_SpikeIn_Second25_MSPlorer.csv")

aqua1 <- bind_rows(aqua1_1, aqua1_2) %>%
  select(PeptideSequence = Label, run_id = `Run ID`,
         ann_label = `Run Label`, Intensity = Area)
aqua2 <- bind_rows(aqua2_1, aqua2_2) %>%
  select(PeptideSequence = Label, run_id = `Run ID`,
         ann_label = `Run Label`, Intensity = Area)
aqua <- bind_rows(
  aqua1 %>% mutate(batch = "BCH1"),
  aqua2 %>% mutate(batch = "BCH2")
)

# Load data set & housekeeping --------------------------------------------

# Parameters
min_len_peptide <- 7  # Minimum acceptable length of peptide

# Modifications of interest
mod_residue <- "K"
mod_symbol <- "\\*"

regex_uniprot_iso <- regex(paste0(
  "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})",
  "([-]\\d{1,}){0,1}"
))

# Load prepared dataset
# load("work/dta_spike_vista.RData")

# Design matrix
design_prot <- design_prot %>%
  rename(group = cond_treatment) %>%
  mutate(TechReplicate = str_c("T", id_injectionset),
         Run = str_c(group, TechReplicate, sep = "-"),
         BioReplicate = id_subject)

design_kgg <- design_kgg %>%
  rename(group = cond_treatment) %>%
  mutate(TechReplicate = str_c("T", id_injectionset),
         Run = str_c(group, TechReplicate, sep = "-"),
         BioReplicate = id_subject)

df_prot <- dta_prot %>%
  filter(vista_confidence_score >= 83, !is.na(ms2_charge), !is.na(ms2_rt)) %>%
  select(
    -vista_confidence_score, -peptide_trypticity,
    -peptide_miscleavages, -peptide_validity
  ) %>%
  mutate(uniprot_iso = str_extract(Reference, pattern = regex_uniprot_iso)) %>%
  rename(peptide_ext = peptide_sequence) %>%
  mutate(PeptideSequence = str_extract(
    peptide_ext, "(?<=\\.)[ACDEFGHIKLMNPQRSTVWY\\*]+(?=\\.)")) %>%
  mutate(
    is_mod = str_detect(PeptideSequence, mod_symbol),
    PrecursorCharge = ms2_charge,
    PSM = str_c(PeptideSequence, PrecursorCharge, sep = "_"),
    Intensity = ifelse(vista_peak_area_light <= 1, 0, vista_peak_area_light)
  ) %>%
  left_join(design_prot %>% select(Run, group, run_id, BioReplicate))

df_kgg <- dta_kgg %>%
  filter(str_detect(Reference, "HUMAN")) %>%
  filter(vista_confidence_score >= 83, !is.na(ms2_charge), !is.na(ms2_rt)) %>%
  select(
    -vista_confidence_score, -peptide_trypticity,
    -peptide_miscleavages, -peptide_validity
  ) %>%
  mutate(uniprot_iso = str_extract(Reference, pattern = regex_uniprot_iso)) %>%
  rename(peptide_ext = peptide_sequence) %>%
  mutate(PeptideSequence = str_extract(
    peptide_ext, "(?<=\\.)[ACDEFGHIKLMNPQRSTVWY\\*]+(?=\\.)")) %>%
  mutate(
    is_mod = str_detect(PeptideSequence, mod_symbol),
    PrecursorCharge = ms2_charge,
    PSM = str_c(PeptideSequence, PrecursorCharge, sep = "_"),
    Intensity = ifelse(vista_peak_area_light <= 1, 0, vista_peak_area_light)
  ) %>%
  left_join(design_kgg %>% select(Run, batch, group, run_id, BioReplicate))

aqua <- aqua %>%
  mutate(ann_label = str_replace(ann_label, "_SpikeIn", "")) %>%
  left_join(design_kgg %>% select(ann_label, group, BioReplicate,
                                  TechReplicate, run_id),
            by = c('ann_label')) %>% rename(Run = run_id.y, run_id = run_id.x)

# Use the max for duplicate features (per run)
df_prot <- df_prot %>%
  filter(!is.na(uniprot_iso)) %>%
  group_by(uniprot_iso, PeptideSequence, PrecursorCharge,
           is_mod, PSM, group, Run, BioReplicate) %>%
  summarise(Intensity = max(Intensity), .groups = "drop")

df_kgg <- df_kgg %>%
  filter(!is.na(uniprot_iso)) %>%
  group_by(uniprot_iso, PeptideSequence, PrecursorCharge,
           is_mod, PSM, batch, group, Run, BioReplicate) %>%
  summarise(Intensity = max(Intensity), .groups = "drop")

aqua <- aqua %>%
  group_by(PeptideSequence, run_id, ann_label, batch, group, TechReplicate,
           BioReplicate, Run) %>%
  summarise(Intensity = max(Intensity), .groups = "drop")

# Use only second batch of the KGG (one run missing in the first batch)
design_kgg <- design_kgg %>% filter(batch == "BCH2")

df_kgg <- df_kgg %>% filter(batch == "BCH2") %>% select(-batch)

aqua <- aqua %>% filter(batch == "BCH2") %>% select(-batch)

# Site representation for modified peptides -------------------------------

df_fasta <- bind_rows(
  tidyFasta("../Sequence/homo_sapiens_all_20160725.fasta"),
  tidyFasta("../Sequence/Swissprot-Ecoli-K12.fasta")
)

peptide_seq <- df_kgg %>%
  distinct(uniprot_iso, PeptideSequence) %>%
  mutate(peptide_unmod = str_remove_all(PeptideSequence, mod_symbol)) %>%
  filter(str_length(peptide_unmod) >= min_len_peptide)

peptide_site <- locatePTM(
  peptide_seq$PeptideSequence, peptide_seq$uniprot_iso, df_fasta,
  mod_residue, mod_symbol, rmConfound = TRUE
)

df_site <- df_kgg %>%
  inner_join(peptide_site, by = c('uniprot_iso', 'PeptideSequence' = 'peptide'))

## MSstatsPTM Format
## Add missing columns
df_site$FragmentIon <- NA
df_prot$FragmentIon <- NA
df_site$ProductCharge <- NA
df_prot$ProductCharge <- NA
df_site$IsotopeLabelType <- 'L'
df_prot$IsotopeLabelType <- 'L'

## Combine proteins and sites so that we can estimate PTMs
df_site$ProteinName <- paste(df_site$uniprot_iso, df_site$site, sep = '_')
drops <- c("uniprot_iso", "site")
df_site <- df_site[ , !(names(df_site) %in% drops)]

## Rename columns
colnames(df_site) <- c('PeptideSequence', 'PrecursorCharge', 'is_mod', 'PSM',
                       'Condition', 'Run', 'BioReplicate', 'Intensity',
                       'FragmentIon', 'ProductCharge', 'IsotopeLabelType',
                       'ProteinName')
colnames(df_prot) <- c('ProteinName', 'PeptideSequence', 'PrecursorCharge',
                       'is_mod', 'PSM', 'Condition', 'Run', 'BioReplicate',
                       'Intensity', 'FragmentIon', 'ProductCharge',
                       'IsotopeLabelType')

# Visualization ----------------------------------------------------------------

# QC plot (boxplot)
df_prot %>%
  ggplot(aes(Run, log2(Intensity))) +
  geom_boxplot() +
  geom_vline(xintercept = c(2.5, 4.5, 6.5), linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Run", y = "Log2-intensity", title = "Global profiling")

df_prot %>%
  semi_join(df_fasta %>% rename(ProteinName = uniprot_iso) %>%
              filter(str_detect(entry_name, "HUMAN"))) %>%
  ggplot(aes(Run, log2(Intensity))) +
  geom_boxplot() +
  geom_vline(xintercept = c(2.5, 4.5, 6.5), linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Run", y = "Log2-intensity",
       title = "Global profiling - human proteins")

df_site %>%
  mutate(is_mod_fac = factor(ifelse(is_mod, "Modified", "Unmodified"))) %>%
  ggplot(aes(Run, log2(Intensity))) +
  geom_boxplot() +
  facet_wrap(~ is_mod_fac) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Run", y = "Log2-intensity", title = "KGG enrichment")

# Number of modified/unmodified peptides
df_site %>%
  distinct(Condition, Run, PeptideSequence, is_mod) %>%
  group_by(Condition, Run) %>%
  summarise(Modified = sum(is_mod), Unmodified = sum(!is_mod)) %>%
  ungroup() %>%
  gather(mod, nb_peptide, Modified:Unmodified) %>%
  ggplot(aes(Condition, nb_peptide)) +
  geom_jitter(aes(shape = mod), width = 0.1, size = 4) +
  theme_bw() +
  labs(x = "Group", y = "Number of peptides", title = "Human proteome")

aqua %>%
  ggplot(aes(Run, log2inty, color = PeptideSequence, group = PeptideSequence)) +
  geom_point() +
  geom_line() +
  labs(x = "", y = "Log2-intensity", title = "AQUA peptides", color = "") +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Apply Normalization ---------------------------------------------------------
prot_med <- df_prot %>% mutate(log2inty = log2(Intensity)) %>%
  group_by(Condition, Run) %>%
  summarise(log2inty_med = median(log2inty, na.rm = TRUE)) %>%
  mutate(
    log2inty_expected = median(log2inty_med),
    adjLog2inty = log2inty_expected - log2inty_med
  ) %>%
  ungroup()

# Modified peptides - applying the same adjustment to AQUA peptides as well
kgg_med <- df_site %>% mutate(log2inty = log2(Intensity)) %>%
  group_by(Run) %>%
  summarise(log2inty_med = median(log2inty, na.rm = TRUE))

log2inty_ref_kgg <- kgg_med$log2inty_med[kgg_med$Run == "mix1-T1"]

kgg_med <- kgg_med %>%
  mutate(log2inty_expected = ifelse(
    str_detect(Run, "mix1|mix2"), log2inty_ref_kgg, log2inty_ref_kgg - 1
  )) %>%
  mutate(adjLog2inty = log2inty_expected - log2inty_med)

refs <- list(
  PTM = kgg_med %>% select(Run, adjLog2inty),
  PROTEIN = prot_med %>% select(Run, adjLog2inty)
)

Reference_Norm <- function(df, ref) {
  df_aug <- left_join(df, ref)
  df_aug$log2inty = log2(df_aug$Intensity)
  df_aug$log2inty <- df_aug$log2inty + df_aug$adjLog2inty
  df_aug$Intensity = 2^df_aug$log2inty
  df_aug = df_aug %>% select(-log2inty, -adjLog2inty)
  return(df_aug)
}

df_site = Reference_Norm(df_site, refs$PTM)
df_prot = Reference_Norm(df_prot, refs$PROTEIN)

## Clean Data ------------------------------------------------------------------
df_site <- df_site %>% filter(is_mod == TRUE)
df_prot <- df_prot %>% filter(!is_mod)

df_site <- df_site[ , !(names(df_site) %in% c('is_mod'))]
df_prot <- df_prot[ , !(names(df_prot) %in% c('is_mod'))]

df_site <- .makeBalancedDesign(df_site, fill_missing = TRUE)
df_prot <- .makeBalancedDesign(df_prot, fill_missing = TRUE)

# AQUA peptides ----------------------------------------------------------------
aqua <- aqua %>% mutate(log2inty = log2(Intensity),
                        Run = paste(group, TechReplicate, sep = "-")) %>%
  left_join(kgg_med %>% select(Run, adjLog2inty)) %>%
  mutate(log2inty = log2inty + adjLog2inty) %>%
  mutate(Intensity = 2^log2inty) %>%
  select(-adjLog2inty, -log2inty)

aqua_seq <- readxl::read_xls("parkincancerKGG50mix_Masses.xls") %>%
  select(entry_name = Uniprot_Symbol, peptide = Sequence) %>%
  mutate(peptide_annot = str_to_upper(peptide) %>%
           str_replace("\\(GG\\)", "\\*"))

aqua_seq <- aqua_seq %>%
  mutate(entry_name = ifelse(entry_name == "EF1A1+_HUMAN", "EF1A1_HUMAN",
                             entry_name)) %>%
  left_join(df_fasta %>% select(uniprot_ac, uniprot_iso, entry_name)) %>%
  filter(uniprot_ac == uniprot_iso)

aqua_site <- locatePTM(
  aqua_seq$peptide_annot, aqua_seq$uniprot_iso,
  df_fasta, mod_residue, mod_symbol) %>%
  mutate(site = ifelse(site == "None", str_c("None_", 1:n()), site)) %>%
  mutate(site = str_c(site, " (heavy)")) %>%
  rename(peptide_annot = peptide) %>%
  left_join(aqua_seq %>% select(peptide, peptide_annot))

aqua_summarized <- aqua %>%
  left_join(aqua_site, by = c('PeptideSequence' = 'peptide')) %>%
  select(protein = uniprot_iso, site, PeptideSequence, group,
         TechReplicate, BioReplicate, Run, Intensity)

aqua_summarized$PrecursorCharge <- 3
aqua_summarized$IsotopeLabelType <- 'L'
aqua_summarized$FragmentIon <- NA
aqua_summarized$ProductCharge <- NA
aqua_summarized$PSM <- paste(aqua_summarized$PeptideSequence,
                             aqua_summarized$PrecursorCharge, sep = '_')
aqua_summarized$Condition <- aqua_summarized$group
aqua_summarized$ProteinName <- paste(aqua_summarized$protein,
                                     aqua_summarized$site, sep = "_")

aqua_summarized <- aqua_summarized %>% select(PeptideSequence,ProteinName,
                                              Condition, BioReplicate, Intensity)

aqua_summarized <- aqua_summarized %>%
  left_join(df_site %>% distinct(BioReplicate, Run), by = "BioReplicate")

## Summarize and model ---------------------------------------------------------
## Combine into one list
aqua_df <- list(
  PTM = aqua_summarized,
  PROTEIN = NULL
)

df_site <- df_site %>% select(-feature, -Fraction)

input_df <- list(
  # PTM = rbindlist(list(df_site,aqua_summarized), use.names=TRUE, fill = TRUE),
  PTM = df_site,
  PROTEIN = df_prot
)

summary_df <- dataSummarizationPTM(input_df,
                                   normalization = FALSE,
                                   normalization.PTM = FALSE,
                                   MBimpute = FALSE, MBimpute.PTM = FALSE)

aqua_summarized$Protein  = aqua_summarized$ProteinName
aqua_summarized$originalRUN = aqua_summarized$Run
aqua_summarized$GROUP = aqua_summarized$Condition
aqua_summarized$SUBJECT = aqua_summarized$BioReplicate
aqua_summarized$LogIntensities = log2(aqua_summarized$Intensity)
aqua_summarized$TotalGroupMeasurements = 4
aqua_summarized$NumMeasuredFeature = 4
aqua_summarized$MissingPercentage = 0
aqua_summarized$more50missing = FALSE
aqua_summarized$NumImputedFeature = 0
aqua_summarized = left_join(aqua_summarized, summary_df$PTM$ProteinLevelData %>% distinct(originalRUN, RUN),
                            by = "originalRUN")
aqua_summarized = aqua_summarized %>% select(-ProteinName, -Run, -Condition,
                                             -BioReplicate, -Intensity,
                                             -PeptideSequence)

summary_df$PTM$ProteinLevelData = rbindlist(
  list(summary_df$PTM$ProteinLevelData,aqua_summarized),
  use.names=TRUE, fill = TRUE)

save(summary_df,
     file = "D:\\Northeastern\\Research\\MSstats\\PTM_Refractor\\Label_Free_PTM_Data\\benchmark-spike\\MSstatsPTM_Summarized.rda")

# summary_df <- temp_summary_df
# summary_df$PTM$ProteinLevelData <- summary_df$PTM$ProteinLevelData %>%
#   mutate(LogIntensities = ifelse(str_detect(originalRUN, "mix3|mix4"),
#                                  LogIntensities -1, LogIntensities))

# Contrast
comparison = matrix(c(-1,1,0,0,
                       -1,0,1,0,
                       -1,0,0,1,
                       0,-1,1,0,
                       0,-1,0,1,
                       0,0,-1,1),nrow=6, byrow=TRUE)
row.names(comparison) = c("mix2 vs mix1", "mix3 vs mix1", "mix4 vs mix1",
                           "mix3 vs mix2", "mix4 vs mix2", "mix4 vs mix3")
colnames(comparison) = c("mix1", "mix2", "mix3", "mix4")

## Model
model_df <- groupComparisonPTM(summary_df, data.type = "LabelFree",
                               contrast.matrix = comparison,
                               use_log_file=FALSE)

save(model_df,
     file = "D:\\Northeastern\\Research\\MSstats\\MSstatsPTM-manuscript\\data\\Spike_in_MSstatsPTM_Model.rda")

# model_df$PTM.Model <- model_df$PTM.Model %>% mutate(
#   Label = ifelse(Label == "mix2 vs mix1", "mix1-mix2",
#                ifelse(Label == "mix3 vs mix1", "mix1-mix3",
#                       ifelse(Label == "mix4 vs mix1", "mix1-mix4",
#                              ifelse(Label == "mix2 vs mix3", "mix2-mix3",
#                                     ifelse(Label == "mix2 vs mix4", "mix2-mix4",
#                                            "mix3-mix4"))))))
#
# model_df$ADJUSTED.Model <- model_df$ADJUSTED.Model %>%
#   mutate(Label = ifelse(Label == "mix2 vs mix1", "mix1-mix2",
#               ifelse(Label == "mix3 vs mix1", "mix1-mix3",
#                      ifelse(Label == "mix4 vs mix1", "mix1-mix4",
#                             ifelse(Label == "mix2 vs mix3", "mix2-mix3",
#                                    ifelse(Label == "mix2 vs mix4", "mix2-mix4",
#                                           "mix3-mix4"))))))
# model_df$PROTEIN.Model <- model_df$PROTEIN.Model %>%
#   mutate(Label = ifelse(Label == "mix2 vs mix1", "mix1-mix2",
#               ifelse(Label == "mix3 vs mix1", "mix1-mix3",
#                      ifelse(Label == "mix4 vs mix1", "mix1-mix4",
#                             ifelse(Label == "mix2 vs mix3", "mix2-mix3",
#                                    ifelse(Label == "mix2 vs mix4", "mix2-mix4",
#                                           "mix3-mix4"))))))


## Model Plots -----------------------------------------------------------------
# Annotation
model_df$ADJUSTED.Model <- model_df$ADJUSTED.Model %>%
  mutate(labeling = ifelse(str_detect(Protein, "heavy"), "heavy", "light"))
model_df$PTM.Model <- model_df$PTM.Model %>%
  mutate(labeling = ifelse(str_detect(Protein, "heavy"), "heavy", "light"))

ground <- tibble(
  Label = c("mix2 vs mix1", "mix3 vs mix1", "mix4 vs mix1", "mix3 vs mix2",
            "mix4 vs mix2", "mix4 vs mix3"),
  trueLog2FC = c(-1, 1, 0, 2, 1, -1)
)

model_df$ADJUSTED.Model <- model_df$ADJUSTED.Model %>%
  left_join(ground) %>%
  mutate(contrast_truth = str_c(Label, " (", trueLog2FC, ")"))
model_df$ADJUSTED.Model$protadj <- "protein adjustment"
model_df$PTM.Model <- model_df$PTM.Model %>%
  left_join(ground) %>%
  mutate(contrast_truth = str_c(Label, " (", trueLog2FC, ")"))
model_df$PTM.Model$protadj <- "no adjustment"

temp_test <- rbind(model_df$ADJUSTED.Model %>% select(-GlobalProtein),
              model_df$PTM.Model %>%
                select(-c(MissingPercentage,ImputationPercentage, issue)))
temp_test$protadj <- factor(temp_test$protadj,
                            levels = c('no adjustment','protein adjustment'))

# Volcano plot
temp_test %>% arrange(desc(labeling)) %>%
  filter(is.finite(log2FC)) %>%
  ggplot(aes(log2FC, -log10(adj.pvalue))) +
  geom_point(aes(color = labeling, alpha = labeling)) +
  scale_color_manual(values = c("#D55E00", "#999999")) +
  scale_alpha_manual(values = c(1, 0.1)) +
  scale_size_manual(values = c(2, 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), color = "darkred") +
  geom_vline(xintercept = c(1, -1), color = "darkred", linetype = "dashed") +
  facet_grid(protadj ~ contrast_truth) +
  theme_bw() +
  theme(legend.position = "bottom")


## Comparison with other methods -----------------------------------------------
## Limma
library(limma)
## Fit limma model (both not and with adjust) for given dataset
fit_limma <- function(summarized_data, conditions, runs){

  ## Convert into format required for limma
  input <- data.frame(summarized_data %>% select(PTM, Run,
                                                 Condition, Abundance.x) %>%
                        pivot_wider(names_from = c(Condition, Run),
                                    values_from = Abundance.x,
                                    names_sort = TRUE))
  rownames(input) <- input$PTM
  input <- input %>% select(-PTM)
  input_adj <- data.frame(summarized_data %>%
                            select(PTM, Run, Condition, Adj_Abundance) %>%
                            pivot_wider(names_from = c(Condition, Run),
                                        values_from = Adj_Abundance,
                                        names_sort = TRUE))
  rownames(input_adj) <- input_adj$PTM
  input_adj <- input_adj %>% select(-PTM)

  ## Create contrast matrix
  class <- c()
  for (x in seq_len(conditions)){
    cond <- rep(paste0("mix", as.character(x)), runs)
    class <- c(class, cond)
  }
  class <- as.factor(class)
  design <- model.matrix(~0+class)
  colnames(design) <- c("mix1", "mix2", "mix3", "mix4")

  input.matrix <- as.matrix(input)
  input_adj.matrix <- as.matrix(input_adj)

  ## Run models
  fit <- lmFit(input.matrix, design=design)
  fit_adj <- lmFit(input_adj.matrix, design=design)

  # contrast.matrix <- design.pairs(colnames(design))
  contrast.matrix <- makeContrasts(mix2-mix1, mix3-mix1,
                                   mix4-mix1, mix3-mix2,
                                   mix4-mix2, mix4-mix3, levels=design)
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

input_df <- list(
  PTM = rbindlist(list(df_site,aqua_summarized), use.names=TRUE, fill=TRUE),
  PROTEIN = df_prot
)

ptm_df <- input_df$PTM
protein_df <- input_df$PROTEIN
df_site_limma <- df_site

summarized_ptm <- data.table()
summarized_proteins <- data.table()
runs <- unique(ptm_df$Run)
for (r in seq_along(runs)){
  sum_runs <- ptm_df %>% filter(Run == runs[[r]]) %>%
    group_by(ProteinName, Condition, Run) %>%
    summarize(Abundance = log2(sum(Intensity, na.rm = TRUE)))
  summarized_ptm <- rbindlist(list(summarized_ptm, sum_runs))

  sum_runs_prot <- protein_df %>% filter(Run == runs[[r]]) %>%
    group_by(ProteinName, Condition, Run) %>%
    summarize(Abundance = log2(sum(Intensity, na.rm = TRUE)))
  summarized_proteins <- rbindlist(list(summarized_proteins, sum_runs_prot))
}

# summarized_ptm <- summarized_ptm %>%
#   mutate(Abundance = ifelse(str_detect(Run, "mix3|mix4"),
#                             Abundance - 1, Abundance))

summarized_ptm$PTM <- summarized_ptm$ProteinName
summarized_ptm$ProteinName <- sapply(summarized_ptm$PTM,
                                     function(x) {str_split(x, "_",2)[[1]][1]})
joined <- merge(summarized_ptm, summarized_proteins,
                by = c("ProteinName", "Run", "Condition"), all.x = TRUE)
joined$Adj_Abundance <- joined$Abundance.x - joined$Abundance.y


limma_test_res <- fit_limma(joined %>% filter(is.finite(Abundance.x)), 4, 2)

limma_test_res$limma_adj_test <- limma_test_res$limma_adj_test %>%
  mutate(Label = ifelse(Label == "mix2 - mix1", "mix2 vs mix1",
      ifelse(Label == "mix3 - mix1", "mix3 vs mix1",
             ifelse(Label == "mix4 - mix1", "mix4 vs mix1",
                    ifelse(Label == "mix3 - mix2", "mix3 vs mix2",
                           ifelse(Label == "mix4 - mix2", "mix4 vs mix2",
                                  "mix4 vs mix3"))))))

limma_test_res$limma_test <- limma_test_res$limma_test %>%
  mutate(Label = ifelse(Label == "mix2 - mix1", "mix2 vs mix1",
                        ifelse(Label == "mix3 - mix1", "mix3 vs mix1",
                               ifelse(Label == "mix4 - mix1", "mix4 vs mix1",
                                      ifelse(Label == "mix3 - mix2", "mix3 vs mix2",
                                             ifelse(Label == "mix4 - mix2", "mix4 vs mix2",
                                                    "mix4 vs mix3"))))))

## Add spike-in id
limma_test_res$limma_test <- limma_test_res$limma_test %>%
  mutate(labeling = ifelse(str_detect(PTM, "heavy"), "heavy", "light"))
limma_test_res$limma_adj_test <- limma_test_res$limma_adj_test %>%
  mutate(labeling = ifelse(str_detect(PTM, "heavy"), "heavy", "light"))

ground <- tibble(
  Label = c("mix2 vs mix1", "mix3 vs mix1", "mix4 vs mix1", "mix3 vs mix2",
            "mix4 vs mix2", "mix4 vs mix3"),
  trueLog2FC = c(-1, 1, 0, 2, 1, -1)
)

limma_test_res$limma_test$adj.pvalue <- 0.
limma_test_res$limma_adj_test$adj.pvalue <- 0.

## adjust pvalue
limma_test_res$limma_test[, adj.pvalue := p.adjust(pvalue, method = "BH"),
                        by = "Label"]
limma_test_res$limma_adj_test[, adj.pvalue := p.adjust(pvalue, method = "BH"),
                          by = "Label"]

limma_test_res$limma_test <- limma_test_res$limma_test %>%
  left_join(ground) %>%
  mutate(contrast_truth = str_c(Label, " (", trueLog2FC, ")"))
limma_test_res$limma_test$protadj <- "no adjustment"
limma_test_res$limma_adj_test <- limma_test_res$limma_adj_test %>%
  left_join(ground) %>%
  mutate(contrast_truth = str_c(Label, " (", trueLog2FC, ")"))
limma_test_res$limma_adj_test$protadj <- "protein adjustment"

temp_test <- rbind(limma_test_res$limma_test, limma_test_res$limma_adj_test)
temp_test$protadj <- factor(temp_test$protadj, levels = c('no adjustment','protein adjustment'))

# Volcano plot
temp_test %>% arrange(desc(labeling)) %>%
  filter(abs(Log2FC) < Inf & Label %in% c("mix1-mix3", "mix1-mix4", "mix2-mix3", "mix2-mix4")) %>%
  ggplot(aes(Log2FC, -log10(adj.pvalue))) +
  geom_point(aes(color = labeling, alpha = labeling)) +
  scale_color_manual(values = c("#D55E00", "#999999")) +
  scale_alpha_manual(values = c(1, 0.1)) +
  scale_size_manual(values = c(2, 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), color = "darkred") +
  geom_vline(xintercept = c(1, -1), color = "darkred", linetype = "dashed") +
  facet_grid(protadj ~ contrast_truth) +
  theme_bw() +
  theme(legend.position = "bottom")

limma_model <- temp_test
save(limma_model, file = "D:\\Northeastern\\Research\\MSstats\\MSstatsPTM-manuscript\\data\\Spike_in_Limma_Model.rda")



## Ttest -----------------------------------------------------------------------
## Run ttest using run level summared data
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
        tryCatch({ttest_ptm <- t.test((temp_joined %>%
                                         filter(Condition == groups[g]) %>%
                                         select(Abundance.x))[[1]],
                                      (temp_joined %>%
                                         filter(Condition == groups[g2]) %>%
                                         select(Abundance.x))[[1]])
                  ttest_temp <- rbindlist(list(ttest_temp, data.table(ptm = ptms[[p]],
                                          label = paste(groups[g], "vs",
                                                        groups[g2], sep = " "),
                                          pval = ttest_ptm$p.value,
                                          tstat = ttest_ptm$statistic[[1]],
                                          SE = ttest_ptm$stderr,
                                          df = ttest_ptm$parameter[[1]],
                                          estimate = ttest_ptm$estimate[[2]] - ttest_ptm$estimate[[1]])))
                  },
                 error=function(e){cat("ERROR :", as.character(p), "\n")}
                 )
        tryCatch({ttest_adj_ptm <- t.test((temp_joined %>%
                                             filter(Condition == groups[g]) %>%
                                             select(Adj_Abundance))[[1]],
                                          (temp_joined %>%
                                             filter(Condition == groups[g2]) %>%
                                             select(Adj_Abundance))[[1]])
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
save(ttest_model, file = "D:\\Northeastern\\Research\\MSstats\\MSstatsPTM-manuscript\\data\\Spike_in_ttest.rda")

ttest_model <- ttest_modelkeep
ttest_model$ttest_adj_temp <- ttest_model$ttest_adj_temp %>% filter(pval != 1)
ttest_model$ttest_temp <- ttest_model$ttest_temp %>% filter(pval != 1)

ttest_model$ttest_adj_temp <- ttest_model$ttest_adj_temp %>%
  mutate(label = ifelse(label == "mix1 vs mix2", "mix1-mix2",
              ifelse(label == "mix1 vs mix3", "mix1-mix3",
                     ifelse(label == "mix1 vs mix4", "mix1-mix4",
                            ifelse(label == "mix2 vs mix3", "mix2-mix3",
                                   ifelse(label == "mix2 vs mix4", "mix2-mix4",
                                          "mix3-mix4"))))))
ttest_model$ttest_temp <- ttest_model$ttest_temp %>%
  mutate(label = ifelse(label == "mix1 vs mix2", "mix1-mix2",
              ifelse(label == "mix1 vs mix3", "mix1-mix3",
                     ifelse(label == "mix1 vs mix4", "mix1-mix4",
                            ifelse(label == "mix2 vs mix3", "mix2-mix3",
                                   ifelse(label == "mix2 vs mix4", "mix2-mix4",
                                          "mix3-mix4"))))))

## Fix comparison direction
ttest_model$ttest_adj_temp[label == "mix2-mix3"]$estimate <- ttest_model$ttest_adj_temp[label == "mix2-mix3"]$estimate*-1
ttest_model$ttest_adj_temp[label == "mix2-mix4"]$estimate <- ttest_model$ttest_adj_temp[label == "mix2-mix4"]$estimate*-1
ttest_model$ttest_adj_temp[label == "mix3-mix4"]$estimate <- ttest_model$ttest_adj_temp[label == "mix3-mix4"]$estimate*-1
ttest_model$ttest_temp[label == "mix2-mix3"]$estimate <- ttest_model$ttest_temp[label == "mix2-mix3"]$estimate*-1
ttest_model$ttest_temp[label == "mix2-mix4"]$estimate <- ttest_model$ttest_temp[label == "mix2-mix4"]$estimate*-1
ttest_model$ttest_temp[label == "mix3-mix4"]$estimate <- ttest_model$ttest_temp[label == "mix3-mix4"]$estimate*-1

ttest_model$ttest_adj_temp <- ttest_model$ttest_adj_temp %>%
  mutate(labeling = ifelse(str_detect(ptm, "heavy"), "heavy", "light"))
ttest_model$ttest_temp <- ttest_model$ttest_temp %>%
  mutate(labeling = ifelse(str_detect(ptm, "heavy"), "heavy", "light"))

ground <- tibble(
  label = c("mix1-mix2", "mix1-mix3", "mix1-mix4", "mix2-mix3",
            "mix2-mix4", "mix3-mix4"),
  trueLog2FC = c(-1, 1, 0, -2, -1, 1)
)

ttest_model$ttest_adj_temp$adj.pvalue <- 0.
ttest_model$ttest_temp$adj.pvalue <- 0.

for (l in seq_along(unique(ttest_model$ttest_adj_temp$label))){
  ttest_model$ttest_adj_temp[label == unique(ttest_model$ttest_adj_temp$label)[l]]$adj.pvalue <- p.adjust(
    ttest_model$ttest_adj_temp[label == unique(ttest_model$ttest_adj_temp$label)[l]]$pval, method = "BH")
  ttest_model$ttest_temp[label == unique(ttest_model$ttest_temp$label)[l]]$adj.pvalue <- p.adjust(
    ttest_model$ttest_temp[label == unique(ttest_model$ttest_temp$label)[l]]$pval, method = "BH")
}

ttest_model$ttest_temp <- ttest_model$ttest_temp %>%
  left_join(ground) %>%
  mutate(contrast_truth = str_c(label, " (", trueLog2FC, ")"))
ttest_model$ttest_temp$protadj <- "no adjustment"
ttest_model$ttest_adj_temp <- ttest_model$ttest_adj_temp %>%
  left_join(ground) %>%
  mutate(contrast_truth = str_c(label, " (", trueLog2FC, ")"))
ttest_model$ttest_adj_temp$protadj <- "protein adjustment"

temp_test <- rbind(ttest_model$ttest_temp, ttest_model$ttest_adj_temp)
temp_test$protadj <- factor(temp_test$protadj, levels = c('no adjustment','protein adjustment'))

ttest_final_model <- temp_test
save(ttest_final_model, file = "D:\\Northeastern\\Research\\MSstats\\MSstatsPTM-manuscript\\data\\Spike_in_ttest_model.rda")

# Volcano plot
temp_test %>% arrange(desc(labeling)) %>%
  filter(abs(estimate) < Inf & label %in% c("mix1-mix3", "mix1-mix4", "mix2-mix3", "mix2-mix4")) %>%
  ggplot(aes(estimate, -log10(adj.pvalue))) +
  geom_point(aes(color = labeling, alpha = labeling)) +
  scale_color_manual(values = c("#D55E00", "#999999")) +
  scale_alpha_manual(values = c(1, 0.1)) +
  scale_size_manual(values = c(2, 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = -log10(0.05), color = "darkred") +
  geom_vline(xintercept = c(1, -1), color = "darkred", linetype = "dashed") +
  facet_grid(protadj ~ contrast_truth) +
  theme_bw() +
  theme(legend.position = "bottom")


## Anova -----------------------------------------------------------------------

fit_anova = function(data, contrast){

  ptms = unique(data$PTM)
  test_results = data.table()
  unadj_test_results = data.table()
  for (i in seq_along(ptms)){#
    print(ptms[[i]])
    ptm_data_adj = data %>% filter(PTM == ptms[[i]] & is.finite(Adj_Abundance))
    ptm_data_unadj = data %>% filter(PTM == ptms[[i]] & is.finite(Abundance.x))

    if (length(unique(ptm_data_adj$Condition)) < 4){
      model_conditions = unique(ptm_data_adj$Condition)
      missing = setdiff(c("mix1", "mix2", "mix3", "mix4"),
                        unique(ptm_data_adj$Condition))
      keep = c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
      for (k in seq_along(missing)){
        check = !str_detect(rownames(contrast), missing[[k]])
        keep = ifelse(!check, check, keep)
        temp_contrast = contrast[keep,, drop = FALSE]
      }
      # print(temp_contrast)
      # print(model_conditions)

      temp_contrast = temp_contrast[, model_conditions]

    } else {
      temp_contrast = contrast
    }

    tryCatch({

      ## Adj model
      res_aov_adj = lm(Adj_Abundance~Condition, data = ptm_data_adj)
      model_summary_adj = summary(res_aov_adj)

      coefs_adj = res_aov_adj$coefficients
      coefs_adj[1] = 0
      df_adj = res_aov_adj$df.residual
      sigma_adj = model_summary_adj$sigma

      ## unadj model
      res_aov_unadj = lm(Abundance.x~Condition, data = ptm_data_unadj)
      model_summary_unadj = summary(res_aov_unadj)

      coefs_unadj = res_aov_unadj$coefficients
      coefs_unadj[1] = 0
      df_unadj = res_aov_unadj$df.residual
      sigma_unadj = model_summary_unadj$sigma

      for (j in seq_len(nrow(temp_contrast))){

        comp = temp_contrast[j,]
        label = rownames(temp_contrast)[j]

        ## Adj
        beta_adj = sum(coefs_adj*comp)
        tstat_adj = beta_adj / sigma_adj

        pval_adj = 2*pt(-abs(tstat_adj),df=df_adj)

        test_results = rbindlist(list(test_results,
                                  data.table(ptm = ptms[[i]], label = label, log2FC = beta_adj,
                        se = sigma_adj, tstat = tstat_adj, pvalue = pval_adj)))

        ## unadj
        beta_unadj = sum(coefs_unadj*comp)
        tstat_unadj = beta_unadj / sigma_unadj

        pval_unadj = 2*pt(-abs(tstat_unadj),df=df_unadj)

        unadj_test_results = rbindlist(list(unadj_test_results,
                                        data.table(ptm = ptms[[i]], label = label, log2FC = beta_unadj,
                        se = sigma_unadj, tstat = tstat_unadj, pvalue = pval_unadj)))

      }

    },
    error=function(e){cat("ERROR : fitting error\n")})
  }

  results = list("anova_adj" = test_results,
                 "anova_unadj" = unadj_test_results)
  return(results)
}

anova_results = fit_anova(joined, comparison)
keep = anova_results
anova_results$anova_adj$protadj <- "protein adjustment"
anova_results$anova_unadj$protadj <- "no adjustment"

anova_results$anova_adj[, adj.pvalue := p.adjust(pvalue, method = "BH"),
            by = "label"]
anova_results$anova_unadj[, adj.pvalue := p.adjust(pvalue, method = "BH"),
                        by = "label"]

anova_results$anova_adj <- anova_results$anova_adj %>%
  mutate(labeling = ifelse(str_detect(ptm, "heavy"), "heavy", "light"))
anova_results$anova_unadj <- anova_results$anova_unadj %>%
  mutate(labeling = ifelse(str_detect(ptm, "heavy"), "heavy", "light"))

anova_results = rbind(anova_results$anova_adj,
                      anova_results$anova_unadj)

save(anova_results, file = "D:\\Northeastern\\Research\\MSstats\\MSstatsPTM-manuscript\\data\\anova_results.rda")
