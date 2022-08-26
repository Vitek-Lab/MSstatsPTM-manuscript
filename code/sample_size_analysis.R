
library(tidyverse)
library(ggthemes)
library(data.table)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Calculations w/o feature variations -----------------------------------------

sigma_model <- function(ptm_var, ptm_n, prot_var, prot_n){
  simga <- sqrt((ptm_var/ptm_n) + (prot_var/prot_n))
  return(simga)
}

## Define adjustable metrics
low <- .3
mid <- .45
high <- .45
n_sample_list <- c(2,5)

## Hard metrics
FDR <- 0.05
m0_m1 <- 99
delta <- seq(1, 2, 0.025)

bind_list <- c()

i <- 1
for (x in seq_along(n_sample_list)){
  for (y in seq_along(n_sample_list)){
    var <- sigma_model(high, n_sample_list[[x]], low, n_sample_list[[y]])
    
    
    t <- delta / var
    
    powerTemp <- seq(0, 1, 0.01)
    power <- numeric(length(t))
    
    for (i in seq_along(t)) {
      diff <- qnorm(powerTemp) + qnorm(1 - powerTemp * FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
      min(abs(diff), na.rm = TRUE)
      power[i] = powerTemp[order(abs(diff))][1]
    }
    
    temp <- data.table("PTM-Protein Replicates" = paste(n_sample_list[[x]],  n_sample_list[[y]], sep = "-"),
               "Log2FC" = delta, "Power" = power)
    
    bind_list <- c(bind_list, list(temp))
    
  }
}

low_var_df <- rbindlist(bind_list)
same_var_df <- rbindlist(bind_list)
high_ptm_var_df <- rbindlist(bind_list)
high_var_df <- rbindlist(bind_list)

## Charts ------------------------------------------------------------------
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#999999", 
               "#0072B2", "#D55E00", "#CC79A7", "#999999")

p1 <- same_var_df %>% mutate(log2FC = log2(Log2FC)) %>% 
  ggplot() + geom_line(aes(x = log2FC, y = Power, color = `PTM-Protein Replicates`), 
                       position = "jitter", size = 1.2) +
  scale_colour_manual(values=cbPalette) +
  labs(title = "Power analysis with same variance") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        title = element_text(size = 22),
        strip.text = element_text(size = 16)) + 
  ylim(0, 1)
p1

p2 <- high_ptm_var_df %>% mutate(log2FC = log2(Log2FC)) %>% 
  ggplot() + geom_line(aes(x = log2FC, y = Power, color = `PTM-Protein Replicates`), size = 1.2) +
  scale_colour_manual(values=cbPalette) +
  labs(title = "Power analysis with high PTM variance") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        title = element_text(size = 22),
        strip.text = element_text(size = 16)) + 
  ylim(0, 1)
p2

png("../supplementary/sim_new/same_var_power.png", width = 750, height = 500)
p1
dev.off()

png("../supplementary/sim_new/high_ptm_var_power.png", width = 750, height = 500)
p2
dev.off()


