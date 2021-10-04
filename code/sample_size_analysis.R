
library(tidyverse)

## Define standard deviation sampling
sample_variance <- function(mod_samples, unmod_samples, mod_features,
                            unmod_features, mod_variance, unmod_variance){

  final_sd <- numeric(100000)

  for (i in seq_len(100000)){

    replicate_means_ptms <- rnorm(mod_samples, mean = 0, mod_variance)
    replicate_observations_ptms <- numeric(length(mod_samples))

    # for (y in seq_along(replicate_means_ptms)){
    #   feature_samples <- rnorm(mod_features, mean = replicate_means_ptms[[y]], 2)
    #   replicate_observations_ptms[[y]] <- mean(feature_samples)
    # }

    replicate_means_prot <- rnorm(unmod_samples, mean = 0, unmod_variance)
    replicate_observations_prot <- numeric(length(unmod_samples))

    # for (y in seq_along(replicate_means_prot)){
    #   feature_samples_prot <- rnorm(unmod_features, mean = replicate_means_prot[[y]], 2)
    #   replicate_observations_prot[[y]] <- mean(feature_samples_prot)
    # }

    std_ptm <- mod_variance^2 / sum((replicate_means_ptms - mean(replicate_means_ptms))^2)
    std_prot <- unmod_variance^2 / sum((replicate_means_prot - mean(replicate_means_prot))^2)

    final_sd[[i]] <- sqrt(std_ptm^2 + std_prot^2)

  }

  return(mean(final_sd))
}



## Code directly from MSstats designSampleSize function
## Power analysis
FDR <- 0.05
m0_m1 <- 99
delta <- seq(.5, 1.5, 0.025)

low <- .1
high <- 3

n_sample_list <- c(3,5,10)
n_feature_list <- c(3,10)


var_list <- numeric(9)
i <- 1
for (x in seq_along(n_sample_list)){
  for (y in seq_along(n_sample_list)){
    var_list[i] <- sample_variance(n_sample_list[[x]], n_sample_list[[y]],
                                   10, 10, low, low)
    i <- i+1
  }
}

var_list <- c(sqrt((.1/3)^2 + (.1/3)^2), sqrt((.1/3)^2 + (.1/5)^2), sqrt((.1/3)^2 + (.1/10)^2),
              sqrt((.1/5)^2 + (.1/3)^2), sqrt((.1/5)^2 + (.1/5)^2), sqrt((.1/5)^2 + (.1/10)^2),
              sqrt((.1/10)^2 + (.1/3)^2), sqrt((.1/10)^2 + (.1/5)^2), sqrt((.1/10)^2 + (.1/10)^2))

t <- delta / var_list[[1]]

powerTemp <- seq(0, 1, 0.01)
power <- numeric(length(t))

for (i in seq_along(t)) {
  diff <- qnorm(powerTemp) + qnorm(1 - powerTemp * FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
  min(abs(diff), na.rm = TRUE)
  power[i] = powerTemp[order(abs(diff))][1]
}
power


p1 <- data.frame("Log2FC" = delta,
                 "Power" = power) %>%
  ggplot() + geom_line(aes(x = Log2FC, y = Power), size = 1.2) +
  #scale_colour_manual(values=cbPalette) +
  labs(title = "Power analysis with low replicate variance") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        title = element_text(size = 16),
        strip.text = element_text(size = 12))



## Original code containing sigma for subject.. can include or assume 0 for now
## numSample = round(2 * (median_sigma_error + median_sigma_subject) / aa, 0)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#999999", "#0072B2", "#D55E00", "#CC79A7", "#999999")

numSample_nomod_1 <- round(2 * (.1) / aa, 0)
numSample_nomod_2 <- round(2 * (.2) / aa, 0)
numSample_nomod_3 <- round(2 * (.4) / aa, 0)

p1 <- data.frame("SD" = c(rep("`.1", length(numSample_nomod_1)),
                      rep(".2", length(numSample_nomod_1)),
                      rep(".4", length(numSample_nomod_1))),
           "Log2FC" = c(desiredFC, desiredFC, desiredFC),
           "num_samples" = c(numSample_nomod_1, numSample_nomod_2, numSample_nomod_3)) %>%
  ggplot() + geom_line(aes(x = Log2FC, y = num_samples, group = SD, color = SD), size = 1.2) +
  scale_colour_manual(values=cbPalette) +
  labs(title = "Replicates Needed to Reach Power of .9", y = "# of replicates per condition") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        title = element_text(size = 16),
        strip.text = element_text(size = 12))
png("supplementary/sim_new/simple_power_analysis.png", width = 750, height = 500)
p1
dev.off()

numSample_nomod_1 <- round((2 * sqrt(.2^2 + .1^2)) /  aa, 0)
numSample_nomod_2 <- round((2 * sqrt(.2^2 + .2^2)) / aa, 0)
numSample_nomod_3 <- round((2 * sqrt(.2^2 + .4^2)) / aa, 0)
numSample_nomod_4 <- round((2 * .2) / aa, 0)

p2 <- data.frame("SD" = c(rep(".2, .1", length(numSample_nomod_1)),
                          rep(".2, .2", length(numSample_nomod_1)),
                          rep(".2, .4", length(numSample_nomod_1)),
                          rep("Baseline: .2", length(numSample_nomod_1))),
                 "Log2FC" = c(desiredFC, desiredFC, desiredFC, desiredFC),
                 "num_samples" = c(numSample_nomod_1, numSample_nomod_2, numSample_nomod_3,
                                   numSample_nomod_4)) %>%
  ggplot() + geom_line(aes(x = Log2FC, y = num_samples, group = SD, color = SD), size = 1.2) +
  scale_colour_manual(values=cbPalette) +
  labs(title = "Replicates Needed to Reach Power of .9 with Adjusted SD", y = "# of replicates per condition") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        title = element_text(size = 16),
        strip.text = element_text(size = 12))

png("supplementary/sim_new/power_analysis_sd_combined.png", width = 750, height = 500)
p2
dev.off()
