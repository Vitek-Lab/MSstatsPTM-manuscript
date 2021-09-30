

## Code directly from MSstats designSampleSize function
## Range of desired FC
delta <- log2(seq(1.5, 2.5, 0.025))
desiredFC <- 2 ^ delta

# delta <- seq(.5, 1.5, 0.025)
# desiredFC <- delta

## Define base metrics

low <- 1
high <- 3

mod_pep_sigma <- low
unmod_pep_sigma <- low

n_sample_list <- c(3,5,7,10)

mod_num_sample <- n_sample_list[[1]]
unmod_num_sample <- n_sample_list[[1]]

#Specify sigma for ptm and unmod pep

samples_mod <- rnorm(n_sample_list[[1]], mean = 0, low)
std_dev_mod <- sum((samples_mod - mean(samples_mod)^2)) / (mod_num_sample)
samples_unmod <- rnorm(n_sample_list[[1]], mean = 0, low)
std_dev_unmod <- sum((samples_unmod - mean(samples_unmod)^2)) / (unmod_num_sample)

t <- delta / sqrt(std_dev_mod^2 + std_dev_unmod^2)

FDR <- 0.05
m0_m1 <- 99
powerTemp <- seq(0, 1, 0.01)
power <- numeric(length(t))

for (i in seq_along(t)) {
  diff <- qnorm(powerTemp) + qnorm(1 - powerTemp * FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
  min(abs(diff), na.rm = TRUE)
  power[i] = powerTemp[order(abs(diff))][1]
}

















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
