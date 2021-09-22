

## Code directly from MSstats designSampleSize function
## Range of desired FC
delta <- log2(seq(1.5, 2.5, 0.025))
desiredFC <- 2 ^ delta

# delta <- seq(.5, 1.5, 0.025)
# desiredFC <- delta

## Define base metrics
FDR <- .05
power <- .8

## function calculations
m0_m1 <- 99 ## it means m0/m1=99, m0/(m0+m1)=0.99
alpha <- power * FDR / (1 + (1 - FDR) * m0_m1)

z_alpha <- qnorm(1 - alpha / 2)
z_beta <- qnorm(power)
aa <- (delta / (z_alpha + z_beta)) ^2

#Specify sigma for ptm and unmod pep


## Original code containing sigma for subject.. can include or assume 0 for now
## numSample = round(2 * (median_sigma_error + median_sigma_subject) / aa, 0)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

numSample_nomod_1 <- round(2 * (.1) / aa, 0)
numSample_nomod_2 <- round(2 * (.2) / aa, 0)
numSample_nomod_3 <- round(2 * (.4) / aa, 0)

p1 <- data.frame("SD" = c(rep(".1", length(numSample_nomod_1)),
                      rep(".2", length(numSample_nomod_1)),
                      rep(".4", length(numSample_nomod_1))),
           "Log2FC" = c(desiredFC, desiredFC, desiredFC),
           "num_samples" = c(numSample_nomod_1, numSample_nomod_2, numSample_nomod_3)) %>%
  ggplot() + geom_line(aes(x = Log2FC, y = num_samples, group = SD, color = SD), size = 1.2) +
  scale_color_discrete(cbPalette) +
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
                          rep("Baseline - .2", length(numSample_nomod_1))),
                 "Log2FC" = c(desiredFC, desiredFC, desiredFC, desiredFC),
                 "num_samples" = c(numSample_nomod_1, numSample_nomod_2, numSample_nomod_3,
                                   numSample_nomod_4)) %>%
  ggplot() + geom_line(aes(x = Log2FC, y = num_samples, group = SD, color = SD)) +
  labs(title = "Replicates Needed to Reach Power of .9", y = "# of replicates per condition") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        title = element_text(size = 16),
        strip.text = element_text(size = 12))

p2
