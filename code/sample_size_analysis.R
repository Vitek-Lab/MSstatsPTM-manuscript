

## Code directly from MSstats designSampleSize function
## Range of desired FC
delta <- log2(seq(.5, 1.5, 0.025))
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
sigma_error_mod <- .2
sigma_error_unmod <- .1

## Original code containing sigma for subject.. can include or assume 0 for now
## numSample = round(2 * (median_sigma_error + median_sigma_subject) / aa, 0)

numSample_nomod <- round(2 * (sigma_error_mod) / aa, 0)
numSample_mod <- round(2 * (sqrt(sigma_error_mod^2 + sigma_error_unmod^2)) / aa, 0)

p1 <- data.frame("type" = c(rep("With adjustment", length(numSample_mod)),
                      rep("Without adjustment", length(numSample_nomod))),
           "Log2FC" = c(desiredFC, desiredFC),
           "num_samples" = c(numSample_mod, numSample_nomod)) %>%
  ggplot() + geom_line(aes(x = Log2FC, y = num_samples, group = type, color = type))

sigma_error_mod <- .2
sigma_error_unmod <- .2

numSample_nomod <- round(2 * (sigma_error_mod + .3) / aa, 0)
numSample_mod <- round(2 * (sqrt(sigma_error_mod^2 + sigma_error_unmod^2) + .3) / aa, 0)

p2 <- data.frame("type" = c(rep("With adjustment", length(numSample_mod)),
                      rep("Without adjustment", length(numSample_nomod))),
           "Log2FC" = c(delta, delta),
           "num_samples" = c(numSample_mod, numSample_nomod)) %>%
  ggplot() + geom_line(aes(x = Log2FC, y = num_samples, group = type, color = type))

sigma_error_mod <- .2
sigma_error_unmod <- .4

numSample_nomod <- round(2 * (sigma_error_mod + .05) / aa, 0)
numSample_mod <- round(2 * (sqrt(sigma_error_mod^2 + sigma_error_unmod^2) + .05) / aa, 0)

p3 <- data.frame("type" = c(rep("With adjustment", length(numSample_mod)),
                      rep("Without adjustment", length(numSample_nomod))),
           "Log2FC" = c(delta, delta),
           "num_samples" = c(numSample_mod, numSample_nomod)) %>%
  ggplot() + geom_line(aes(x = Log2FC, y = num_samples, group = type, color = type))

