
library(tidyverse)

## Calculations w/o feature variations -----------------------------------------

sigma_model <- function(ptm_var, ptm_n, prot_var, prot_n){
  simga <- sqrt((ptm_var/(.25*ptm_n)) + (prot_var/(.25*prot_n)))
  return(simga)
}

## Define adjustable metrics
low <- .1
high <- .3
n_sample_list <- c(3,5,10)

## Hard metrics
FDR <- 0.05
m0_m1 <- 99
delta <- seq(1, 2, 0.025)

var_list <- numeric(9)
i <- 1
for (x in seq_along(n_sample_list)){
  for (y in seq_along(n_sample_list)){
    var_list[i] <- sigma_model(low, n_sample_list[[x]], low, n_sample_list[[y]])
    i <- i+1
  }
}

t <- delta / var_list[[9]]

powerTemp <- seq(0, 1, 0.01)
power <- numeric(length(t))

for (i in seq_along(t)) {
  diff <- qnorm(powerTemp) + qnorm(1 - powerTemp * FDR / (1 + (1 - FDR) * m0_m1) / 2) - t[i]
  min(abs(diff), na.rm = TRUE)
  power[i] = powerTemp[order(abs(diff))][1]
}
power


## Simulate feature effect -----------------------------------------------------

n_feature_list <- 1:10

final_var <- c()

for (n in seq_along(n_feature_list)){
  test <- c()
  for (i in 1:1000) {
    r1 <- c()
    r2 <- c()

    for (x in seq_len(3)){

      t1 <- rnorm(1, mean = 0, .5)
      t2 <- rnorm(1, mean = 2, .5)

      for (y in seq_len(n_feature_list[[n]])){
        f1 <- rnorm(y, mean = t1, .001)
        f2 <- rnorm(y, mean = t2, .001)
      }

      r1 <- c(r1, log(sum(f1^2)))
      r2 <- c(r2, log(sum(f2^2)))

    }

    combined <- c(r1, r2)
    test_mean <- mean(combined)
    std_test <- sum((combined - test_mean)^2) / (length(combined) - 2)
    test <- c(test, std_test)
  }
  final_var <- c(final_var, mean(test))
}

a1 <- rnorm(1000, mean = 0, .5)
a2 <- rnorm(1000, mean = 2, .5)
acombined <- c(a1, a2)
atest_mean <- mean(acombined)
astd_test <- sum((acombined - atest_mean)^2) / (length(acombined) - 2)

percent_change <- c()
for (i in 2:length(final_var)){

  percent_change <- c(percent_change, (final_var[[i]] - final_var[[i-1]])/ final_var[[i-1]])
}

keep_1 <- percent_change

## Old charts ------------------------------------------------------------------


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


set.seed(1)
x<-c(0,0,2)
y<-3*x+5+rnorm(3,0,1)
lr<-linear.regression(x,y)
y.hat.x1 = lr[1]+lr[2]*1
halfWidth.c<-t95*sqrt(lr[3]*(1/10+(y.hat.x1-mean(x))ˆ2)/sum((x-mean(x))ˆ2))
halfWidth.p<-t95*sqrt(lr[3]*(1+1/10+(y.hat.x1-mean(x))ˆ2)/sum((x-mean(x))ˆ2))
l95c<-y.hat.x1-halfWidth.c
l95p<-y.hat.x1-halfWidth.p
u95c<-y.hat.x1+halfWidth.c
u95p<-y.hat.x1+halfWidth.p
dt.x1.95<-rbind(data.frame(lower = l95c,upper = u95c,type = "CI"),
                data.frame(lower = l95p,upper = u95p,type = "PI"))
ggplot()+geom_point(data = data.frame(x,y),aes(x = x,y = y))+
  geom_errorbar(data = dt.x1.95,mapping = aes(x = c(1,1.05),ymin = lower.y))

linear.regression <- function(x,y){
  mx<-mean(x)
  my<-mean(y)
  b1<-sum((x-mx)*(y-my))/sum((x-mx)ˆ2)
  b0<-my-mx*b1
  e<-y - (b0+b1*x)
  mse<-sum(eˆ2)/(length(x)-2)
  var.b1<-mse/sum((x-mx)ˆ2)
  t.b1<-b1/sqrt(var.b1)
  p.b1<-(1-pt(abs(t.b1),df = (length(y)-2)))*2 #two tails
  return(c(b0,b1,mse,var.b1,t.b1,p.b1))
  }
