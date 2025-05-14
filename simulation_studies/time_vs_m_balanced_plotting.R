################################################################################
# Plot runtime results for parameter estimation for a balanced random effects design.
################################################################################

library(ggplot2)
library(ggh4x)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

Gaussian_results <- readRDS("./time_vs_m_balanced_Gaussian.rds")
Bernoulli_results <- readRDS("./time_vs_m_balanced_Bernoulli.rds")

Gaussian_results$likelihood <- "Gaussian"
Bernoulli_results$likelihood <- "Bernoulli"
results <- rbind(Gaussian_results, Bernoulli_results)

#Remove not evaluated combinations
results <- results[!is.na(results$time_estimation),]

results$method <- as.factor(results$method)
results$likelihood <- as.factor(results$likelihood)

ggplot(data=results, aes(x=m, y=time_estimation, color=method, shape=method)) +
  geom_line(linewidth=1) + geom_point(size=2) + 
  facet_nested_wrap(~likelihood, nrow=1, drop=T, scales = "free") +
  xlab("m") + ylab("Time (s)") + theme_bw() + 
  theme(legend.position = "top", legend.title=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_brewer(type = "qual", palette=6) +
  facetted_pos_scales(y = list(likelihood == "Gaussian" ~ scale_y_continuous(trans = "log1p", labels=function(x) format(x, big.mark = "'", scientific = FALSE)),
                               likelihood == "Bernoulli" ~ scale_y_continuous(trans = "log1p", labels=function(x) format(x, big.mark = "'", scientific = FALSE))),
                      x = list(likelihood == "Gaussian" ~ scale_x_log10(labels=function(x) format(x, big.mark = "'", scientific = FALSE)),
                               likelihood == "Bernoulli" ~ scale_x_log10(labels=function(x) format(x, big.mark = "'", scientific = FALSE))))
