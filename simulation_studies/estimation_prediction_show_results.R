################################################################################
# Create plots and tables for results from repeated parameter estimation and prediction 
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(xtable)

sigma2 <- 0.5^2
sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
true_beta <- c(0,rep(1,5))

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

Gaussian_results <- readRDS("./estimation_prediction_Gaussian.rds")
Bernoulli_results <- readRDS("./estimation_prediction_Bernoulli.rds")

Gaussian_results$likelihood <- "Gaussian"
Bernoulli_results$likelihood <- "Bernoulli"
Bernoulli_results$sigma2 <- NA

results <- rbind(Gaussian_results, Bernoulli_results)
results$method <- as.factor(results$method)
results$likelihood <- as.factor(results$likelihood)

###Plot variances###############################################################
c <- c("sigma2", "sigma2_1", "sigma2_2")
var_results_long <- gather(results[,c("method", "likelihood", c)], key="var_type", value="var_value", c, na.rm = T)

var_results_long$var_type <- factor(var_results_long$var_type, 
                                    levels=c("sigma2", "sigma2_1", "sigma2_2"),
                                    labels=c("sigma^2", "sigma[1]^2", "sigma[2]^2"))

ggplot(data=var_results_long, aes(x=method, y=var_value))+
  geom_hline(yintercept=sigma2, linetype="dashed") + 
  geom_boxplot() +
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE) +
  facet_nested_wrap(~likelihood + var_type, nrow=1, drop=T, scales = "free_x", labeller = label_parsed) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("") 

###Plot RMSE and log-score######################################################
c <- c("RMSE_latent_mu", "mean_log_score")
pred_results_long <- gather(results[,c("method", "likelihood", c)], key="pred_type", value="pred_value", c, na.rm = T)

pred_results_long$pred_type <- factor(pred_results_long$pred_type, 
                                    levels=c("RMSE_latent_mu", "mean_log_score"),
                                    labels=c("RMSE", "LS"))

ggplot(data=pred_results_long, aes(x=method, y=pred_value))+
  geom_boxplot() +
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE) +
  facet_nested_wrap(~likelihood + pred_type, nrow=1, drop=T, scales = "free", labeller = label_parsed) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("") 

###Table: Bias and RMSE of estimates############################################

#add true values
results$sigma2_true <- sigma2
results$sigma2_1_true <- sigma2_1
results$sigma2_2_true <- sigma2_2
results$beta_0_true <- true_beta[1]
results$beta_1_true <- true_beta[2]
results$beta_2_true <- true_beta[3]
results$beta_3_true <- true_beta[4]
results$beta_4_true <- true_beta[5]
results$beta_5_true <- true_beta[6]

rmse_bias <- results %>%
  group_by(across(all_of(c("likelihood", "method")))) %>%
  summarize(rmse_sigma2 = sqrt(mean((sigma2-sigma2_true)^2)),
            bias_sigma2 = mean(sigma2-sigma2_true),
            rmse_sigma2_1 = sqrt(mean((sigma2_1-sigma2_1_true)^2)),
            bias_sigma2_1 = mean(sigma2_1-sigma2_1_true),
            rmse_sigma2_2 = sqrt(mean((sigma2_2-sigma2_2_true)^2)),
            bias_sigma2_2 = mean(sigma2_2-sigma2_2_true),
            rmse_beta_0 = sqrt(mean((beta_0-beta_0_true)^2)),
            bias_beta_0 = mean(beta_0-beta_0_true),
            rmse_beta_1 = sqrt(mean((beta_1-beta_1_true)^2)),
            bias_beta_1 = mean(beta_1-beta_1_true),
            rmse_beta_2 = sqrt(mean((beta_2-beta_2_true)^2)),
            bias_beta_2 = mean(beta_2-beta_2_true),
            rmse_beta_3 = sqrt(mean((beta_3-beta_3_true)^2)),
            bias_beta_3 = mean(beta_3-beta_3_true),
            rmse_beta_4 = sqrt(mean((beta_4-beta_4_true)^2)),
            bias_beta_4 = mean(beta_4-beta_4_true),
            rmse_beta_5 = sqrt(mean((beta_5-beta_5_true)^2)),
            bias_beta_5 = mean(beta_5-beta_5_true))

print(xtable(rmse_bias))

###Table: Average RMSE and LS for predicton with standard errors################

average_sd_RMSE_LS <- results %>%
  group_by(across(all_of(c("likelihood", "method")))) %>%
  summarize(average_RMSE = mean(RMSE_latent_mu),
            sd_average_RMSE = sqrt(mean((RMSE_latent_mu-mean(RMSE_latent_mu))^2))/sqrt(n()),
            average_LS = mean(mean_log_score),
            sd_average_LS = sqrt(mean((mean_log_score-mean(mean_log_score))^2))/sqrt(n()))

print(xtable(average_sd_RMSE_LS))
