################################################################################
# Evaluate accuracy of predictive variances when using Algorithm 1.
################################################################################
library(gpboost)

###Generate data################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/gen_data.R")

set.seed(1)

m1 <- 5000
n <- 2*(2*m1)*10
sigma2 <- 0.5^2
sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
true_covpars <- c(sigma2,sigma2_1,sigma2_2)

the_data <- make_data(n=n, 
                      m1=m1,
                      sigma2=sigma2,
                      sigma2_1=sigma2_1,
                      sigma2_2=sigma2_2,
                      randef="Two_randomly_crossed_random_effects", 
                      likelihood="gaussian", 
                      has_F=F)

###Sample train and test set####################################################
set.seed(123)
i_train <- sample(n, n/2, replace = FALSE)
i_test <- 1:n
i_test <- i_test[!(i_test %in% i_train)]

group_train <- the_data$group_data[i_train,]
y_train <- the_data$y[i_train]
group_test <- the_data$group_data[i_test,]
y_test <- the_data$y[i_test]

###Cholesky#####################################################################
CHOLmodel <- GPModel(group_data = group_train,
                     likelihood="gaussian",
                     matrix_inversion_method = "cholesky")

time_chol <- system.time(pred_chol <- predict(CHOLmodel,
                                              y=y_train,
                                              group_data_pred=group_test,
                                              cov_pars = true_covpars,
                                              predict_var = TRUE,
                                              predict_cov_mat = FALSE,
                                              predict_response = FALSE))[3]

###Predictive variances using simulation########################################
nSim <- c(50, 100, 200, 500, 1000, 2000, 3000, 4000)
n_rep <- 100

ITresults <- data.frame(matrix(nrow=n_rep*length(nSim), ncol = 3))
colnames(ITresults) <- c("nSim", "RMSE", "time")

i <- 1
for(s in 1:length(nSim)){
  for(r in 1:n_rep){
    ITmodel <- GPModel(group_data = group_train,
                       likelihood="gaussian",
                       matrix_inversion_method = "iterative")
    
    ITmodel$set_optim_params(params = list(trace=TRUE,
                                           seed_rand_vec_trace=i,
                                           cg_preconditioner_type="symmetric_successive_over_relaxation"))
    ITmodel$set_prediction_data(nsim_var_pred = nSim[s])
    
    time <- system.time(pred_iterative <- predict(ITmodel, 
                                                  y=y_train,
                                                  group_data_pred=group_test,
                                                  cov_pars = true_covpars,
                                                  predict_var = TRUE, 
                                                  predict_cov_mat = FALSE,
                                                  predict_response = FALSE))[3]
    ITresults$time[i] <- time
    ITresults$nSim[i] <- nSim[s]
    ITresults$RMSE[i] <- sqrt(mean((pred_iterative$var-pred_chol$var)^2))
    
    i <- i + 1
    gc()
  }
}

################################################################################
# Plotting
################################################################################
library(ggplot2)

ITresults$nSim <- as.factor(ITresults$nSim)
mean <- aggregate(cbind(RMSE, time) ~ nSim, data=ITresults, FUN = mean)
twice_sd_mean <- 2 * aggregate(RMSE ~ nSim, data=ITresults, FUN = sd)[,2]/sqrt(n_rep)
agg_ITresults <- cbind(mean, RMSE_twice_sd_mean = twice_sd_mean)
my_colors <- RColorBrewer::brewer.pal(6,"Set1")

ggplot(agg_ITresults, aes(x=time, y=RMSE, label=nSim)) + 
  geom_line(linewidth=1) + 
  geom_point(size=2.5) + 
  scale_y_log10(n.breaks = 8) +
  scale_x_log10(n.breaks = 10) +
  ylab("RMSE") + xlab("Time (s)") + 
  theme_bw() +
  theme(legend.position = "top") +
  geom_errorbar(aes(ymin=RMSE-RMSE_twice_sd_mean, ymax=RMSE+RMSE_twice_sd_mean), linewidth=1, width = 0.01) + 
  geom_text(nudge_y=0.1, nudge_x=0.01, size=3.5) +
  scale_color_manual(values = my_colors, name="") +
  annotate("text", label = paste0("Cholesky: ", round(time_chol, digits = 1), "s"), x = 10, y = min(agg_ITresults$RMSE)+0.001, size=4)
