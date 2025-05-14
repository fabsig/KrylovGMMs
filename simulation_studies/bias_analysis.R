################################################################################
# Variance parameter estimation with a Laplace approximation for different 
# numbers of repeated observations per group level d=n/m_1.
################################################################################

library(gpboost)
library(lme4)
library(glmmTMB)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/gen_data.R")
set.seed(1)

m1 <- 2000
N <- c(20000, 40000, 100000, 200000, 400000)

n_rep <- 100

sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
true_covpars <- c(sigma2_1,sigma2_2)

res_cols <- c("n", "d", "sigma2_1", "sigma2_2")
results <- data.frame(matrix(nrow=n_rep*length(N), ncol = length(res_cols)))
colnames(results) <- res_cols

i <- 1
for(n in 1:length(N)){
  for(r in 1:n_rep){
    print(i)
    
    ###Generate data############################################################
    set.seed(i)
    the_data <- make_data(n=N[n], 
                          m1=m1,
                          sigma2=sigma2_1, #for signal-to-noise ratio
                          sigma2_1=sigma2_1,
                          sigma2_2=sigma2_2,
                          randef="Two_randomly_crossed_random_effects", 
                          likelihood="bernoulli_logit", 
                          has_F=TRUE,
                          num_covariates=5)

    ###Estimation###############################################################
    it_model <- GPModel(group_data = the_data$group_data,
                        likelihood="bernoulli_logit",
                        matrix_inversion_method = "iterative")
    
    it_model$set_optim_params(params = list(maxit=1000,
                                            trace=TRUE,
                                            seed_rand_vec_trace=r,
                                            cg_preconditioner_type="symmetric_successive_over_relaxation"))
    
    it_model$fit(y=the_data$y, X=the_data$X)
    results$sigma2_1[i] <- it_model$get_cov_pars()[1]
    results$sigma2_2[i] <- it_model$get_cov_pars()[2]
    results$n[i] <- N[n]
    results$d[i] <- N[n]/m1

    i <- i + 1
    gc()
  }
}

############################################################################
# Plotting
############################################################################
library(ggplot2)

results$d <- as.factor(results$d)

ggplot(results, aes(x=d, y=sigma2_1, group=d)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               linewidth=0.5,
               width = 0.8,
               fun.max = function(x) mean(x) + 2*sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - 2*sd(x) / sqrt(length(x))) + 
  stat_summary(fun=mean, colour="darkred", geom="point",shape=18, size=3, show.legend = FALSE)  +
  theme_bw() + ylab(expression(sigma[1]^2)) + geom_hline(yintercept=sigma2_1, linetype="dashed")
