################################################################################
# Runtime analysis for parameter estimation for an unbalanced random effects design
# (m_2=m_1*0.5) and Bernoulli likelihoods.
################################################################################

library(gpboost)
library(lme4)
library(glmmTMB)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/gen_data.R")
set.seed(1)

M <- c(999, 1998, 3000, 4998, 9999, 19998, 49998, 99999, 199998, 499998, 999999)

sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
true_covpars <- c(sigma2_1,sigma2_2)
init_cov_pars <- c(0.5,0.5)
num_covariates <- 5
init_betas <- c(0.0062, rep(0, num_covariates))

res_cols <- c("m","method", "time_estimation")
results <- data.frame(matrix(nrow=length(M)*4, ncol = length(res_cols)))
colnames(results) <- res_cols

i <- 1
for(m in 1:length(M)){
  ###Generate data##############################################################
  set.seed(i)
  the_data <- make_data(n=M[m]*10, 
                        m1=M[m]/3*2,
                        factor_m2=0.5,
                        sigma2=sigma2_1, #for signal-to-noise ratio
                        sigma2_1=sigma2_1,
                        sigma2_2=sigma2_2,
                        randef="Two_randomly_crossed_random_effects", 
                        likelihood="bernoulli_logit", 
                        has_F=TRUE,
                        num_covariates=num_covariates)

  ###Cholesky###################################################################
  results$m[i] <- M[m]
  results$method[i] <- "Cholesky (GPBoost)"

  if(M[m] <= 19998){
    ##Estimation
    chol_model <- GPModel(group_data = the_data$group_data,
                          likelihood="bernoulli_logit",
                          matrix_inversion_method = "cholesky")

    chol_model$set_optim_params(params = list(maxit=1000,
                                              init_cov_pars=init_cov_pars,
                                              init_coef=init_betas))

    results$time_estimation[i] <- system.time(chol_model$fit(y=the_data$y, X=the_data$X))[3]
  }
  
  i <- i + 1
  
  ###Krylov#####################################################################
  results$m[i] <- M[m]
  results$method[i] <- "Krylov (GPBoost)"
  
  it_model <- GPModel(group_data = the_data$group_data,
                      likelihood="bernoulli_logit",
                      matrix_inversion_method = "iterative")
  
  it_model$set_optim_params(params = list(maxit=1000,
                                         init_cov_pars=init_cov_pars,
                                         init_coef=init_betas,
                                         seed_rand_vec_trace=i,
                                         cg_preconditioner_type="symmetric_successive_over_relaxation"))
  
  results$time_estimation[i] <- system.time(it_model$fit(y=the_data$y, X=the_data$X))[3]
  
  i <- i + 1
  
  ###lme4#######################################################################
  results$m[i] <- M[m]
  results$method[i] <- "lme4"

  if(M[m] <= 3000){
    formula = as.formula(paste0("y ~ -1 + ",paste0(colnames(the_data$X), collapse = ' + ')," + ",
                                paste0("(1|",colnames(the_data$group_data),")", collapse = ' + ')))
    results$time_estimation[i] <- system.time(lme4_model <- glmer(formula, family=binomial, data = data.frame(y = the_data$y, cbind(the_data$X, the_data$group_data)),
                                                                  start=list(theta = init_cov_pars, fixef = init_betas)))[3]
  }
    
  i <- i + 1
  
  ###glmmTMB####################################################################
  results$m[i] <- M[m]
  results$method[i] <- "glmmTMB"

  if(M[m] <= 9999){
    formula = as.formula(paste0("y ~ -1 + ",paste0(colnames(the_data$X), collapse = ' + ')," + ",paste0("(1|",colnames(the_data$group_data),")", collapse = ' + ')))
    results$time_estimation[i] <- system.time(glmmTMB_model <- glmmTMB(formula, family=binomial, data=data.frame(y = the_data$y, cbind(the_data$X, the_data$group_data)),
                                                                       start=list(theta = init_cov_pars, beta = init_betas)))[3]
  }
    
  i <- i + 1
  
  saveRDS(results, "./time_vs_m_unbalanced_Bernoulli.rds")
  gc()
}
