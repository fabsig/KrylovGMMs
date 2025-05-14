################################################################################
# Preconditioner comparison (SSOR, ZIC) for Bernoulli likelihoods
################################################################################
library(gpboost)

###Generate data################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("./../data/simulated/gen_data.R")
set.seed(1)

sigma2_1 <- 0.5^2
sigma2_2 <- 0.5^2
true_covpars <- c(sigma2_1,sigma2_2)

m1 <- 50000
n <- (2*m1) * 10

the_data <- make_data(n=n, 
                      m1=m1,
                      sigma2_1=sigma2_1,
                      sigma2_2=sigma2_2,
                      randef="Two_randomly_crossed_random_effects", 
                      likelihood="bernoulli_logit", 
                      has_F=F,
                      factor_m2=1)

####Krylov######################################################################
NUM_RAND_VEC_TRACE <- c(10, 20, 50, 100)
PRECONDITIONER <- c("symmetric_successive_over_relaxation", "zero_infill_incomplete_cholesky", "diagonal")
n_rep <- 100

CHOLresult <- NA
CHOLtime <- NA

Itresults <- data.frame(matrix(nrow=length(NUM_RAND_VEC_TRACE)*length(PRECONDITIONER)*n_rep,ncol = 4))
colnames(Itresults) <- c("preconditioner", "t", "negLL", "time")

i = 1
for(t in 1:length(NUM_RAND_VEC_TRACE)){
  for(p in 1:length(PRECONDITIONER)){
    for(r in 1:n_rep){
      Itmodel <- GPModel(group_data = the_data$group_data,
                         likelihood="bernoulli_logit",
                         matrix_inversion_method = "iterative")
      
      Itmodel$set_optim_params(params = list(init_cov_pars=true_covpars,
                                             num_rand_vec_trace=NUM_RAND_VEC_TRACE[t],
                                             cg_preconditioner_type = PRECONDITIONER[p],
                                             seed_rand_vec_trace=i))

      Itresults$preconditioner[i] <- PRECONDITIONER[p]
      Itresults$t[i] <- NUM_RAND_VEC_TRACE[t]
      Itresults$time[i] <- system.time(Itresults$negLL[i] <- Itmodel$neg_log_likelihood(cov_pars=true_covpars, y=the_data$y))[3]
      
      i = i+1
      gc()
    }
  }
}

###Cholesky#####################################################################
CHOLmodel <- GPModel(group_data = the_data$group_data,
                     likelihood="bernoulli_logit",
                     matrix_inversion_method = "cholesky")

CHOLmodel$set_optim_params(params = list(init_cov_pars=true_covpars))

CHOLtime <- system.time(CHOLresult <- CHOLmodel$neg_log_likelihood(cov_pars=true_covpars, y=the_data$y))[3]

################################################################################
#Plotting
################################################################################
library(ggplot2)
library(grid)

Itresults$preconditioner <- as.factor(Itresults$preconditioner)
Itresults$t <- as.factor(Itresults$t)

p1 <- ggplot(Itresults, aes(x=t, y=negLL, fill=preconditioner)) + 
      geom_hline(yintercept=CHOLresult, linetype = "dashed") + 
      theme_bw() +
      geom_boxplot() + labs(fill  = "") + ylab("log-likelihood") +
      scale_fill_brewer(type = "qual", palette=6, labels = scales::parse_format()) + 
      guides(fill = guide_legend(nrow = 1, byrow = TRUE)) + 
      theme(axis.title.x=element_blank(), 
            axis.text.x=element_blank(), 
            axis.ticks.x=element_blank(), 
            legend.position = "top", 
            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
  
p2 <- ggplot(Itresults, aes(x=t, y=time, color=preconditioner, shape=preconditioner)) +
      stat_summary(aes(group = preconditioner), fun = mean, geom = 'line', size=1, alpha=0.9) + 
      stat_summary(aes(group = preconditioner), fun = mean, geom = 'point', size=2) +
      scale_color_brewer(type = "qual", palette=6) + labs(color = "") + ylab("Time (s)") +
      theme_bw() + theme(legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
      annotate("text", label = paste0("Cholesky: ", round(CHOLtime, digits = 1), "s"), x = 4, y = 1.5, size=4)

grid.newpage()
grid.draw(rbind(ggplotGrob(p1), 
                ggplotGrob(p2), size = "first"))
