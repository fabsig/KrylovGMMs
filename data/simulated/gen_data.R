
make_data <- function(n, m1, 
                      sigma2=0.5^2, sigma2_1=0.5^2, sigma2_2=0.5^2, sigma2_3=0.5^2,
                      randef, likelihood, has_F = FALSE, factor_m2=1, num_covariates=5){
  group <- rep(1,n)
  for(i in 1:m1) group[((i-1)*n/m1+1):(i*n/m1)] <- i
  b1 <- sqrt(sigma2_1) * rnorm(m1)
  if(randef == "One_random_effect"){
    eps <- b1[group]
    group_data <- group
  } else if(randef == "Two_completely_crossed_random_effects"){
    n_obs_gr <- n/m1
    group2 <- rep(1,n)
    for(i in 1:m1) group2[(1:n_obs_gr)+n_obs_gr*(i-1)] <- 1:n_obs_gr
  } else if(randef == "Two_randomly_crossed_random_effects"){
    if(factor_m2 != 1){
      m2 <- factor_m2*m1
      group2 <- rep(1,n)
      for(i in 1:m2) group2[((i-1)*n/m2+1):(i*n/m2)] <- i
      group2 <- group2[sample.int(n=n, size=n, replace=FALSE)]
    }
    else{
      group2 <-  group[sample.int(n=n, size=n, replace=FALSE)]  
    }
  } else if(randef == "Two_nested_random_effects"){
    m_nested <- m1*2
    group2 <- rep(1,n)
    for(i in 1:m_nested) group2[((i-1)*n/m_nested+1):(i*n/m_nested)] <- i
  }
  else if(randef == "Three_randomly_crossed_random_effects"){
    group2 <-  group[sample.int(n=n, size=n, replace=FALSE)]
    group3 <-  group[sample.int(n=n, size=n, replace=FALSE)]
  }
  if(randef != "One_random_effect"){
    if(randef == "Three_randomly_crossed_random_effects") {
      b2 <- sqrt(sigma2_2) * rnorm(length(unique(group2)))
      b3 <- sqrt(sigma2_3) * rnorm(length(unique(group3)))
      if (length(unique(group2)) != max(group2)) stop("not all levels samples -> gives index problems")
      if (length(unique(group3)) != max(group3)) stop("not all levels samples -> gives index problems")
      eps <- b1[group] + b2[group2] + b3[group3]
      group_data <- cbind(group,group2,group3)
    } else {
      b2 <- sqrt(sigma2_2) * rnorm(length(unique(group2)))
      if (length(unique(group2)) != max(group2)) stop("not all levels samples -> gives index problems")
      eps <- b1[group] + b2[group2]
      group_data <- cbind(group,group2)
    }
  }
  #Simulate fixed effects
  if (has_F) {
    #The covariates X are sampled from a normal distribution with mean 0 and variance chosen such that 
    #the signal-to-noise ratio between the fixed and random effects (sigma2) is one, 
    #and the true regression coefficients are all 1 except for the intercept which is 0.
    beta <- c(0,rep(1,num_covariates))
    X <- cbind(rep(1,n),matrix(rnorm(num_covariates*n),nrow=n))
    f <- as.vector(X%*%beta)
    delta_var_f <- sqrt(sigma2/var(f))
    X[,-1] <- X[,-1] * delta_var_f
    f <- as.vector(X%*%beta)
    colnames(X) <- c("1",paste0("Cov_",1:num_covariates))
  }else{
    X <- matrix(rep(0,(num_covariates+1)*n),ncol=num_covariates+1)
    colnames(X) <- c("1",paste0("Cov_",1:num_covariates))
    f <- 0
  }
  #Simulate response variable
  if (likelihood == "bernoulli_probit") {
    probs <- pnorm(f+eps)
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "bernoulli_logit") {
    probs <- 1/(1+exp(-(f+eps)))
    y <- as.numeric(runif(n) < probs)
  } else if (likelihood == "poisson") {
    mu <- exp(f+eps)
    y <- qpois(runif(n), lambda = mu)
  } else if (likelihood == "gamma") {
    mu <- exp(f+eps)
    y <- qgamma(runif(n), scale = mu, shape = 1)
  } else if (likelihood == "gaussian") {
    mu <- f + eps
    y <- sqrt(sigma2) * rnorm(n) + mu
  }
  list(y=y, X=X, group_data=group_data, eps=eps, f=f)
}
