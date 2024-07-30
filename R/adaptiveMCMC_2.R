# Install and load necessary packages
install.packages("foreach")
install.packages("doParallel")
library(foreach)
library(doParallel)

adaptiveMCMC_2 <- function(dat, num.params, s2, num.iterations, num.burnin, rseed=1, info=0, prev.state=NA) {
  
  # Seed for random number generator
  set.seed(rseed)
  
  # dat is list, each element contains one or more X matrices with node information
  # save all samples (also from burnin iterations) for diagnostic purposes
  sample.params <- matrix(0, nrow=num.iterations, ncol=num.params)
  sample.postli <- rep(0, num.iterations)
  
  S1 <- matrix(0, nrow=num.params, ncol=1)
  S2 <- matrix(0, nrow=num.params, ncol=num.params)
  ns1s2 <- 0
  
  if (length(prev.state) == 1) {
    # Start from scratch
    k <- rep(1, num.params)
    ac <- rep(0, num.params) 
    cc <- rep(0, num.params)
    eig.sqrtval <- matrix(1, nrow=num.params, ncol=1)
    eig.vec <- diag(1, num.params, num.params)
    params <- matrix(0, nrow=num.params, ncol=1)
  } else {
    # Continue MCMC from given state
    k <- prev.state$k
    ac <- prev.state$ac
    cc <- prev.state$cc
    eig.sqrtval <- prev.state$eig.sqrtval
    eig.vec <- prev.state$eig.vec
    params <- prev.state$params
  }
  
  old.logli <- logprior(params, s2)
  cl <- makeCluster(detectCores() - 1) # Use all available cores except one
  registerDoParallel(cl)
  
  old.logli <- foreach(di = 1:length(dat), .combine = '+', .export = c('logprior', 'loglikelihood')) %dopar% {
    loglikelihood(dat[[di]], params)
  } + old.logli
  
  for (iter in 1:num.iterations) {
    if (info) {
      cat(".")
      if (iter %% 50 == 0) {
        cat(sprintf(" %d\n", iter))
      }
    }
    
    # Rotate SCAM 
    for (i in 1:num.params) {
      if (info > 1) {
        cat(sprintf("  i %d\n", i))
      }
      
      # Propose new parameter values
      new.params <- params + rnorm(1) * k[i] * eig.sqrtval[i] * eig.vec[, i, drop=F]
      
      # Likelihood
      new.logli <- logprior(new.params, s2)
      new.logli <- foreach(di = 1:length(dat), .combine = '+', .export = c('logprior', 'loglikelihood')) %dopar% {
        loglikelihood(dat[[di]], new.params)
      } + new.logli
      
      # Metropolis acceptance (priors already included in new and old logli)
      cc[i] <- cc[i] + 1
      
      if (info > 1) {
        print(new.params)
        cat(sprintf("iter %d new.logli %g old.logli %g\n", iter, new.logli, old.logli))
      }
      
      if (runif(1) < exp(new.logli - old.logli)) {
        ac[i] <- ac[i] + 1
        params <- new.params
        old.logli <- new.logli
      }
    } # End of rotate SCAM: all eigenvectors utilized
    
    