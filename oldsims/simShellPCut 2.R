rm(list = ls()) # Clean the computing environment
require(tidyverse)

nmlist = list(c(5000, 10000), c(500, 10)) # A pair of (n, m)
# n: number of obs
# m: number of causal snps
p = 10000 # number snps
mlist = 100 # A list number causal snps
# errsd = 0 # error for generating true y values 0=deterministic model
n_test = 100 # number of observations in test dataset
iter = 100 # number of iterations 
plist = c(1, 0.1, 0.01, 10^(-3), 10^(-5), 10^(-8)) # a list of p-value cutoffs

## custom function performing linear regression & returning summary statistics
regression <- function(y, x){ 
  x = as.matrix(x) # Design Matrix
  n = length(y) 
  p = dim(x)[2] 
  xTx.inv = solve(t(x) %*% x)
  betahat = as.vector(xTx.inv %*% t(x) %*% y) # OLS solution for beta
  resid = y - x %*% betahat # residual
  varhat = (sum(resid^2))*diag(xTx.inv)/(n - p) # estimated variance for beta
  se = sqrt(varhat)
  tval = betahat/se # t statistic
  pval = 2*pt(-abs(tval), df = n - p)
  out = list(betahat, se, tval , n - p , pval) 
  names(out) = c("betahat", "se", "t", "df", "pval")
  return(out)
}

## custom function to perform p-value thresholding on beta
p_thresh_beta <- function(y, x, pmax) {
   lm = regression(y = y, x = x)
   betahat = lm[["betahat"]] * (lm[["pval"]] < pmax)
   return(betahat)
}

## The loop for calculating the predicted PRS and compare it with the true PRS 
PRS_corr <- function(n, m, correlated) {
  # initiate the final correlation data frame 
  corr_df = tibble(p_cuts = plist)
  corr_mat = matrix(nrow = length(plist), ncol = iter)
  for (j in 1:length(plist)) {
    pmax = plist[j] # The p-cutoff
    for(i in 1:iter){
      set.seed(i) 
      maf = runif(p,.05,.45) # minor allele frequency for each snp
      
      if (correlated == 0){
        train = matrix(rbinom(p*n, 1, rep(maf, n)) + rbinom(p*n, 1, rep(maf, n)), 
                       ncol = p, 
                       byrow = TRUE) # creates genotypes for each snp for all observations for training set
        test = matrix(rbinom(p*n_test, 1, rep(maf,n_test)) + rbinom(p*n_test, 1, rep(maf,n_test)), 
                      ncol = p, 
                      byrow = TRUE) # same as above but for testing set
      }
      
      if (correlated == 1){
        
        size=5
        locs=seq(1,p,size)
        blocks=cbind(locs,locs+round(runif(p/size,size-2,size-1)),runif(p/size,.8,.95))
        sig=matrix(0,nrow=p,ncol=p)
        
        for(i in 1:as.integer(p/size)){
          sig[blocks[i,1]:blocks[i,2],blocks[i,1]:blocks[i,2]]=blocks[i,3]
        }
    
        diag(sig)=1
        
        mysnips = rmvbin(n, margprob = maf, sigma = sig) + rmvbin(n, margprob = maf, sigma = sig)
        
        train = matrix(rbinom(p*n, 1, rep(maf, n)) + rbinom(p*n, 1, rep(maf, n)), 
                       ncol = p, 
                       byrow = TRUE) # creates genotypes for each snp for all observations for training set
        test = matrix(rbinom(p*n_test, 1, rep(maf,n_test)) + rbinom(p*n_test, 1, rep(maf,n_test)), 
                      ncol = p, 
                      byrow = TRUE) # same as above but for testing set
      }

      causalSNPS = sample(1:p, m) # choose which snps are causal snps
      beta1 = rnorm(m, 0, 1) # generate true coefficient for m causal snps
      beta = rep(0, p)# create empty beta vector
      beta[causalSNPS] = beta1 # add betas for causal snps in correct locations in beta vector, non-causal snps are zero
      y = train %*% beta #+ rnorm(m, 0, errsd) # generate training set true outcomes
      y_test = test %*% beta #+ rnorm(m, 0, errsd) # generate testing set true outcomes
    
      # perform linear regression for each snp to predict its effects with p-value thresholding
      betahat = apply(train, 2, p_thresh_beta, y = y, pmax = pmax) 
    
      yhat = test %*% betahat # predicted y's for testing set
      corr_mat[j, i] =  cor(yhat, y_test) # correlation between true testing set y's and predicted y'
    }
  }
  # combine the correlations to the p-value cutoffs
  corr_df = cbind(corr_df, corr_mat)
  # write the csv file
  write_csv(corr_df, paste0("m = ", m,"_", "n = ", n, "_", "p = ", p, "_", "corr.csv"))
}

## The loop for executing the function
for (i in nmlist) {
  n = i[1]
  m = i[2]
  PRS_corr(n, m)
}



