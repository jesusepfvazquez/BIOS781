rm(list = ls()) # Clean the computing environment
require(tidyverse)
library(bindata)
nmlist = list(c(100,25,1,5,1)) # A pair of (n, m, correlated, size, clumping)
# n: number of obs
# m: number of causal snps
# correlated: 1 correlated snps or 0 for non-correlated snps
# size: size of block for correlated snips
# clumping: 1 cumpling 0 no cumpling
p = 500 # number snps
#mlist = 100 # A list number causal snps
# errsd = 0 # error for generating true y values 0=deterministic model
n_test = 1000 # number of observations in test dataset
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
   pvalues = lm[["pval"]] * (lm[["pval"]] < pmax)
   return(cbind(betahat,pvalues))
}

p_thresh_beta_clumping <- function(y, x, pmax) {
  lm = regression(y = y, x = x)
  betahat = lm[["betahat"]]
  pvalues = lm[["pval"]]
  return(cbind(betahat,pvalues))
}

## The loop for calculating the predicted PRS and compare it with the true PRS 
PRS_corr <- function(n, m, correlated, size, clumping) {
  # initiate the final correlation data frame 
  corr_df = tibble(p_cuts = plist)
  corr_mat = matrix(nrow = length(plist), ncol = iter)
  for (j in 1:length(plist)) {
    pmax = plist[j] # The p-cutoff
    for(i in 1:iter){
      set.seed(i) 
      maf = runif(p,.05,.45) # minor allele frequency for each snp
      
      # Non-correlated SNIPS
      if (correlated == 0){
        train = matrix(rbinom(p*n, 1, rep(maf, n)) + rbinom(p*n, 1, rep(maf, n)), 
                       ncol = p, 
                       byrow = TRUE) # creates genotypes for each snp for all observations for training set
        test = matrix(rbinom(p*n_test, 1, rep(maf,n_test)) + rbinom(p*n_test, 1, rep(maf,n_test)), 
                      ncol = p, 
                      byrow = TRUE) # same as above but for testing set
      }
      
      # Correlated Snips
      if (correlated == 1){
        
        locs=seq(1,p,size) #defining blocks
        blocks=cbind(locs,locs+round(runif(p/size,size-2,size-1)),runif(p/size,.8,.95)) # define correlation within blocks
        sig=matrix(0,nrow=p,ncol=p) #varinace covariance matrix
        
        #create var-cov matrix: off-diagonals
        for(k in 1:as.integer(p/size)){
          sig[blocks[k,1]:blocks[k,2],blocks[k,1]:blocks[k,2]]=blocks[k,3]
        }
        
        #assing 1 to diagonal values
        diag(sig)=1
        # creates genotypes for each snp for all observations for training set
        train = rmvbin(n, margprob = maf, sigma = sig) + rmvbin(n, margprob = maf, sigma = sig) 
        # same as above but for testing set
        test = rmvbin(n_test, margprob = maf, sigma = sig) + rmvbin(n_test, margprob = maf, sigma = sig) 
        
      }
      
      causalSNPS = sample(1:p, m) # choose which snps are causal snps
      beta1 = rnorm(m, 0, 1) # generate true coefficient for m causal snps
      beta = rep(0, p)# create empty beta vector
      beta[causalSNPS] = beta1 # add betas for causal snps in correct locations in beta vector, non-causal snps are zero
      y = train %*% beta #+ rnorm(m, 0, errsd) # generate training set true outcomes
      y_test = test %*% beta #+ rnorm(m, 0, errsd) # generate testing set true outcomes

      if (clumping == 0){
        # perform linear regression for each snp to predict its effects with p-value thresholding
        betahat = apply(train, 2, p_thresh_beta, y = y, pmax = pmax)
        yhat = test %*% betahat[1,] # predicted y's for testing set
        corr_mat[j, i] =  cor(yhat, y_test) # correlation between true testing set y's and predicted y'
      }
      
      #clumping, selecting lowest p-value per bin before tresholding then redoing the analysis
      if (clumping == 1){
        # perform linear regression for each snp to predict its effects without p-value thresholding initially
        betahat = apply(train, 2, p_thresh_beta_clumping, y = y, pmax = pmax)
        # need to retrieve location of min p-value per block
        clumping_set = rep(0,dim(blocks)[1])
        for(tt in 1:dim(blocks)[1]){
          clumping_set[tt] <- which(betahat[2,blocks[tt,1]:blocks[tt,2]] == 
                                       min(betahat[2,blocks[tt,1]:blocks[tt,2]])) + blocks[tt,1] -1
        } 
        
        # Recalculating Analysis using only selected SNPS from clumping, now with p-value tresholding 
        y = train[,clumping_set] %*% beta[clumping_set] 
        y_test = test[,clumping_set] %*% beta[clumping_set] 
        betahat = apply(train[,clumping_set], 2, p_thresh_beta, y = y, pmax = pmax)
        yhat = test[,clumping_set] %*% betahat[1,] 
        corr_mat[j, i] =  cor(yhat, y_test)
      }

    }
    
    print(i)
  }
  # combine the correlations to the p-value cutoffs
  corr_df = cbind(corr_df, corr_mat)
  # write the csv file
  write_csv(corr_df, paste0("m = ", m,"_", "n = ", n, "_", "p = ", p,
                            " correlated = ", correlated, " size = ", size,
                            "_", "clumping = ",clumping,  "corr.csv"))
}

PRS_corr(100,25,1,10,1)
PRS_corr(100,25,0,0,0)




