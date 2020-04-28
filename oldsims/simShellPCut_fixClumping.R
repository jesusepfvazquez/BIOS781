rm(list = ls()) # Clean the computing environment
library(bindata)
set.seed(12)

# n: number of obs
# m: number of causal snps
# correlated: 1 correlated snps or 0 for non-correlated snps
# size: size of block for correlated snips
# clumping: 1 cumpling 0 no cumpling
p = 10000 # number snps
#mlist = 100 # A list number causal snps
# errsd = 0 # error for generating true y values 0=deterministic model
n_test = 500 # number of observations in test dataset
iter = 100 # number of iterations 
plist = c(1, 0.1, 0.01, 10^(-3), 10^(-5), 10^(-8)) # a list of p-value cutoffs
r_thresh=.3

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
  corr_df = plist
  corr_mat = matrix(nrow = length(plist), ncol = iter)
  for (j in 1:length(plist)) {
    pmax = plist[j] # The p-cutoff
    for(i in 1:iter){
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
        blocks=cbind(locs,locs+round(runif(p/size,size-2,size-1)),runif(p/size,.6,.9),runif(p/size,.05,.15)) # define correlation within blocks
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
      
      train=apply(train, 2, function(x) (x-mean(x))/sd(x))
      test=apply(test, 2, function(x) (x-mean(x))/sd(x))
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
        temp = betahat
        # need to retrieve location of min p-value per block
        corr=cor(test, method="spearman")
        clumping_set=0
        mini=which.min(temp[2,])
        k=min(temp[2,])
        if(k<pmax){
          clumping_set[1]=mini
        }
        temp=temp[,-mini]
        while(k<pmax & dim(temp)[2]>0){
           k=min(temp[2,])
           mini=which.min(temp[2,])
           if(k<pmax){
             mcorr=max(corr[mini,clumping_set])
             if(mcorr<r_thresh){
               clumping_set=append(clumping_set,mini)
             }
            
           }
           temp=as.matrix(temp[,-mini],nrow=2)
           
        } 
        
        # Recalculating Analysis using only selected SNPS from clumping, now with p-value tresholding
        betahat[1,-c(clumping_set)]=0
        yhat = test %*% betahat[1,] 
        corr_mat[j, i] =  cor(yhat, y_test)
      }
    }
  }
  # combine the correlations to the p-value cutoffs
  corr_df = cbind(corr_df, corr_mat)
  # write the csv file
  write.csv(corr_df, paste0("\pine\scr\a\y\ayoung31\BIOS781\\", "m = ", m,"_", "n = ", n, "_", "p = ", p,
                            " correlated = ", correlated, " size = ", size,
                            "_", "clumping = ",clumping,  "corr.csv"))
}

PRS_corr(500,100,1,5,1)
