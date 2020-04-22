n=1000 #number of observations in train dataset
p=10000 #number snps
m=100 #number causal snps
errsd=0 #error for generating true y values 0=deterministic model
n_test=100 #number of observations in test dataset
iter=100 #number of iterations 
corr=rep(0,iter) #initialize correlation matrix

# Population Stratification Settings 
# 2 Groups: Cases and Controls


regression=function(y,x){ #performs slr, need to add code to calculate p-values for thresholding
  x=as.matrix(x)
  n=length(y)
  p=dim(x)[2]
  xTx.inv=solve(t(x)%*%x)
  betahat=as.vector(xTx.inv%*%t(x)%*%y)
  return(betahat)
}

for(i in 1:iter){
  set.seed(i)
  maf=runif(p,.05,.45) #minor allele frequency for each snp
  train=matrix(rbinom(p*n, 1, rep(maf,n))+rbinom(p*n, 1, rep(maf,n)),ncol=p,byrow=TRUE) #creates genotypes for each snp for all observations for training set
  test=matrix(rbinom(p*n_test, 1, rep(maf,n_test))+rbinom(p*n_test, 1, rep(maf,n_test)),ncol=p,byrow=TRUE) #same as above but for testing set
  causalSNPS=sample(1:p,m)#choose which snps are causal snps
  beta1=rnorm(m,0,1)#generate true coefficient for m causal snps
  beta=rep(0,p)#create empty beta vector
  beta[causalSNPS]=beta1 #add betas for causal snps in correct locations in beta vector, non-causal snps are zero
  y=train%*%beta+rnorm(m,0,errsd) #generate training set true outcomes
  y_test=test%*%beta+rnorm(m,0,errsd) #generate testing set true outcomes
  
  betahat=apply(train, 2, regression,y=y) #perform slr for each causal snp (if we want to use ridge we would do that here instead)
  yhat=test%*%betahat #predicted y's for testing set
  corr[i]=cor(yhat,y_test) #correlation between true testing set y's and predicted testing set y's
}

boxplot(corr)#useful visualization

write.csv(corr,"descriptivename.csv")
