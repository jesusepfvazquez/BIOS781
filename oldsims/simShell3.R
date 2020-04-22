# Population Stratification Settings 
# 2 Groups: Cases and Controls

regression=function(y,x){
  x=as.matrix(x)
  n=length(y)
  p=dim(x)[2]
  xTx.inv=solve(t(x)%*%x)
  betahat=as.vector(xTx.inv%*%t(x)%*%y)
  return(betahat)
}

# k number of populations, 1 or 2, 
agreement = function(k,iter){
  
  n=1000 # number of observations 
  p=10000 # number of snips
  m=p*0.01 # number of causal SNIPS
  errsd=0 # Variance of the error term, 0 if deterministic 
  n_test=100 
  thold = c(0.1, 0.15, 0.20, 0.5, 0.8, 1) # thresholds
  corr=matrix(ncol = length(thold),nrow=iter) #empty matrix vector 
  
  for(i in 1:iter){
    set.seed(i)
    
    # For one population
    if(k==1){
      set.seed(i)
      maf=runif(p,.05,.45)
      train=matrix(rbinom(p*n, 1, rep(maf,n)), ncol=p,byrow=TRUE)
      test=matrix(rbinom(p*n_test, 1, rep(maf,n_test)),ncol=p,byrow=TRUE)
    }
    
    #For 2 Populations
    ##Balding-Nichols model
    
    if (k == 2){
      MAF.hist=runif(p,0.1, 0.9)
      Fst= 0.001 #Default, generational effect
      a = MAF.hist*(1-Fst)/Fst
      b = (1-MAF.hist)*(1-Fst)/Fst
      MAF1=MAF2=MAF.hist
      for(z in 1:p){
        MAF1[z] = rbeta(1, a[z], b[z])
        MAF2[z] = rbeta(1, a[z],b[z])
      }
      ##which has mean p and variance Fst*p*(1-p)
      MAF1 = MAF1*(MAF1<(1-MAF1))+ (1-MAF1)*(MAF1>(1-MAF1))
      MAF2 = MAF2*(MAF2<(1-MAF2))+ (1-MAF2)*(MAF2>(1-MAF2))
      
      # selecting MAF for cases from pop1 and pop2
      which.pop.cases=rbinom(n/2,1,.55)+1
      MAF.cases=MAF1*(which.pop.cases==1)+MAF2*(which.pop.cases==2)
      
      # selecting MAF for controls from pop1 and pop2 
      which.pop.controls=rbinom(n/2,1,.45)+1
      MAF.controls=MAF1*(which.pop.controls==1)+MAF2*(which.pop.controls==2)
      
      #generating geno for controls and cases for training set
      geno.cases=matrix(rbinom(p*n/2,2,rep(MAF.cases,each=n/2)),nrow=n/2,ncol=p, byrow=T)
      geno.controls=matrix(rbinom(p*n/2,2,rep(MAF.controls,each=n/2)),nrow=n/2,ncol=p,byrow=T)
      train=rbind(geno.cases, geno.controls)
      
      #generating geno for controls and cases for testingset
      geno.cases=matrix(rbinom(p*n_test,2,rep(MAF.cases,each=n_test)),nrow=n_test,ncol=p, byrow=T)
      geno.controls=matrix(rbinom(p*n_test,2,rep(MAF.controls,each=n_test)),nrow=n_test,ncol=p,byrow=T)
      test=rbind(geno.cases, geno.controls)
    }
    
    
    beta1_original=rnorm(m,0,1) # generating the beta estimates
    beta1_pvalues = 1 - pnorm(abs(beta1_original)) + pnorm(-abs(beta1_original)) #get p-value from beta estimate. z = (beta-0)/1 = beta
    causalSNPS=sample(1:p,m)#choose causal snps
    beta1=rnorm(m,0,1)
    beta=rep(0,p)
    beta[causalSNPS]=beta1
    
    counter = 1
    for (t in thold) {
      beta_mod = beta[beta1_pvalues <= t] # apply treshhold 
      y=train[,beta1_pvalues <= t]%*%beta_mod+rnorm(dim(train)[1],0,errsd)
      y_test=test[,beta1_pvalues <= t]%*%beta_mod+rnorm(dim(test)[1],0,errsd)
      betahat=apply(train, 2, regression,y=y)
      yhat=test%*%betahat
      corr[i,counter]=cor(yhat,y_test)
      counter = counter + 1
    }
  
    print(i)  
  }
  
  return(corr)
}

# Results for 10 iterations 

## Two population results
corr_data = as.data.frame(agreement(2,10))
colnames(corr_data) = thold
means = colMeans(corr_data)
plot(thold,means, xlab = "P-value Treshold", ylab = "Correlation", main = "Correlation Level by Treshold Level,Two Population")

## One Population Results
corr_data_1 = as.data.frame(agreement(1,10))
colnames(corr_data_1) = thold
means_1 = colMeans(corr_data_1)
plot(thold,means_1, xlab = "P-value Treshold", ylab = "Correlation", main = "Correlation Level by Treshold Level, One Population")




