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
agreement = function(k){
  
  n=1000 # number of observations 
  p=10000 # number of snips
  m=p*.01 # number of causal SNIPS
  errsd=0 # Variance of the error term, 0 if deterministic 
  n_test=100 # 
  iter=100
  corr=rep(0,iter)
  
  for(i in 1:iter){
    set.seed(i)
    
    # For one population
    if(k==1){
      set.seed(i)
      maf=runif(p,.05,.45) #minor allele frequency for each snp
      train=matrix(rbinom(p*n, 1, rep(maf,n)), ncol=p,byrow=TRUE) #creates genotypes for each snp for all observations for training set
      test=matrix(rbinom(p*n_test, 1, rep(maf,n_test)),ncol=p,byrow=TRUE)  #same as above but for testing set
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
    
    causalSNPS=sample(1:p,m) #choose causal snps
    beta1=rnorm(m,0,1) #generate true coefficient for m causal snps
    beta=rep(0,p) #create empty beta vector
    beta[causalSNPS]=beta1 #add betas for causal snps in correct locations in beta vector, non-causal snps are zero
    y = train%*%beta+rnorm(m,0,errsd) #generate training set true outcomes
    y_test=test%*%beta+rnorm(m,0,errsd) #generate testing set true outcomes
    
    betahat=apply(train, 2, regression,y=y) #perform slr for each causal snp (if we want to use ridge we would do that here instead)
    yhat=test%*%betahat  #predicted y's for testing set
    corr[i]=cor(yhat,y_test) #correlation between true testing set y's and predicted testing set y's
  }
  
  output = as.data.frame(cbind(corr, rep(k, iter)))
  return(output)
}

corr_data = rbind(agreement(1), agreement(2))

# Histogram overlaid with kernel density curve
library(ggplot2)
ggplot(corr_data, aes(x=corr, fill=as.factor(V2))) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity") +
  labs(xlab='Correlation', )

# summary stats
tapply(corr_data$corr, corr_data$V2, summary)


