n=1000
p=10000
m=p*.01
errsd=0
n_test=100
iter=100
corr=rep(0,iter)

regression=function(y,x){
  x=as.matrix(x)
  n=length(y)
  p=dim(x)[2]
  xTx.inv=solve(t(x)%*%x)
  betahat=as.vector(xTx.inv%*%t(x)%*%y)
  return(betahat)
}

for(i in 1:iter){
  set.seed(i)
  maf=runif(p,.05,.45)
  train=matrix(rbinom(p*n, 1, rep(maf,n))+rbinom(p*n, 1, rep(maf,n)),ncol=p,byrow=TRUE)
  test=matrix(rbinom(p*n_test, 1, rep(maf,n_test))+rbinom(p*n_test, 1, rep(maf,n_test)),ncol=p,byrow=TRUE)
  causalSNPS=sample(1:p,m)#choose causal snps
  beta1=rnorm(m,0,1)
  beta=rep(0,p)
  beta[causalSNPS]=beta1
  y=train%*%beta+rnorm(m,0,errsd)
  y_test=test%*%beta+rnorm(m,0,errsd)
  
  betahat=apply(train, 2, regression,y=y)
  yhat=test%*%betahat
  corr[i]=cor(yhat,y_test)
}
write.csv(corr,"descriptivename.csv")