library(bindata)

p=100
n=10000
size=5
locs=seq(1,p,size)
blocks=cbind(locs,locs+round(runif(p/size,size-2,size-1)),runif(p/size,.8,.95))

sig=matrix(0,nrow=p,ncol=p)
for(i in 1:as.integer(p/size)){
  sig[blocks[i,1]:blocks[i,2],blocks[i,1]:blocks[i,2]]=blocks[i,3]
}
diag(sig)=1

maf=runif(p,.05,.45)
mysnips = rmvbin(n, margprob = maf, sigma = sig) + rmvbin(n, margprob = maf, sigma = sig)
cor(mysnips)[1:5,1:5]
sig[1:5,1:5]

cbind(as.matrix(colMeans(mysnips))/2, as.matrix(maf))
