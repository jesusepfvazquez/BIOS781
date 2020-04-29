rm(list = ls()) # Clean the computing environment
set.seed(1000)
require(tidyverse)

n = 500 
m = 10
p = 10000 # number snps 
# errsd = 0 # error for generating true y values 0=deterministic model
n_test = 100 # number of observations in test dataset
# plist = c(1, 0.1, 0.01, 10^(-3), 10^(-5), 10^(-8)) # a list of p-value cutoffs

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

getPval <- function(y, x){
  lm = regression(y = y, x = x)
  return(lm[["pval"]])
}

## The loop for calculating the predicted PRS and compare it with the true PRS 
PRS_detect <- function(n, m, iter) {
  set.seed(1000)
  plots = list()
  # initiate the final correlation data frame 
  for (i in 1:iter) {
    maf = runif(p,.05,.45) # minor allele frequency for each snp
    train = matrix(rbinom(p*n, 1, rep(maf, n)) + rbinom(p*n, 1, rep(maf, n)), 
                  ncol = p, 
                  byrow = TRUE) # creates genotypes for each snp for all observations for training set
    test = matrix(rbinom(p*n_test, 1, rep(maf,n_test)) + rbinom(p*n_test, 1, rep(maf,n_test)), 
                ncol = p, 
                byrow = TRUE) # same as above but for testing set
    train = apply(train, 2, function(x) {(x - mean(x))/sd(x)})
    test = apply(test, 2, function(x) {(x - mean(x))/sd(x)})
    causalSNPS = 1:m # choose which snps are causal snps
    beta1 = rnorm(m) # generate true coefficient for m causal snps
    beta1 = beta1[order(abs(beta1), decreasing = TRUE)]
    beta = rep(0, p)# create empty beta vector
    beta[causalSNPS] = beta1 # add betas for causal snps in correct locations in beta vector, non-causal snps are zero
    y = train %*% beta
    betahat_p = apply(train, 2, getPval, y = y)
    p_val_ranks = rank(betahat_p)
    df = tibble(beta_ID = 1:p,
                beta_hat_pval = betahat_p,
                p_val_ranks = p_val_ranks,
                neg_log_10_pval = -log10(betahat_p))[1:m, ]
      
    plots[[i]] = ggplot(data = df, aes(x = beta_ID, y = neg_log_10_pval)) +
      geom_point(size = 1.5) +
      geom_text(data = df, aes(label = signif(beta_hat_pval, 3)), 
                vjust = -1, hjust = 0.1) +
      geom_text(data = df, aes(label = paste0("(", p_val_ranks, ")")), 
                vjust = 1.25) +
      geom_hline(yintercept = 3, col = "red") + 
      geom_hline(yintercept = 8, col = "blue") +
      scale_x_continuous(breaks = 1:10) + 
      xlab("Beta ID") + 
      ylab(expression(paste("-log"["10"], "(p-value)"))) 
  }
  return(plots)
}

## The loop for calculating the predicted PRS and compare it with the true PRS 
PRS_detect_oneplot <- function(n, m, iter) {
  set.seed(1000)
  true_m = m
  m = 2*m
  df = tibble('0' = 1:m)
  # initiate the final correlation data frame 
  for (i in 1:iter) {
    maf = runif(p,.05,.45) # minor allele frequency for each snp
    train = matrix(rbinom(p*n, 1, rep(maf, n)) + rbinom(p*n, 1, rep(maf, n)), 
                   ncol = p, 
                   byrow = TRUE) # creates genotypes for each snp for all observations for training set
    test = matrix(rbinom(p*n_test, 1, rep(maf,n_test)) + rbinom(p*n_test, 1, rep(maf,n_test)), 
                  ncol = p, 
                  byrow = TRUE) # same as above but for testing set
    train = apply(train, 2, function(x) {(x - mean(x))/sd(x)})
    test = apply(test, 2, function(x) {(x - mean(x))/sd(x)})
    causalSNPS = 1:(m/2) # choose which snps are causal snps
    beta1 = rnorm(m/2) # generate true coefficient for m causal snps
    beta1 = beta1[order(abs(beta1), decreasing = TRUE)]
    beta = rep(0, p)# create empty beta vector
    beta[causalSNPS] = beta1 # add betas for causal snps in correct locations in beta vector, non-causal snps are zero
    y = train %*% beta
    betahat_p = apply(train, 2, getPval, y = y)[1:m]
    neg_log_10_p = -log10(betahat_p) 
    df = cbind(df, neg_log_10_p)
  }
  colnames(df)[2:dim(df)[2]] = 1:iter
  
  pplot = ggplot(data = df, aes_(x = as.name(0), y = as.name(1))) +
    geom_point(size = 1) +
    geom_hline(yintercept = 3, col = "red") +
    geom_hline(yintercept = 8, col = "blue") +
    geom_vline(xintercept = (true_m + 0.5), col = "grey") + 
    scale_x_continuous(breaks = 1:10) +
    xlab("Beta ID") +
    ylab(expression(paste("-log"["10"], "(p-value)"))) +
    theme_classic()
  
  for (i in 2:iter) {
    pplot = pplot +  geom_point(data = df,
                                aes_(x = as.name(0), y = as.name(i)),
                                size = 1,
                                position = 'jitter') 
  }
  
  return(pplot)

}

dot_plot_2 = PRS_detect_oneplot(n, m, 100)

dot_plot_3 = dot_plot_2 + 
  annotate("text", x = 5, y = 200, label = "Causal SNPs") +
  annotate("text", x = 15, y = 200, label = "Null SNPs")
  

# dot_plot = PRS_detect_oneplot(n, m, 100)
# dot_plot + title(main = "P-values of The Causal SNPs' Effects in The GWAS Model")
ggsave("dot_plot_3.jpg", width = 7, height = 7, units = "in")


# ps = PRS_detect(n = n, m = m)
# 
# 
# library(gridExtra)
# 
# png(filename = "detect_beta.png", width = 1280, height = 960, units = "px")
# do.call(grid.arrange, c(ps, list(ncol = 3, nrow = 2)))
# dev.off()

