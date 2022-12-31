##################################################
# BAYESIAN
# TRUNCATED NORMAL DISTRIBUTIONS
##################################################
library(MASS)
library(ggplot2)
library(dplyr)
library(rjags)
library(rstan)
options(mc.cores = parallel::detectCores())

##################################################
### Truncated Normal I[tra < Z]
rnormleft <- function(tra,mu,sig){
  rp <- pnorm(tra, mean=mu, sd=sig)
  u <- rp+(1-rp)*runif(1)
  q <- qnorm(u, mean=mu, sd=sig)
  if(!is.finite(q)){ q = tra }
  return(q)
}
################################################## 
### generate data

set.seed(1234)

m1 <- 500 # sujetos en la clase I

outcome1 <- rep(NA,m1)
mu1 <- 5
sigma1 = 5 ### desviacion estandar
tau1 = 1/sigma1^2 ### precision
cota1 <- 0
for(j in 1:m1){		
  outcome1[j] <- rnormleft(cota1, mu1, sigma1)
}

data2 <- list(outcome=outcome1,n=m1,cota=cota1)
hist(outcome1,nclass=20,freq=FALSE)
lines(seq(mu1-3*sigma1,mu1+3*sigma1,length.out=101), dnorm(seq(mu1-3*sigma1,mu1+3*sigma1,length.out=101),mean=mu1,sd=sigma1),col="blue")

##################################################
### generate data using Walker and Damien (2001)
### Sampling Truncated Normal, Beta, and Gamma densitites
### Journal of Computational and Graphical Statistics, 10(2): 206-215

### Normal(mu,sigma^2)
V =  rgamma(m1,shape=3/2,scale=2)
Xa = mu1-sigma1*sqrt(V)
Xb = mu1+sigma1*sqrt(V) 
X = rep(NA,m1)
for(i in 1:m1){
  X[i] = runif(1,Xa[i],Xb[i]) # X|V  
}
hist(X,freq=FALSE,nclass=20)
lines(seq(mu1-3*sigma1,mu1+3*sigma1,length.out=101), dnorm(seq(mu1-3*sigma1,mu1+3*sigma1,length.out=101),mean=mu1,sd=sigma1),col="blue")


### Normal(mu,sigma^2)I[cota,infinity]
V =  rgamma(m1,shape=3/2,scale=2)
X = Xa = Xb = rep(NA,m1)
for(i in 1:m1){
  Xa[i] = max(mu1-sigma1*sqrt(V[i]),cota1)
  Xb[i] = mu1+sigma1*sqrt(V[i]) 
  X[i] = runif(1,Xa[i],Xb[i]) # X|V  
}
hist(X,freq=FALSE,nclass=20)
lines(seq(mu1-3*sigma1,mu1+3*sigma1,length.out=101), dnorm(seq(mu1-3*sigma1,mu1+3*sigma1,length.out=101),mean=mu1,sd=sigma1),col="blue")

##################################################
setwd("~/Documents/Articulos/Osvaldo.Ruth/HMMordinalLatentClass/")
##################################################

##################################################
### JAGS
##################################################
### Fit the model
param1 = c("mu","tau","sigma","pred")

fit1 <- jags.model("LCLMM0_normtrunc.bug", data2, n.chains=2)

update(fit1,1000)
sample1 <- coda.samples(fit1, param1, n.iter=10000, thin=2)


plot(sample1)
summary(sample1)
post1 <- summary(sample1)


### No funciona estimar con mixture of scale uniform distributions
inits3 <- list("tau"=1/1000)
fit3 <- jags.model("LCLMM0_normtrunc_scalemix.bug", data2, n.chains=2,inits=inits3)

update(fit3,50000)
sample3 <- coda.samples(fit3, param1, n.iter=10000, thin=2)

plot(sample3)
summary(sample3)

##################################################

##################################################
### STAN
##################################################

param2 = c("mu","sigma","pred")

fit2 <- stan("LCLMM0_normtrunc.stan",
            data=data2,
            chains=2,warmup=500,iter=1000,thin=1,cores=4)

print(fit2, pars=param2)
plot(fit2, pars=param2)

stan_plot(fit2,pars=param2)
stan_plot(fit2,pars=param2, point_est = "mean", show_density = TRUE)
stan_trace(fit2,pars=param2)




fit4 <- stan("LCLMM0_normtrunc_target.stan",
             data=data2,
             chains=2,warmup=1000,iter=2000,thin=1,cores=4)

print(fit4, pars=param2)
plot(fit4, pars=param2)

stan_plot(fit4,pars=param2)
stan_plot(fit4,pars=param2, point_est = "mean", show_density = TRUE)
stan_trace(fit4,pars=param2)





fit5 <- stan("LCLMM0_normtrunc_scalemix.stan",
             data=data2,
             chains=2,warmup=1000,iter=2000,thin=1,cores=4)

print(fit5, pars=param2)
plot(fit5, pars=param2)

stan_plot(fit5,pars=param2)
stan_plot(fit5,pars=param2, point_est = "mean", show_density = TRUE)
stan_trace(fit5,pars=param2)

##################################################

##################################################

