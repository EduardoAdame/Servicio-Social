##################################################
### PROBIT
##################################################

### Packages:
library(rjags)
library(mnormt)
library(MCMCpack) ### MCMC
library(LearnBayes) ### Bayes
library(coda)
library(boa)

setwd("~/Documents/_NotasCiencias/LibroBayesiana/Ejemplos/")
getwd()


##################################################
##################################################
### Simular datos
n = N = 100   # numero de sujetos
K = 3   # numero de parametros de regresion
betas = as.vector(c(2,2,-2))

set.seed(12345)
x1 = runif(n)
x2 = runif(n)
X = cbind(x1,x2,1)
eta = betas[3] + betas[1]*x1 + betas[2]*x2

p = pnorm(eta)
#u = runif(n)
#y = ifelse(p>u,1,0)
y = rep(NA,n)
for(i in 1:n){
  y[i] = rbinom(1,size=1,prob=p[i])
}
Y = as.vector(y)
plot(eta,y)
points(eta,p,col="blue")


##################################################
### Ajuste Frecuentista

mod <- glm(y ~ x1+x2, family=binomial(link="probit"))
plot(mod)
summary(mod)


##################################################
### Usando paquete de R: LearnBayes

modLearn <- bayes.probit(y,cbind(x1,x2,1),1000,prior=list(beta=rep(0,3),P=diag(1/1000,3)))
summary(modLearn)
plot(as.mcmc(modLearn$beta))
summary(modLearn$beta)

##################################################
### Usando paquete de R: MCMCpack

modMCMC <- MCMCprobit(y~x1+x2, burnin=1000, mcmc=1000, thin=1) 
summary(modMCMC)
plot(modMCMC)


##################################################
### JAGS

data <- list(
	y = y ,
	x1 = x1 ,
	x2 = x2 ,
	n = n  
)

param <- c("Beta")


### Probit

inits <- function(){	list(
  "Beta" = rnorm(3,0,0.1) 
)	}

fit <- jags.model("Bayes9_2Probit.bug", data, inits,  n.chains=3)

#update(fit,2000)

sample <- coda.samples(fit, param, n.iter=20000, thin=1)

plot(sample)
summary(sample)

par(mfrow=c(2,2))
traceplot(sample)


dic.samples(fit, n.iter=2000,thin=1, type="pD")

dic.samples(fit, n.iter=2000,thin=1, type="popt")




### Probit Augmented Data

initsA <- function(){	list(
  "Beta" = rnorm(3,0,0.1) , 
  "z" = rep(0,n) 
)	}

fitA <- jags.model("Bayes9_2ProbitAug.bug", data, initsA,  n.chains=4)

update(fitA,2000)

sampleA <- coda.samples(fitA, param, n.iter=20000, thin=1)

plot(sampleA)
summary(sampleA)

par(mfrow=c(2,2))
traceplot(sampleA)



##################################################
##################################################

##################################################
### MCMC 
### Gibbs Sampling

### Truncated Normal I[Z < trb]
normalderecha <- function(trb,mu,sig){
	rp <- pnorm(trb, mean=mu, sd=sig)
	u <- rp*runif(1)
	qnorm(u, mean=mu, sd=sig)
}

### Truncated Normal I[tra < Z]
normalizquierda <- function(tra,mu,sig){
	rp <- pnorm(tra, mean=mu, sd=sig)
	u <- rp+(1-rp)*runif(1)
	qnorm(u, mean=mu, sd=sig)
}

SIM = 1000

BETA21 = matrix(0,SIM,K)
ZZ21 = matrix(0,SIM,N)   # variables latentes

### Iniciales
Z = as.vector(rnorm(N,0,1))#matrix(0,N)
betas <- as.vector(rnorm(K,0,1))#solve(t(X)%*%X)%*%t(X)%*%Y   ### EMV

### Prior
### betas ~ Normal(be0m,be0v)
be0m = rep(0,K)### as.vector(c(2,2,-2))
be0v = diag(10000,K) ### diag(0.1,K)###

for(sim in 1:SIM){

### Z
	for(i in 1:N){
		if(Y[i]==1)Z[i] = normalizquierda(0,X[i,]%*%betas,1)
		else Z[i] = normalderecha(0,X[i,]%*%betas,1)
	}

### beta
	bm = solve(solve(be0v)+t(X)%*%X)%*%(solve(be0v)%*%be0m+t(X)%*%Z)
	bv = solve(solve(be0v)+t(X)%*%X)
	betas <- t(rmnorm(n=1, mean=bm, varcov=bv))

### sample
	BETA21[sim,] = betas
	ZZ21[sim,] = Z
}

plot(as.mcmc(BETA21))
summary(as.mcmc(BETA21))
pairs(BETA21)

plot(as.mcmc(ZZ21[,1:10]))

##################################################
##################################################

##################################################
### MCMC 
### Metropolis-Hasting

MHstepLogBeta <- function(betasold,betasnew,Y,X,bv){
  f1old = pmnorm(t(betasold), mean=rep(0,K),varcov=bv)
  f1new = pmnorm(t(betasnew), mean=rep(0,K),varcov=bv)
  f2old = pnorm(X%*%betasold)
  f2new = pnorm(X%*%betasnew)
  f3old = f2old*Y + (1-f2old)*(1-Y)
  f3new = f2new*Y + (1-f2new)*(1-Y)
  f1 = f1old[1]/f1new[1]
  f0 = ifelse(is.finite(f1),f1,1)
  f3 = f3new/f3old
  f4 = prod(ifelse(is.finite(f3),f3,1))
  ff = f4
  return(ff)
}



SIM = 1000

BETA31 = matrix(0,SIM,K)

### Iniciales
betas <- solve(t(X)%*%X)%*%t(X)%*%Y   ### EMV
countB = 0

### Prior
### betas ~ Normal(be0m,be0v)
be0m = rep(0,K)
be0v = diag(10000,K)

for(sim in 1:SIM){
  
  ### beta
  bv = diag(0.1,K)
  betasProposal <- t(rmnorm(n=1, mean=betas, varcov=bv))
  if(MHstepLogBeta(betas,betasProposal,Y,X,bv) >runif(1)){
    betas = betasProposal
    countB = countB+1
  }
  
  ### sample
  BETA31[sim,] = betas
}

summary(as.mcmc(BETA31))
plot(as.mcmc(BETA31))

##################################################
##################################################


