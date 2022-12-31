##################################################
### SIMULACIONES
##################################################

##################################################
### Packages:
library(rjags)
library(mnormt)
library(MCMCpack) ### MCMC
library(LearnBayes) ### Bayes

##################################################

setwd("/Users/EstadisticaCiencias/Documents/_NotasCiencias/LibroBayesiana/Ejemplos/")
getwd()

##################################################

##################################################
### Simular datos
n = 100
betas = as.vector(c(1,2,3))
sigma2 = 4
tau = 1/sigma2

set.seed(12345)
x1 = runif(n)
x2 = runif(n)
epsilon = rnorm(n,0,sqrt(sigma2))

eta = betas[3] + betas[1]*x1 + betas[2]*x2
y = eta + epsilon
plot(eta,y)
pairs(cbind(x1,x2,y))


##################################################
### Ajuste Frecuentista

mod <- lm(y ~ x1+x2)
plot(mod)
summary(mod)
anova(mod)

##################################################
##################################################
### Exacto: Familia Conjugada

### datos
Y = matrix(y)
X = cbind(x1,x2,1)

### Prior: hyperparameters 
b0 = matrix(0,3)
B0 = diag(1000,3)
a0 = 0.001
d0 = 0.001
shape0 = a0/2
scale0 = d0/2

### Posterior: hyperparameters 
b1 = solve( t(X)%*%X + solve(B0))%*%(t(X)%*%Y+solve(B0)%*%b0)
B1 = solve( t(X)%*%X + solve(B0))
a1 = n+a0
d1 = t(Y-X%*%b1)%*%(Y-X%*%b1) + t(b1-b0)%*%solve(B0)%*%(b1-b0) + d0


d1/(a1-2)   ### Estimador de \sigma^2
S2var = 2*(d1^2)/(((a1-2)^2)*(a1-4))   ###?Var(\sigma^2|Y)
b1   ### Estimador de \beta
Bvar = B1*d1[1,1]/(a1-2)   ### Var(\beta|Y)

##################################################
### Usando paquete de R: LearnBayes

modLearn <- blinreg(y,cbind(x1,x2,1),1000)
summary(modLearn)
plot(as.mcmc(modLearn$beta))
plot(as.mcmc(modLearn$sigma))
summary(modLearn$beta)
summary(modLearn$sigma)

##################################################
### Usando paquete de R: MCMCpack

modMCMC <- MCMCregress(y~x1+x2, burnin=1000, mcmc=1000, thin=1) 
summary(modMCMC)
plot(modMCMC)

##################################################
### Gibbs sampling

### datos
Y = matrix(y)
X = cbind(x1,x2,1)

###?hyperparameters prior
### \beta ~ Normal(b0,B0)
### \sigma2 ~ InversaGamma(shape0,scale0)
b0 = matrix(0,3)
B0 = diag(1000,3)
shape0 = 0.001
scale0 = 0.001

###?MCMC
SIM =1000
SIG2 = matrix(0,SIM)
BETA = matrix(0,SIM,3)

### iniciales (estimador maximo verosimil)
betas <- as.vector(rnorm(3,mean=0,sd=10))#solve(t(X)%*%X)%*%t(X)%*%Y
mu = X%*%betas
sig2 = rgamma(1,shape=1,scale=100)#t(Y-mu)%*%(Y-mu)/n

  BETA[1,] = betas
  SIG2[1] = sig2
  
for(sim in 1:(SIM-1)){
  
### beta
  Bf = solve(t(X)%*%X/c(sig2)+solve(B0))
  bf = Bf %*% (t(X)%*%Y/c(sig2)+solve(B0)%*%b0)
  betas <- t(rmnorm(n=1 , mean=bf, varcov =Bf))
  
### sigma2
  shapef = shape0+n/2
  scalef = scale0+sum((Y-X%*%betas)^2)/2
  sig2 = rinvgamma(1, shape =shapef, scale=scalef)
  
### MCMC sample
  BETA[sim+1,] = betas
  SIG2[sim+1] = sig2
}

plot(as.mcmc(BETA))
plot(as.mcmc(BETA[-c(1:5),]))
plot(as.mcmc(SIG2))
plot(as.mcmc(SIG2[-c(1:5)]))
summary(as.mcmc(BETA))
summary(as.mcmc(SIG2))


##################################################
### Metropolis-Hastings algorithm

### datos
Y = matrix(y)
X = cbind(x1,x2,1)

### hyperparameters prior
### \beta ~ Normal(b0,B0)
### \sigma ~ InversaGamma(shape0,scale0)
b0 = matrix(0,3)
B0 = diag(1000,3)
shape0 = 0.001
scale0 = 0.001

### MCMC
SIM =2000
SIG2 = matrix(0,SIM)
BETA = matrix(0,SIM,3)

### iniciales (estimadores maximo verosimiles)
betas <- solve(t(X)%*%X)%*%t(X)%*%Y
mu = X%*%betas
sig2 = t(Y-mu)%*%(Y-mu)/n
#betas <- as.vector(rnorm(3,mean=0,sd=10))
#sig2 = rgamma(1,shape=1,scale=100)

countB = 0
countS = 0

### Distribucion de Transicion
alphaMH <- function(betasOld,sig2Old,betasNew,sig2New, n,Y,X){
  ### verosimilitud
  like <- ((sig2New/sig2Old)^(-n/2)) * exp((-t(Y-X%*%betasNew)%*%(Y-X%*%betasNew)/(2*sig2New))+(t(Y-X%*%betasOld)%*%(Y-X%*%betasOld)/(2*sig2Old))) 
  ### iniciales
  prior <- 1
  ### Q transicion
  QB = 1
  QS <- 1
  ### posterior
  posterior <- like*prior*QB*QS
  return(posterior)
}


for(sim in 1:SIM){
  
### beta
  betasProposal <- t(rmnorm(n=1 , mean=betas, varcov=diag(0.1,3)))
  if(alphaMH(betas,sig2,betasProposal,sig2, n,Y,X) >runif(1)){
    betas = betasProposal
    countB = countB+1 
  }
### sigma2
  sig2Proposal = rnorm(1, sig2,1)
  if(alphaMH(betas,sig2,betas,sig2Proposal, n,Y,X) >runif(1)){
    sig2 = sig2Proposal
    countS = countS+1 
  }

### MCMC sample
  BETA[sim,] = betas
  SIG2[sim] = sig2
}

plot(as.mcmc(BETA))
plot(as.mcmc(SIG2))
countB
countS
summary(as.mcmc(BETA))
summary(as.mcmc(SIG2))


##################################################
### JAGS

datos <- list(
	y = y ,
	x1 = x1 ,
	x2 = x2 ,
	n = n , 
	x1.new = x1 ,
	x2.new = x2 ,
	nnew = n , 
	zeros = c(0,0,0) ,
	diagonal = diag(1/1000,3) ,
	a = 0.001 ,
	b = 0.001
)

param <- c("Beta","Sigma2")

inicial <- function(){	list(
  "Beta" = rnorm(3,0,1) , 
  "Tau" = rgamma(1,1,1)
)	}

fit <- jags.model(file="Bayes8_1LinearRegression.bug", data=datos, inits=inicial, n.chains=2,n.adapt=100)

update(fit,200)

sample <- coda.samples(m=fit, var=param, n.iter=100, thin=1)

plot(sample)
summary(sample)

par(mfrow=c(2,2))
traceplot(sample)


##################################################
##################################################

###Stan
data("iris")
versicolor = which(iris$Species == "versicolor")
x = iris$Sepal.Length[versicolor]
y = iris$Petal.Length[versicolor]
data = list(N = length(x), x = x , y = y)
require(rstan)
fit = stan(file = "Bayes8_1LinearRegression.stan",data = data)

summary(fit)
plot(x,y)
params = extract(fit)
alpha = mean(params$alpha)
beta = mean(params$beta)
abline(a = alpha, b = beta)
xr = seq(4,7.5,0.1)
yCI = sapply(xr ,function(k) quantile(params$beta * k + params$alpha , probs = c(0.05,0.95)))
lines(xr,yCI[1,], col = "red")
lines(xr, yCI[2,], col = "red")

fit = stan(file = "Bayes8_1LinearRegression.stan",data = data)

plot(density(y), xlim= c(2,6.5), ylim = c(0,1.4))
params=extract(fit)
for(i in 1:10){lines(density(params$y_sim[i,]),col = "red")}

y_new = params$y_sim[20,]
data_new = list(N=length(x), x=x, y= y_new)
fit_new = stan(file = "Bayes8_1LinearRegression.stan", data = data_new)
params_new = extract(fit_new)
plot(density(params$alpha))
lines(density(params_new$alpha),col = "red")
plot(density(params$beta))
lines(density(params_new$beta),col = "red")
####
