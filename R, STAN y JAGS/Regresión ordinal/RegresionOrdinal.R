##################################################
# Ordinal regression
##################################################
library(dplyr)
library(rstan)
library(rjags)
library(ggplot2)
options(mc.cores = parallel::detectCores())

##################################################
### generate data
Kappa = c(-Inf,0,2,4,Inf)
J = 4 ### categories 
set.seed(1234) 
n <- 200 # sujetos 

theta0 <- 1
theta1 <- 3
theta2 <- -1
X1 <- rnorm(n, 0, 1)
X2 <- rnorm(n, 0, 1)
eta <- theta0 + theta1*X1 + theta2*X2
y <- eta+rnorm(n,0,1)
outcome <- rep(NA,n) 
for(i in 1:n){  for(j in 1:J){
      if(Kappa[j]< y[i] & y[i]<=Kappa[j+1]){
        outcome[i] = j
} } } 
table(outcome)


data1 <- cbind(X1,X2,outcome,eta,y)
data1 <- data.frame(data1)
str(data1)
data2 <- list("X1"=X1,"X2"=X2,"outcome"=outcome-1,"n"=n)
data3 <- list("X1"=X1,"X2"=X2,"outcome"=outcome,"n"=n,"J"=J)

p1 <- ggplot(data1, aes(x=eta, y=y)) +
  geom_point() 
p1
p2 <- ggplot(data1, aes(x=eta, y=outcome)) +
  geom_point() 
p2

##################################################
setwd("~/Documents/_NotasCiencias/LibroBayesiana/Ejemplos/Bayes9_3Ordinal/")
##################################################

##################################################
### JAGS
##################################################
### Fit the model
param1 = c("theta0","theta1",
           "theta2","kappa"
) 

inits <- function(){	list(
  "theta0"=rnorm(1,0,0.001), 
  "theta1"=rnorm(1,0,0.001) ,  
  "theta2"=rnorm(1,0,0.001) , 
  "gam" = c(0:(J-2)*1) ,
  "Wout" = data2$outcome-0.5
)}

fit1 <- jags.model("RegresionOrdinal.bug", data2, inits,n.chains=2)

update(fit1,10000)
sample1 <- coda.samples(fit1, param1, n.iter=10000, thin=2)

plot(sample1)
summary(sample1)

##################################################
### STAN
##################################################
library(rstan)
library(MASS)

param2 = c("theta0","theta1","theta2",
           "kappa","kappa1")

fit2 <- stan("RegresionOrdinal.stan",
            data=data3,
            chains=2,warmup=1000,iter=2000,thin=1,cores=4) 

print(fit2, pars=param2)
plot(fit2, pars=param2)

stan_plot(fit2,pars=param2)
stan_plot(fit2,pars=param2, point_est = "mean", show_density = TRUE)
stan_trace(fit2,pars=param2)

##################################################
##################################################

