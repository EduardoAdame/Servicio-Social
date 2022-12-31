##################################################
### Ejemplo 7.3 (pag. 202)
### Robert & Casella (2010)
### Introducing Monte Carlo Methods with R
### Springer
##################################################

##################################################
### Packages:
library(rjags)
library(MCMCpack) ### MCMC

##################################################

setwd("/Users/EstadisticaCiencias/Documents/_NotasCiencias/LibroBayesiana/Ejemplos/")
getwd()

##################################################

##################################################
### Simular datos
n = 100
theta = 3   ### media
sigma2 = 1   ### varianza

set.seed(12345)
X = rnorm(n,theta,sqrt(sigma2))

hist(X,nclass=10)

##################################################
### Ajuste Frecuentista

mod <- lm(X ~ 1)
summary(mod)

mean(X)   ### Estimador de la Media
sum((X-mean(X))^2)/(n-1)   ###?Estimador de la varianza

##################################################

##################################################
### Gibbs sampling JAGS

### hyperparameters prior No informative 
theta0 = 0
tau2 = 100
a = 0.1
b = 0.1

### MCMC
SIM =1000   # chain MCMC
SIGMA2 = matrix(0,SIM)
THETA = matrix(0,SIM)

### iniciales
theta = mean(X)   ### Estimador de la Media
sigma2 = sum((X-mean(X))^2)/(n-1)   ### Estimador de la varianza


for(sim in 1:SIM){
  
### theta
  theta.var = sigma2*tau2/(sigma2+n*tau2)
  theta.mean = (sigma2*theta0 + n*tau2*mean(X))/(sigma2+n*tau2)
  theta <- rnorm(1,theta.mean,sqrt(theta.var))
  
### sigma2
  sigma2.shape = n/2+a
  sigma2.scale = sum((X-theta)^2)/2+b
  sigma2 = rinvgamma(1,sigma2.shape,sigma2.scale)
  
### MCMC sample
  THETA[sim,] = theta
  SIGMA2[sim] = sigma2
}

plot(as.mcmc(THETA))
plot(as.mcmc(SIGMA2))
summary(as.mcmc(THETA))
summary(as.mcmc(SIGMA2))

##################################################
### JAGS

data <- list(
	X = X,
	n = n , 
	theta0 = 0 ,
	tau2inv = 1/100 ,
	a = 0.1 ,
	b = 0.1 
)

param <- c("theta","sigma2")

inits <- function(){	list(   ### opcional
  "theta" = rnorm(1,0,1) , 
  "sigma2inv" = rgamma(1,1,1)
)	}

fit <- jags.model("Bayes7_3DistribucionNormal.bug", data, inits,n.chains=2) 

update(fit,1000)

sample <- coda.samples(fit, param, n.iter=100, thin=1)

plot(sample)
summary(sample)

traceplot(sample)


##################################################
##################################################
#### STAN

#Este es el paso de compliación del modelo
stan_norm_cpp <- stan_model("Bayes7_3DistribucionNormal.stan")

#El modelo ya esta especificado, pero aun falta indicar los datos observados.
N <- 20 # Observamos 20 realizaciones
set.seed(122) #Semila para las realizaciones
#Distribucion muestral
y <- rnorm(N, 0, 1) #Mecanismo aleatorio de los datos dado el estado de informacion previo 

#Ponemos los datos en una lista
data_list <- list(y = y, N = N )
#corremos las muestras
norm_fit <- sampling(object = stan_norm_cpp, data = data_list,
                     chains = 3 , iter = 1000 , warmup = 500)

#Podemos ver los intervalos de la distribucipon posterior de los parámetros 
norm_fit

#Podemos realizar graficas de las cadenas con más detalle (y con facilidad)
#usando el paquete bayesplot.

library(bayesplot)
#Modelo con las 3 cadenas de markov, 100 iteraciones y una warmup de 500
norm_posterior_inc_warmup <- extract(norm_fit, inc_warmup = TRUE, 
                                     permuted = FALSE)

color_scheme_set("mix-blue-red")
#Monte Carlo Markov Chains trazo
mcmc_trace(norm_posterior_inc_warmup,  pars = c("mu", "sigma2"), n_warmup = 500,
           facet_args = list(nrow = 2, labeller = label_parsed)) + 
  facet_text(size = 15)

#Y podemos graficar la distribución posterior de los parámetros, con intervalos de probabilidad.
norm_posterior <- as.array(norm_fit)

mcmc_areas(norm_posterior, pars = c("mu", "sigma2"), prob = 0.95, 
           point_est = "median", adjust = 1.4) 
