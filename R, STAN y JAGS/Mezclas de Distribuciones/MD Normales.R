
###MetodosAproximacionS
library(actuar)
library(fitdistrplus)
library(vcd)
library(ADGofTest)
library(mixtools)
library(mclust)
library(MASS)

### Famoso ejemplo de la literatura

###measurements give time in minutes between eruptions of the Old Faithful geyser in Yellowstone National Park

data(faithful)
attach(faithful)

faithful

hist(waiting, main = "Tiempo entre las erupciones del Old Faithful", xlab = "Minutos", ylab = "", cex.main = 1.5, 
     cex.lab = 1.5, cex.axis = 1.4,col="darkgreen",border="white",prob=T)
lines(density(waiting),col="red",lwd=2)

###Algoritmo "a mano"

### Valores para iniciar el algoritmo

s = c(0.5,50,80,5,5)

EM = function(x,s) {
  pi.hat= s[1]*dnorm(x, s[2], sqrt(s[4]))/(s[1]*dnorm(x, s[2], sqrt(s[4])) +(1-s[1])*dnorm(x, s[3], sqrt(s[5])))
  s[1] = mean(pi.hat)
  s[2] = sum(pi.hat*x) / sum(pi.hat)
  s[3] = sum((1-pi.hat)*x) / sum(1-pi.hat)
  s[4] = sqrt(sum(pi.hat*(x-s[2])^2) / sum(pi.hat))
  s[5] = sqrt(sum((1-pi.hat)*(x-s[3])^2) / sum((1-pi.hat)))
  s
}

iter = function(x, s) {
  s1 = EM(x,s)
  for (i in 1:5) {
    if (abs(s[i]-s1[i]) > 0.00001) {
      s=s1
      iter(x,s)
    }
    else s1
  }
  s1
}

x<-waiting

iter(x,s)

###Probabilidades de clasificación aposteriori

pi.hat<-function(s,x){s[1]*dnorm(x, s[2], sqrt(s[4]))/(s[1]*dnorm(x, s[2], sqrt(s[4])) +(1-s[1])*dnorm(x, s[3], sqrt(s[5])))}

prob.apos<-pi.hat(c(0.3510847, 54.2247274, 79.9173437,  5.4662678^2,  5.9875619^2),waiting)

head(prob.apos,20)

options(scipen=999)
cbind(waiting[1:20],head(prob.apos,20),1-head(prob.apos,20))

f<-function(x,s){s[1]*dnorm(x, s[2], sqrt(s[4])) +(1-s[1])*dnorm(x, s[3], sqrt(s[5]))}

x<-seq(min(waiting),max(waiting),0.01)
plot(x,f(x,c(0.3510847, 54.2247274, 79.9173437,  5.4662678^2,  5.9875619^2)),col="darkmagenta",ylab="",xlab="",lwd=2)
lines(density(waiting),col="darkblue",main="Mezcla de distribuciones normales",col.main="darkviolet",ylab="",xlab="",lwd=2)


######################################################################################################################
###Con mixtools

mezajuste <-normalmixEM(waiting, lambda = 0.5, mu = c(55, 80), sigma = 5) #Se supone que las normales en la mezcla tienen varianzas iguales

plot(mezajuste, density = TRUE, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.8, main2 = "Ajuste mezcla: tiempo entre las erupciones del Old Faithful", 
     xlab2 = "Minutos")


mezajuste[c("lambda", "mu", "sigma")]


### Y ajusta bien?

x<-seq(min(waiting),max(waiting),length.out=272)

p<-mezajuste$lambda
m<-mezajuste$mu
s<-mezajuste$sigma


### Ajuste de densidad

fest<-(p[1]*dnorm(x,m[1],s[1])+p[2]*dnorm(x,m[2],s[2]))

hist(waiting,col="darkgreen",border="white",prob=T,ylim=c(0,0.045))

lines(x,fest,col="blue",lwd=2)

### Ajuste de distribución

edest<-(p[1]*pnorm(x,m[1],s[1])+p[2]*pnorm(x,m[2],s[2]))

plot(ecdf(waiting),verticals=TRUE,do.points=FALSE,col.hor="red", col.vert="bisque",ylab="")
lines(x,edest,col="darkblue")

###Prueba formal de bondad de ajuste

mix.norm<-function(x, mu, sigma, pmix) {
  pmix[1]*pnorm(x,mu[1],sigma[1])+pmix[2]*pnorm(x,mu[2],sigma[2])
}

ks.test(waiting,mix.norm,mu=m,sigma=s,pmix=p)
ad.test(waiting,mix.norm,mu=m,sigma=s,pmix=p)