// Regresion lineal simple 

data {
  int<lower=0> N; //Número de datos (positivo)
  vector[N] x; //Covariables
  vector[N] y; //variables
}

parameters {
  real alpha; //Intercepto eje y (ordenada al orígen)
  real beta; // Pendiente
  real<lower=0> sigma; //Dispersin o desviacion estandar
}


model {
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  sigma ~ normal(0,1);
  
  y ~ normal(alpha + beta * x, sigma); //
}  

generated quantities{
    vector[N] y_sim; // Simulación de los datos de la posterior

    for(i in 1:N)
        y_sim[i] =normal_rng(alpha + beta * x[i], sigma);
}


