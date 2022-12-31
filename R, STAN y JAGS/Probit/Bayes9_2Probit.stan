// Modelo de regresión logística con un predictor y una intersección
data {
  int<lower=0> N; //Número de defectos en el conjunto de datos de prueba 
  vector <lower = 0>[N] x; // Describe el comportamiento 
  int<lower=0,upper=1> y [N]; // Una variable que describe si se detectó cada defecto (1) o no (0)
  
  int <lower = 0> K; //Número de predicciones
  vector <lower=0> [K] k_pred; 
}

parameters { //Los parámetros del modelo no observados que queremos recuperar
  real alpha;
  real beta;
}

model {
  y ~ bernoulli_logit(alpha + beta * log(x));
  
  //Modelo a priori para los parámetros
  
  alpha ~ normal(0,1);
  beta ~ normal(1,1);
}

generated quantities{

vector [N] y_gorrito;

    for(i in 1:N){
      y_gorrito[i] = bernoulli(phi(alpha + beta * x[i]));
    }

