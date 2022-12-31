// Modelo de regresi�n log�stica con un predictor y una intersecci�n
data {
  int<lower=0> N; //N�mero de observaciones en el conjunto de datos de prueba 
  
  int<lower=0,upper=1> y [N]; // Variable dependiente como binaria
  
  vector <lower = 0> [N] x1; // Variable independiente 1
  
  vector <lower = 0> [N] x2; // Variable independiente 2
  
}


parameters { //Los par�metros del modelo no observados que queremos recuperar
  
  real alpha; //Intercepto en ordenadas
  
  real beta1;
  
  real beta2;
}


model {//Modelo a priori para los par�metros

  alpha ~ normal(0,100); // Intercepto
  
  beta1 ~ normal(0,100); // Beta 1
  
  beta2 ~ normal(0,100); // Beta 2
  
  y ~ bernoulli_logit(alpha + beta1 * x1 + beta2 * x2); //Modelo
}

generated quantities{//Simulacion de las cantidades de inter�s

vector [N] y_gorrito;

    for(i in 1:N){
      
      y_gorrito[i] = inv_logit(alpha + beta1 * x1[i] + beta2 * x2[i]); // modelo
    
    }
    
}


