//
// Este programa de Stan define un modelo simple, con un 
// vector de valores 'X[N]' modelado como la distribuci�n normal
// con esperanza 'mu' y desviaci�n est�ndar 'sigma'.

// Los datos de entrada son un vector 'X[N]' de longitud 'N'.

data {
  int<lower = 0> N; // N�mero de datos
  vector [N] y; // Valores de los datos distribucion 
}

// Los par�metros aceptados por el modelo. Nuestro modelo
// acepta dos par�metros 'mu' and 'sigma'.
parameters {
  real mu; // Esperanza
  real<lower = 0> sigma2; // Desviaci�n est�ndar
}

// El modelo que ser� estimado. Modelamos la salida
// 'X' como una distribucipon normal con media 'mu'
// y desviaci�n est�ndar 'sigma'.
model {
  y ~ normal(mu,sqrt(sigma2));
  mu ~ normal(0, 1);
  sigma2 ~ normal(0, 1);
}


