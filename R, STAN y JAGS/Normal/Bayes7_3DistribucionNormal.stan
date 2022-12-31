//
// Este programa de Stan define un modelo simple, con un 
// vector de valores 'X[N]' modelado como la distribución normal
// con esperanza 'mu' y desviación estándar 'sigma'.

// Los datos de entrada son un vector 'X[N]' de longitud 'N'.

data {
  int<lower = 0> N; // Número de datos
  vector [N] y; // Valores de los datos distribucion 
}

// Los parámetros aceptados por el modelo. Nuestro modelo
// acepta dos parámetros 'mu' and 'sigma'.
parameters {
  real mu; // Esperanza
  real<lower = 0> sigma2; // Desviación estándar
}

// El modelo que será estimado. Modelamos la salida
// 'X' como una distribucipon normal con media 'mu'
// y desviación estándar 'sigma'.
model {
  y ~ normal(mu,sqrt(sigma2));
  mu ~ normal(0, 1);
  sigma2 ~ normal(0, 1);
}


