/*
Bayesian 
Ordinal Regression
*/
 
data {
  int<lower=0> n; // number of subjects
  real X1[n]; // predictor variables
  real X2[n]; // predictor variables
  int<lower=0> J; // categories
  int outcome[n]; // observed outcome
}

parameters {
real theta0; 
real theta1; 
real theta2;  
ordered[J-1] kappa; // (internal) cut points, thresholds
}

transformed parameters {vector[J-1] kappa1; 
vector[n] mu; 
kappa1 = kappa - kappa[1];
for(i in 1:n){
    mu[i] = theta0 +theta1*X1[i] +theta2*X2[i] ; 
}
}


model{ 
kappa ~ normal(0,1000);
for (j in 1:J-1){
  kappa[j] ~ normal(j-1 + .05, 100);
}theta0 ~ normal(0,1000); 
theta1 ~ normal(0,1000); 
theta2 ~ normal(0,1000); 

for(i in 1:n){
  target += ordered_probit_lpmf(outcome[i]| mu[i] , kappa1); // 
} 
}

