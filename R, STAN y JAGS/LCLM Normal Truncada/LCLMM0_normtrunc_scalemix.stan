/*
Bayesian 
Truncated Normal Distribution
*/

functions{
// NO funciona
real normal_lub_rng(real mu, real sigma, real lb, real ub) {
  real p_lb = normal_cdf(lb, mu, sigma);
  real p_ub = normal_cdf(ub, mu, sigma);
  real u = uniform_rng(p_lb, p_ub);
  real y = mu + sigma * Phi(u);
  return y;
}
// SI funciona bien
real normal_lb_ub_rng(real mu, real sigma, real lb, real ub) {
    real p1 = normal_cdf(lb, mu, sigma);  // cdf with lower bound
    real p2 = normal_cdf(ub, mu, sigma);  // cdf with upper bound
    real u = uniform_rng(p1, p2);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf 
}
}


data {
int<lower=0> n; // number of subjects
real outcome[n]; // observed outcome
real cota; // cota for truncated distribution, I[cota,infinity]
}

parameters {
real mu; 
real<lower=0> sigma; 
real<lower=0> U[n];
}

transformed parameters{
real outa[n];
real outb[n];
for(i in 1:n){
      outa[i] = max([mu-sigma*sqrt(U[i]),cota]); 
      outb[i] = mu+sigma*sqrt(U[i]); 
} 
}

model{
// PRIOR
mu ~ normal(0,1000); 
sigma ~ gamma(0.01,0.01); 
// LIKELIHOOD
for(i in 1:n){
      U[i] ~ gamma(1.5,0.5);
      outcome[i] ~ uniform(outa[i],outb[i]); } 
}

// PREDICTIONS
generated quantities{
real pred[3];
real<lower=0> v;
real preda; 
real predb; 
pred[1]= normal_lub_rng( mu , sigma, cota, 1000) ; // No funciona
pred[2]= normal_lb_ub_rng( mu , sigma, cota, 1000) ; // SI funciona bien 
v = gamma_rng(1.5,0.5); 
preda = max([mu-sigma*sqrt(v),cota]); 
predb = mu+sigma*sqrt(v); 
pred[3]= uniform_rng(preda,predb); // SI funciona bien 
}

