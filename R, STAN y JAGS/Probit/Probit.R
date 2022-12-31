# Load necessary packages
library(rstan)  # Stan interface for R
library(Zelig)  # Data

# Load data set
data(turnout)

df = list(N = nrow(turnout), vote = turnout$vote, 
          educate = turnout$educate, income = turnout$income,
          age = turnout$age, agesq = turnout$age^2,
          nu = 4)     

# adjust nu as desired degree of freedom
fit1 = stan(model_code = "Probit.stan", data = df, chains = 4, iter = 1000)
print(fit1)