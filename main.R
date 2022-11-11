# analysis of blackfish age-length data with generalised catch curve model

# load packages
library(dplyr)
library(rstan)

# load helpers


# load data
## NOTES: need to define covars based on year of recruitment
## WIll need to define an assumed cohort for each individual, then
##   use this to define cohort and X
## NEED TO THINK ABOUT model structure: catch curve has counts of each age class? Can also use PP model to set it up on continuous line?


N <- 200
K <- 4
nage <- 60
no_age <- N - nage
alpha <- 5
beta <- rnorm(K, sd = 0.1)
zeta <- -2
X <- matrix(rnorm(N * K), nrow = N)
age_idx <- sample(seq_len(N), size = nage, replace = FALSE)
age_idy <- seq_len(N)[-age_idx]
length_mm <- sample(seq_len(250)[51:250], size = N, replace = TRUE)
nsystem <- 3
ncohort <- 10
system <- sample(seq_len(nsystem), size = N, replace = TRUE)
cohort <- sample(seq_len(ncohort), size = N, replace = TRUE)
gamma_sys <- rnorm(nsystem, sd = 0.1)
gamma_coh <- rnorm(ncohort, sd = 0.1)
log_effort <- log(rnorm(N, mean = 900, sd = 100))
age <- round(10 * log(length_mm) - 38)

mu <- alpha + X %*% beta + zeta * age + log_effort #+ gamma_sys[system] + gamma_coh[cohort]
y <- rpois(N, lambda = exp(mu))

# format data
dat <- list(
  N = N,
  K = K,
  nage = nage,
  noage = no_age,
  y = y,
  X = X,
  age = age[age_idx],
  age_idx = age_idx,
  age_idy = age_idy,
  length_mm = length_mm,
  nsystem = nsystem,
  system = system,
  ncohort = ncohort,
  cohort = cohort,
  log_effort = log_effort
)

# prepare analysis inputs

# compile model
mod <- stan_model("src/catch-curve.stan")

# sample
draws <- sampling(
  mod,
  data = dat,
  chains = 4,
  iter = 1000,
  warmup = 500,
  control = list(max_treedepth = 20)
)


# summarise