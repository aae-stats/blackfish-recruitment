data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> nage;
  int<lower=0> noage;
  int<lower=1> ncohort;
  int<lower=1> nsystem;
  int y[N];
  vector[N] length_mm;
  matrix[N, K] X;
  vector[nage] age;
  int<lower=0,upper=N> age_idx[nage];          // obs with age known
  int<lower=0,upper=N> age_idy[noage];         // obs with age unknown
  int<lower=1,upper=nsystem> system[N];
  int<lower=1,upper=ncohort> cohort[N];
  vector[N] log_effort;
}

transformed data {
  matrix[N, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  vector[N] log_length_mm;
  
  // thin and scale the QR decomposition of covariates
  Q_ast = qr_Q(X)[, 1:K] * sqrt(N - 1);
  R_ast = qr_R(X)[1:K, ] / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
  
  // log transform length data
  log_length_mm = log(length_mm);
  
}

parameters {
  
  // std normal intercept and coef effects
  real zalpha;
  vector[K] ztheta;
  real zdelta;
  real zepsilon;
  
  // and catch curve slope (mortality term)
  vector<upper=0>[nsystem] zeta;

  // std normal random effects
  vector[nsystem] zgamma_system;
  vector[ncohort] zgamma_cohort;
  
  // and their SDs
  real<lower=0> sigma_system;
  real<lower=0> sigma_cohort;
  
  // and SD of age observations
  real<lower=0> sigma_age;
  
}

transformed parameters {
  vector[N] age_est;
  real alpha;
  vector[K] theta;
  vector[N] theta_term;
  real delta;
  real epsilon;
  vector[nsystem] gamma_system;
  vector[ncohort] gamma_cohort;
  vector[N] mu;
  vector[nage] mu_age;

  // scale main effects
  alpha = 5. * zalpha;
  theta = 5. * ztheta;
  delta = 2. * zdelta;
  epsilon = 2. * zepsilon;
  
  // scale random effects
  gamma_system = sigma_system * zgamma_system;
  gamma_cohort = sigma_cohort * zgamma_cohort;
  
  // precalculate coef effects for all obs
  theta_term =  Q_ast * theta;
  
  // estimate age from lengths based on regression of known-age 
  //   individuals against their length
  mu_age = delta + epsilon * log_length_mm[age_idx];
  age_est[age_idx] = age;
  age_est[age_idy] = delta + epsilon * log_length_mm[age_idy];
  
  // define linear predictor
  for (i in 1:N) {
    mu[i] = alpha + 
    zeta[system[i]] * age_est[i] +
    gamma_system[system[i]] +
    gamma_cohort[cohort[i]] +
//    gamma_site[site[i]] +   // need this to get correct slope estimate? Should be based on counts of binned data...
//    gamma_survey_year[survey_year[i]] + // and this?
    theta_term[i] + 
    log_effort[i];
  }

}

model {

  // priors for main params
  zalpha ~ std_normal();
  ztheta ~ std_normal();
  zdelta ~ std_normal();
  zepsilon ~ std_normal();

  // and for mortality parameter
  zeta ~ std_normal();
  
  // priors for random effects
  zgamma_system ~ std_normal();
  zgamma_cohort ~ std_normal();
  
  // and their SDs
  sigma_system ~ std_normal();
  sigma_cohort ~ std_normal();
  
  // and SD of age observations
  sigma_age ~ std_normal();

  // increment loglik for all count obs
  y ~ poisson_log(mu);
  
  // increment loglik for age-length obs
  age ~ normal(mu_age, sigma_age);
  
}

generated quantities {
  vector[K] beta;
  
  // coefficients rotated back to X space
  beta = R_ast_inverse * theta;
  
}
