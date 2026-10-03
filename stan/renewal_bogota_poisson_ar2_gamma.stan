// =============================================================================
//  Ecuación de Renovación Bayesiana — COVID-19 Bogotá
//  Modelo: Poisson con AR(2) — baseline para LOO-CV vs NB-2
//  Exógeno: Gamma(alpha_mu, beta_mu)
// =============================================================================

data {
  int<lower=1> T;
  int<lower=1> n_exogeno;
  int<lower=1> max_si;
  array[T] int<lower=0> y;
  vector[max_si] g;
  real<lower=0> alpha_mu;
  real<lower=0> beta_mu;
}

parameters {
  real mu_rt;
  real<lower=0, upper=1> rho1;
  real<lower=0, upper=1> rho2;
  real<lower=0> sigma_epsilon;
  vector[T] epsilon_raw;
  vector<lower=0>[n_exogeno] mu_exo;
}

transformed parameters {
  vector[T] epsilon;
  vector<lower=0>[T] Rt;
  vector<lower=0>[T] mu;
  vector<lower=0>[T] f;

  {
    real denom = 1 - rho1^2 - rho2^2 -
                 2 * rho1^2 * rho2 / (1 - rho2);
    real sigma_stat = sigma_epsilon / sqrt(fmax(denom, 1e-6));

    epsilon[1] = sigma_stat * epsilon_raw[1];
    epsilon[2] = rho1 * epsilon[1] + sigma_epsilon * epsilon_raw[2];
    for (t in 3:T)
      epsilon[t] = rho1*epsilon[t-1] + rho2*epsilon[t-2] +
                   sigma_epsilon * epsilon_raw[t];
  }

  for (t in 1:T)
    Rt[t] = exp(mu_rt + epsilon[t]);

  for (t in 1:T)
    mu[t] = (t <= n_exogeno) ? mu_exo[t] : 0.0;

  for (t in 1:T) {
    real conv = 0.0;
    int lag_max = min(t - 1, max_si);
    for (tau in 1:lag_max)
      conv += f[t - tau] * g[tau];
    f[t] = fmax(mu[t] + Rt[t] * conv, 1e-6);
  }
}

model {
  mu_rt         ~ normal(0.0, 0.3);
  rho1          ~ normal(0.7, 0.15);
  rho2          ~ normal(0.1, 0.10);
  sigma_epsilon ~ normal(0.0, 0.2) T[0,];
  epsilon_raw   ~ std_normal();

  for (t in 1:n_exogeno)
    mu_exo[t] ~ gamma(alpha_mu, beta_mu);

  for (t in 1:T)
    y[t] ~ poisson(f[t]);
}

generated quantities {
  vector[T] log_lik;
  vector[T] y_pred;
  array[T] int y_rep;
  vector[T] log_Rt;

  for (t in 1:T) {
    real f_safe = fmin(f[t], 1e7);
    log_lik[t] = poisson_lpmf(y[t] | f[t]);
    y_pred[t]  = f[t];
    y_rep[t]   = poisson_rng(f_safe);
    log_Rt[t]  = mu_rt + epsilon[t];
  }
}
