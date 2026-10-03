// =============================================================================
//  Ecuación de Renovación Bayesiana — COVID-19 Bogotá
//  Modelo: NB-2 con AR(1) sobre log(Rt)
//  Exógeno: Gamma(alpha_mu, beta_mu) calibrado por MoM
//
//  Propósito: comparar AR(1) vs AR(2) en análisis de sensibilidad Fase 2
//
//  Diferencia respecto a renewal_bogota_nb2_ar2_gamma.stan:
//    - Eliminado: rho2, epsilon[t-2]
//    - AR(1): epsilon[t] = rho1*epsilon[t-1] + sigma_epsilon*epsilon_raw[t]
//    - Varianza estacionaria: sigma_stat = sigma_epsilon / sqrt(1 - rho1^2)
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
  // --- AR(1) sobre log(Rt) ---
  real mu_rt;
  real<lower=0, upper=1> rho1;
  real<lower=0> sigma_epsilon;
  vector[T] epsilon_raw;

  // --- Exógeno: Gamma ---
  vector<lower=0>[n_exogeno] mu_exo;

  // --- Sobredispersión NB-2 ---
  real<lower=0> phi;
}

transformed parameters {
  vector[T] epsilon;
  vector<lower=0>[T] Rt;
  vector<lower=0>[T] mu;
  vector<lower=0>[T] f;

  // --- AR(1) estacionario (parametrización no centrada) ---
  {
    real sigma_stat = sigma_epsilon / sqrt(fmax(1 - rho1^2, 1e-6));
    epsilon[1] = sigma_stat * epsilon_raw[1];
    for (t in 2:T)
      epsilon[t] = rho1 * epsilon[t-1] + sigma_epsilon * epsilon_raw[t];
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
  // --- Priors AR(1) ---
  mu_rt         ~ normal(0.0, 0.3);
  rho1 ~ normal(0.7, 0.15);  // prior centrado en autocorrelación alta
  sigma_epsilon ~ normal(0.0, 0.2) T[0,];
  epsilon_raw   ~ std_normal();

  // --- Prior exógeno: Gamma(alpha_mu, beta_mu) ---
  for (t in 1:n_exogeno)
    mu_exo[t] ~ gamma(alpha_mu, beta_mu);

  // --- Sobredispersión NB-2 ---
  phi ~ normal(0.0, 5.0) T[0,];

  // --- Verosimilitud ---
  for (t in 1:T)
    y[t] ~ neg_binomial_2(f[t], phi);
}

generated quantities {
  vector[T] log_lik;
  vector[T] y_pred;
  array[T] int y_rep;
  vector[T] log_Rt;

  for (t in 1:T) {
    real f_safe = fmin(f[t], 1e7);
    log_lik[t] = neg_binomial_2_lpmf(y[t] | f[t], phi);
    y_pred[t]  = f[t];
    y_rep[t]   = neg_binomial_2_rng(f_safe, phi);
    log_Rt[t]  = mu_rt + epsilon[t];
  }
}
