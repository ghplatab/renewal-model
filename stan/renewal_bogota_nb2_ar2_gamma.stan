// =============================================================================
//  Ecuación de Renovación Bayesiana — COVID-19 Bogotá
//  Modelo: Negative Binomial tipo 2 (NB-2) con AR(2) sobre log(Rt)
//  Exógeno: Gamma(alpha_mu, beta_mu) — calibrado por MoM sobre datos reales
//
//  Basado en: Mishra et al. (2020), arXiv:2006.16487
//  Período: 14-Mar-2020 → 09-Sep-2020 (T = 180 días)
//
//  Cambios respecto a renewal_bogota_nb2_ar2.stan (Rayleigh):
//    - Exógeno: mu_exo[t] ~ gamma(alpha_mu, beta_mu) en lugar de rayleigh(sigma_mu)
//    - alpha_mu = 1.401, beta_mu = 0.168 (MoM: media=8.33, var=49.52)
//    - Eliminado: sigma_mu (parámetro Rayleigh)
//    - Moda Gamma = (alpha-1)/beta = 0.401/0.168 = 2.39 > 0 (flujo continuo)
//
//  Priors AR(2):
//    rho1 ~ N(0.7, 0.15) T[0,1]
//    rho2 ~ N(0.1, 0.10) T[0,1]
// =============================================================================

data {
  int<lower=1> T;
  int<lower=1> n_exogeno;
  int<lower=1> max_si;

  array[T] int<lower=0> y;
  vector[max_si] g;

  // Hiperparámetros Gamma exógeno (pasados como datos para flexibilidad)
  real<lower=0> alpha_mu;   // shape = 1.401
  real<lower=0> beta_mu;    // rate  = 0.168
}

parameters {
  // --- AR(2) sobre log(Rt) ---
  real mu_rt;
  real<lower=0, upper=1> rho1;
  real<lower=0, upper=1> rho2;
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

  // --- AR(2) estacionario (parametrización no centrada) ---
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
  // --- Priors AR(2) ---
  mu_rt         ~ normal(0.0, 0.3);
  rho1          ~ normal(0.7, 0.15);
  rho2          ~ normal(0.1, 0.10);
  sigma_epsilon ~ normal(0.0, 0.2) T[0,];
  epsilon_raw   ~ std_normal();

  // --- Prior exógeno: Gamma(alpha_mu, beta_mu) ---
  // Hiperparámetros calibrados por MoM sobre datos reales (14-25 mar 2020):
  //   media_obs = 8.33, var_obs = 49.52
  //   alpha = media^2/var = 1.401, beta = media/var = 0.168
  // Moda = (alpha-1)/beta = 2.39 > 0 — coherente con hub aeroportuario
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
