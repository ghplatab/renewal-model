// =============================================================================
//  Ecuacion de Renovacion Bayesiana - COVID-19 Bogota
//  Modelo: NB-2 con AR(1) sobre log(Rt)
//  Exogeno: Exponencial(lambda_mu)
//
//  Proposito: Fase 2 con IS China-Gamma. Comparar la familia del componente
//  exogeno (Gamma vs Exponencial) SOBRE EL PROCESO DEL MODELO DEFINITIVO,
//  cambiando un solo componente.
//
//  Construido a partir de renewal_bogota_nb2_ar1_gamma.stan. UNICAS diferencias:
//    - data: se reemplazan alpha_mu, beta_mu por lambda_mu
//    - model: mu_exo[t] ~ exponential(lambda_mu) en lugar de gamma(alpha_mu, beta_mu)
//  El exogeno Exponencial es identico al de renewal_bogota_nb2_ar2_exponencial.stan
//  (lambda_mu = 0.120 = 1/8.33). Todo lo demas es identico al AR(1) + Gamma.
// =============================================================================

data {
  int<lower=1> T;
  int<lower=1> n_exogeno;
  int<lower=1> max_si;
  array[T] int<lower=0> y;
  vector[max_si] g;
  real<lower=0> lambda_mu;   // rate Exponencial = 1/media = 0.120
}

parameters {
  // --- AR(1) sobre log(Rt) ---
  real mu_rt;
  real<lower=0, upper=1> rho1;
  real<lower=0> sigma_epsilon;
  vector[T] epsilon_raw;

  // --- Exogeno: Exponencial ---
  vector<lower=0>[n_exogeno] mu_exo;

  // --- Sobredispersion NB-2 ---
  real<lower=0> phi;
}

transformed parameters {
  vector[T] epsilon;
  vector<lower=0>[T] Rt;
  vector<lower=0>[T] mu;
  vector<lower=0>[T] f;

  // --- AR(1) estacionario (parametrizacion no centrada) ---
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
  // --- Priors AR(1) (identicos al AR(1) + Gamma) ---
  mu_rt         ~ normal(0.0, 0.3);
  rho1          ~ normal(0.7, 0.15);
  sigma_epsilon ~ normal(0.0, 0.2) T[0,];
  epsilon_raw   ~ std_normal();

  // --- Prior exogeno: Exponencial (unico cambio) ---
  for (t in 1:n_exogeno)
    mu_exo[t] ~ exponential(lambda_mu);

  // --- Sobredispersion NB-2 ---
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
