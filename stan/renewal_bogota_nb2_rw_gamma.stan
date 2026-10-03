// =============================================================================
//  Ecuacion de Renovacion Bayesiana - COVID-19 Bogota
//  Modelo: NB-2 con RANDOM WALK sobre log(Rt)
//  Exogeno: Gamma(alpha_mu, beta_mu), identico al modelo AR(1)
//
//  Proposito: responder la observacion 4 del revisor. Comparar el AR(1)
//  estacionario contra un random walk (rho1 = 1), por capacidad predictiva.
//
//  Diferencias respecto a renewal_bogota_nb2_ar1_gamma.stan:
//
//   1. NO hay mu_rt. En un random walk el nivel global no es identificable
//      por separado del estado inicial: log Rt = log Rt_1 + suma de
//      innovaciones. No existe media de largo plazo a la que revertir.
//      Esa ausencia ES el punto del modelo, no una omision.
//
//   2. NO hay inicializacion estacionaria (sigma_stat). Un random walk no
//      tiene distribucion estacionaria, asi que el estado inicial necesita
//      su propia previa.
//
//   3. NO hay fmax sobre el radicando, porque no hay radicando: desaparece
//      la fuente de gradiente discontinuo del AR(1).
//
//  Lo que se mantiene IGUAL para que la comparacion sea limpia:
//   - parametrizacion no centrada (epsilon_raw ~ std_normal)
//   - previa sigma_epsilon ~ N+(0, 0.2)
//   - previa phi ~ N+(0, 5)
//   - componente exogeno Gamma con los mismos hiperparametros
//   - verosimilitud NB-2
//   - bloque generated quantities identico
//
//  Sobre la previa del estado inicial: en el AR(1),
//    log Rt[1] = mu_rt + sigma_stat * z,  con mu_rt ~ N(0, 0.3)
//    y sigma_stat = 0.099 / sqrt(1 - 0.962^2) ~ 0.36
//  lo que implica log Rt[1] ~ N(0, ~0.47). Se usa N(0, 0.5) para que el
//  estado inicial tenga una previa comparable y la comparacion sea justa.
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
  // --- Random walk sobre log(Rt) ---
  real log_Rt1;                    // nivel inicial (reemplaza a mu_rt)
  real<lower=0> sigma_epsilon;     // escala de las innovaciones
  vector[T - 1] epsilon_raw;       // innovaciones estandarizadas

  // --- Exogeno: Gamma ---
  vector<lower=0>[n_exogeno] mu_exo;

  // --- Sobredispersion NB-2 ---
  real<lower=0> phi;
}

transformed parameters {
  vector[T] log_Rt;
  vector<lower=0>[T] Rt;
  vector<lower=0>[T] mu;
  vector<lower=0>[T] f;

  // --- Random walk (parametrizacion no centrada, sin clamps) ---
  log_Rt[1] = log_Rt1;
  for (t in 2:T)
    log_Rt[t] = log_Rt[t-1] + sigma_epsilon * epsilon_raw[t-1];

  for (t in 1:T)
    Rt[t] = exp(log_Rt[t]);

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
  // --- Previa del estado inicial ---
  log_Rt1       ~ normal(0.0, 0.5);

  // --- Innovaciones (identicas al AR(1)) ---
  sigma_epsilon ~ normal(0.0, 0.2) T[0,];
  epsilon_raw   ~ std_normal();

  // --- Exogeno Gamma (identico al AR(1)) ---
  for (t in 1:n_exogeno)
    mu_exo[t] ~ gamma(alpha_mu, beta_mu);

  // --- Sobredispersion NB-2 (identica al AR(1)) ---
  phi ~ normal(0.0, 5.0) T[0,];

  // --- Verosimilitud ---
  for (t in 1:T)
    y[t] ~ neg_binomial_2(f[t], phi);
}

generated quantities {
  vector[T] log_lik;
  vector[T] y_pred;
  array[T] int y_rep;

  for (t in 1:T) {
    real f_safe = fmin(f[t], 1e7);
    log_lik[t] = neg_binomial_2_lpmf(y[t] | f[t], phi);
    y_pred[t]  = f[t];
    y_rep[t]   = neg_binomial_2_rng(f_safe, phi);
  }
}
