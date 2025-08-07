data {
  int<lower=1> N;      // días observados
  int<lower=1> N2;     // días para exógeno
  int cases[N];        // casos diarios observados
  real SI[N];          // intervalos seriales fijos (arreglo de longitud N)
}
transformed data {
  vector[N] SI_rev;
  for (i in 1:N) {
    SI_rev[i] = SI[N - i + 1];
  }
}
parameters {
  real<lower=0> phi;
  vector[N+1] weekly_effect;
  real<lower=0, upper=1> weekly_rho;
  real<lower=0, upper=1> weekly_rho1;
  real<lower=0> weekly_sd;

  // --- Nueva sección para mu usando Rayleigh ---
  real<lower=0> sigma_mu;      // escala de la Rayleigh
  vector<lower=0>[N2] mu;      // importaciones diarias
}
transformed parameters {
  vector[N] prediction = rep_vector(1e-5, N);
  vector<lower=0>[N] Rt;

  // Construcción de Rt y la predicción
  Rt[1:N] = exp(weekly_effect[1:N]);

  // Inicializa predicción con parte exógena en los primeros N2 días
  prediction[1:N2] = prediction[1:N2] + mu[1:N2];

  for (i in 2:N) {
    real convolution = dot_product(prediction[1:(i-1)], tail(SI_rev, i-1));
    prediction[i] = prediction[i] + Rt[i] * convolution;
  }
}
model {
  // --- Priors para el componente AR(2) semanales (igual que antes) ---
  weekly_sd    ~ normal(0, 1);
  weekly_rho   ~ normal(0.8, 0.05);
  weekly_rho1  ~ normal(0.1, 0.05);
  phi          ~ normal(0, 5);

  // Priors AR(2) sobre weekly_effect
  weekly_effect[3:(N+1)] ~ normal(
    weekly_effect[2:N] * weekly_rho + weekly_effect[1:(N-1)] * weekly_rho1,
    weekly_sd * sqrt(
      1 - pow(weekly_rho, 2) - pow(weekly_rho1, 2)
      - 2 * pow(weekly_rho, 2) * weekly_rho1 / (1 - weekly_rho1)
    )
  );
  weekly_effect[2] ~ normal(-1, weekly_sd * sqrt(
    1 - pow(weekly_rho, 2) - pow(weekly_rho1, 2)
    - 2 * pow(weekly_rho, 2) * weekly_rho1 / (1 - weekly_rho1)
  ));
  weekly_effect[1] ~ normal(-1, 0.1);

  // --- NUEVO: Priors para el componente exógeno mu[t] usando Rayleigh ---
  sigma_mu ~ normal(0, 2);       // half-normal para el parámetro de escala
  for (t in 1:N2) {
    mu[t] ~ rayleigh(sigma_mu);
  }

  // --- Likelihood de los casos diarios con Neg-Binomial ---
  cases ~ neg_binomial_2(prediction, phi);
}

generated quantities {
  vector[N] y_rep;
  for (t in 1:N) {
    y_rep[t] = neg_binomial_2_rng(prediction[t], phi);
  }
}
