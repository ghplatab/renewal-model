# =============================================================================
#  FASE 2 - AJUSTE CON EL INTERVALO SERIAL CHINA-GAMMA (Bi et al., 2020)
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  El modelo definitivo adopta el intervalo serial China-Gamma,
#  Gamma(2.29, 0.36), de media 6.36 dias, tomado directamente de la
#  distribucion reportada por Bi et al. (2020).
#
#  Diseno: el modelo definitivo es la base y cada alternativa modifica un
#  unico componente, de modo que las diferencias observadas sean atribuibles
#  a ese componente.
#
#      M0  AR(1)         + Gamma        modelo definitivo
#      M1  AR(2)         + Gamma        orden del proceso latente
#      M2  Paseo aleat.  + Gamma        reversion a la media
#      M3  AR(1)         + Exponencial  familia del componente exogeno
#
#  Configuracion comun a los cuatro: adapt_delta = 0.99, max_treedepth = 12,
#  4 cadenas x (1500 calentamiento + 1500 muestreo), seed = 42.
#
#  Ejecucion:
#    - La configuracion de rutas va siempre primero.
#    - La verificacion predictiva a priori no usa Stan (un par de minutos).
#    - Los cuatro ajustes tardan unos 35-50 minutos cada uno y son
#      independientes: si un ajuste ya esta guardado se carga en vez de
#      reajustarse, de modo que la sesion puede interrumpirse y retomarse.
#    - Los dos ultimos bloques comparan los cuatro modelos y no usan Stan.
#
#  Los CSV crudos de CmdStan se escriben en C:/stan_out: las rutas con tildes
#  corrompen su salida en Windows.
# =============================================================================


# =============================================================================
#  PASO 0 - CONFIGURACION
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(cmdstanr)
  library(posterior); library(loo); library(bayesplot); library(ggplot2)
})

CMDSTAN_DIR <- Sys.getenv("CMDSTAN", unset = "C:/cmdstan/cmdstan-2.39.0")
if (!dir.exists(CMDSTAN_DIR)) stop("No encuentro CmdStan en: ", CMDSTAN_DIR)
set_cmdstan_path(CMDSTAN_DIR)

BASE      <- file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Escritorio",
                       "Tesis Maestría", "AJUSTES_2026", "AJUSTES_PROFE_JCS")
STAN_ROOT <- file.path(BASE, "Documentos George", "MODELO_STAN")
SEP       <- file.path(BASE, "AJUSTES_SEP_2026")
DEF_ORIG  <- file.path(STAN_ROOT, "MODELO_REPO_GITHUB", "MODELO_STAN_DEFINITIVO")

OUT_OLD <- file.path(SEP, "resultados")                    # corridas viejas (10.5 d)
OUT     <- file.path(SEP, "resultados", "fase2_chinagamma")
FIG     <- file.path(SEP, "figuras",    "fase2_chinagamma")
STAN    <- file.path(SEP, "stan")
for (p in c(OUT, FIG, STAN)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

BUILD <- "C:/stan_build"; CSVOUT <- "C:/stan_out"
dir.create(BUILD,  showWarnings = FALSE, recursive = TRUE)
dir.create(CSVOUT, showWarnings = FALSE, recursive = TRUE)

# --- Modelos Stan -------------------------------------------------------------
# AR(1)+Gamma y RW ya estan en AJUSTES_SEP_2026/stan (copiados en el script 03).
# AR(2)+Gamma se copia sin modificar desde tu carpeta del modelo definitivo.
# AR(1)+Exponencial es nuevo: el AR(1)+Gamma con solo la linea del exogeno cambiada.
STAN_M0 <- file.path(STAN, "renewal_bogota_nb2_ar1_gamma.stan")
STAN_M1 <- file.path(STAN, "renewal_bogota_nb2_ar2_gamma.stan")
STAN_M2 <- file.path(STAN, "renewal_bogota_nb2_rw_gamma.stan")
STAN_M3 <- file.path(STAN, "renewal_bogota_nb2_ar1_exponencial.stan")
if (!file.exists(STAN_M0))
  file.copy(file.path(DEF_ORIG, "renewal_bogota_nb2_ar1_gamma.stan"), STAN_M0)
if (!file.exists(STAN_M1))
  file.copy(file.path(DEF_ORIG, "renewal_bogota_nb2_ar2_gamma.stan"), STAN_M1)

# --- Datos (identico a analisis_fase2_definitivo.R) --------------------------
DATOS  <- file.path(STAN_ROOT, "MODELO_STAN_4", "datos_agregados.xlsx")
df_raw <- as.data.frame(read_excel(DATOS, sheet = "NO_IMPORTADOS"))
names(df_raw) <- c("fecha", "casos")
df_raw <- df_raw[!grepl("nan", as.character(df_raw$fecha), ignore.case = TRUE), ]
df_raw <- df_raw[!is.na(df_raw$fecha), ]
df_raw$fecha <- as.Date(as.numeric(df_raw$fecha), origin = "1899-12-30")

FECHA_INICIO <- as.Date("2020-03-14"); FECHA_FIN <- as.Date("2020-09-09")
fechas <- seq(FECHA_INICIO, FECHA_FIN, by = "day")
dfp    <- left_join(data.frame(fecha = fechas), df_raw, by = "fecha")
dfp$casos[is.na(dfp$casos)] <- 0L
y_obs  <- as.integer(dfp$casos)

T_PERIODO <- 180L; N_EXOGENO <- 12L; MAX_SI <- 30L
ALPHA_MU  <- 1.401; BETA_MU <- 0.168     # exogeno Gamma (MoM)
LAMBDA_MU <- 0.120                       # exogeno Exponencial (1/8.33)

# --- Intervalo serial: CHINA-GAMMA (el unico cambio respecto a la version vieja)
discret_si <- function(pfun, ...) {
  g <- numeric(MAX_SI); g[1] <- pfun(1.5, ...) - pfun(0, ...)
  for (s in 2:MAX_SI) g[s] <- pfun(s + 0.5, ...) - pfun(s - 0.5, ...)
  g[g < 0] <- 0; g / sum(g)
}
g_cg  <- discret_si(pgamma, shape = 2.29, rate = 0.36)
g_old <- discret_si(pgamma, shape = 6.5,  rate = 0.62)   # solo para verificar

media_g <- function(g) sum(seq_along(g) * g)
de_g    <- function(g) sqrt(sum((seq_along(g) - media_g(g))^2 * g))
if (abs(media_g(g_cg) - 6.36) > 0.02)
  stop("El IS no tiene la media esperada (6.36). Revisa g_cg antes de seguir.")

datos_gamma <- list(T = T_PERIODO, n_exogeno = N_EXOGENO, max_si = MAX_SI,
                    y = y_obs, g = g_cg, alpha_mu = ALPHA_MU, beta_mu = BETA_MU)
datos_exp   <- list(T = T_PERIODO, n_exogeno = N_EXOGENO, max_si = MAX_SI,
                    y = y_obs, g = g_cg, lambda_mu = LAMBDA_MU)

# --- Configuracion comun del muestreador --------------------------------------
N_CHAINS <- 4L; N_WARMUP <- 1500L; N_SAMPLING <- 1500L
ADAPT_DELTA <- 0.99; MAX_TREEDEPTH <- 12L; SEED <- 42

# --- Valores iniciales (los mismos de tus scripts) ----------------------------
init_ar1 <- function() list(mu_rt = 0.0, rho1 = 0.8, sigma_epsilon = 0.1,
                            epsilon_raw = rep(0.0, T_PERIODO),
                            mu_exo = rep(ALPHA_MU / BETA_MU, N_EXOGENO), phi = 4.0)
init_ar2 <- function() list(mu_rt = 0.0, rho1 = 0.7, rho2 = 0.1, sigma_epsilon = 0.1,
                            epsilon_raw = rep(0.0, T_PERIODO),
                            mu_exo = rep(ALPHA_MU / BETA_MU, N_EXOGENO), phi = 4.0)
init_rw  <- function() list(log_Rt1 = 0.0, sigma_epsilon = 0.1,
                            epsilon_raw = rep(0.0, T_PERIODO - 1),
                            mu_exo = rep(ALPHA_MU / BETA_MU, N_EXOGENO), phi = 4.0)
init_exp <- function() list(mu_rt = 0.0, rho1 = 0.8, sigma_epsilon = 0.1,
                            epsilon_raw = rep(0.0, T_PERIODO),
                            mu_exo = rep(1 / LAMBDA_MU, N_EXOGENO), phi = 4.0)

# --- Funcion comun: ajusta (o carga si ya existe) y reporta diagnosticos -------
correr <- function(nombre, stan_file, datos, init_fn, pars) {
  rds <- file.path(OUT, paste0(nombre, ".rds"))
  if (file.exists(rds)) {
    cat(sprintf("  [%s] ya existe: se carga, no se reajusta\n", nombre))
    fit <- readRDS(rds)
  } else {
    cat(sprintf("  [%s] ajustando...\n", nombre))
    mod <- cmdstan_model(stan_file, dir = BUILD)
    t0  <- proc.time()
    fit <- mod$sample(
      data = datos, iter_warmup = N_WARMUP, iter_sampling = N_SAMPLING,
      chains = N_CHAINS, parallel_chains = N_CHAINS,
      adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH,
      init = init_fn, seed = SEED, refresh = 250,
      show_messages = FALSE, output_dir = CSVOUT)
    cat(sprintf("  Tiempo: %.1f min\n", (proc.time() - t0)[["elapsed"]] / 60))
    fit$save_object(rds)
    cat(sprintf("  Guardado: %s\n", basename(rds)))
  }
  s  <- fit$summary(pars)
  dg <- fit$diagnostic_summary(quiet = TRUE)
  print(s[, c("variable", "mean", "median", "q5", "q95", "rhat", "ess_bulk", "ess_tail")])
  cat(sprintf("\n  Divergencias : %d de %d (%.2f%%)\n", sum(dg$num_divergent),
              N_CHAINS * N_SAMPLING, 100 * sum(dg$num_divergent) / (N_CHAINS * N_SAMPLING)))
  cat(sprintf("  Treedepth max: %d\n", sum(dg$num_max_treedepth)))
  cat(sprintf("  Rhat max     : %.4f   (criterio <= 1.01)\n", max(s$rhat)))
  cat(sprintf("  ESS bulk min : %.0f     (criterio > 400)\n", min(s$ess_bulk)))
  cat(sprintf("  E-BFMI       : %s   (< 0.3 es problematico)\n",
              paste(sprintf("%.2f", dg$ebfmi), collapse = "  ")))
  invisible(fit)
}

PARS_AR1 <- c("mu_rt", "rho1", "sigma_epsilon", "phi")
PARS_AR2 <- c("mu_rt", "rho1", "rho2", "sigma_epsilon", "phi")
PARS_RW  <- c("log_Rt1", "sigma_epsilon", "phi")

cat("=== PASO 0 OK ===\n")
cat(sprintf("  R %s | CmdStan %s\n", getRversion(), cmdstan_version()))
cat(sprintf("  Datos: %d dias, %s casos\n", length(y_obs), format(sum(y_obs), big.mark = ",")))
cat(sprintf("  IS China-Gamma : media %.2f  DE %.2f   (debe ser 6.36 / 4.17)\n",
            media_g(g_cg), de_g(g_cg)))
cat(sprintf("  IS viejo       : media %.2f  DE %.2f   (solo referencia, NO se usa)\n",
            media_g(g_old), de_g(g_old)))
cat(sprintf("  Stan M0 AR1+Gamma: %s | M1 AR2+Gamma: %s | M2 RW: %s | M3 AR1+Exp: %s\n",
            file.exists(STAN_M0), file.exists(STAN_M1), file.exists(STAN_M2), file.exists(STAN_M3)))


# =============================================================================
#  PASO 1 - PRIOR PREDICTIVE CHECK DEL DEFINITIVO   (obs. 3a)     ~2 min, sin Stan
#
#  Igual al PASO 8 del script 02, pero con el IS China-Gamma. Simula de las
#  previas del AR(1)+Gamma y mira que epidemias implican, antes de ver datos.
# =============================================================================

cat("\n=== PASO 1: PRIOR PREDICTIVE CHECK (IS China-Gamma) ===\n")

rnorm_trunc <- function(mean, sd, lo, hi) {
  repeat { v <- rnorm(1, mean, sd); if (v > lo && v < hi) return(v) }
}

set.seed(2026)
N_SIM <- 1000
filas <- vector("list", N_SIM); n_explotan <- 0

for (s in 1:N_SIM) {
  mu_rt  <- rnorm(1, 0, 0.3)
  rho1   <- rnorm_trunc(0.7, 0.15, 0, 1)
  sigma  <- abs(rnorm(1, 0, 0.2))
  phi    <- abs(rnorm(1, 0, 5))
  mu_exo <- rgamma(N_EXOGENO, shape = ALPHA_MU, rate = BETA_MU)

  eps <- numeric(T_PERIODO)
  eps[1] <- (sigma / sqrt(1 - rho1^2)) * rnorm(1)
  for (t in 2:T_PERIODO) eps[t] <- rho1 * eps[t - 1] + sigma * rnorm(1)
  Rt <- exp(mu_rt + eps)

  f <- numeric(T_PERIODO)
  for (t in 1:T_PERIODO) {
    mut  <- if (t <= N_EXOGENO) mu_exo[t] else 0
    lagm <- min(t - 1, MAX_SI)
    conv <- if (lagm >= 1) sum(f[t - (1:lagm)] * g_cg[1:lagm]) else 0
    f[t] <- max(mut + Rt[t] * conv, 1e-6)
  }
  if (!all(is.finite(f)) || max(f) > 1e9) { n_explotan <- n_explotan + 1; next }
  y <- rnbinom(T_PERIODO, mu = f, size = max(phi, 1e-3))
  if (!all(is.finite(y))) { n_explotan <- n_explotan + 1; next }

  filas[[s]] <- data.frame(Rt_ini = Rt[1], Rt_med = Rt[90], Rt_fin = Rt[180],
                           Rt_max = max(Rt), casos_tot = sum(y), casos_max = max(y))
}
sim <- do.call(rbind, filas)

cat(sprintf("  Simulaciones validas: %d de %d (%d explotaron)\n", nrow(sim), N_SIM, n_explotan))
q <- function(v) quantile(v, c(.025, .25, .5, .75, .975), names = FALSE)
cat("\n  Cuantiles 2.5 / 25 / 50 / 75 / 97.5 implicados por las previas:\n")
for (cl in names(sim))
  cat(sprintf("    %-10s %s\n", cl,
              paste(formatC(q(sim[[cl]]), format = "f", digits = 1, big.mark = ","),
                    collapse = "  ")))
cat(sprintf("\n  Observado: casos_tot = %s | casos_max = %s\n",
            format(sum(y_obs), big.mark = ","), format(max(y_obs), big.mark = ",")))
cat(sprintf("  Percentil previo del observado: casos_tot %.1f%% | casos_max %.1f%%\n",
            100 * mean(sim$casos_tot <= sum(y_obs)), 100 * mean(sim$casos_max <= max(y_obs))))

write.csv(sim, file.path(OUT, "f2cg_prior_predictive.csv"), row.names = FALSE)
cat("  Guardado: f2cg_prior_predictive.csv\n")


# =============================================================================
#  PASO 2 - M0: AR(1) + GAMMA   (MODELO DEFINITIVO)                 ~35-50 min
#
#  Para comparar: la corrida vieja (nucleo 10.5 d, adapt_delta 0.99) tuvo
#  58 divergencias (0.97%), rho1 = 0.962, sigma_epsilon = 0.099, phi = 3.26.
# =============================================================================

cat("\n=== PASO 2: M0 AR(1) + Gamma (definitivo) ===\n")
fit_m0 <- correr("f2cg_M0_ar1_gamma", STAN_M0, datos_gamma, init_ar1, PARS_AR1)

# Las figuras van dentro de try(): si una falla, el script sigue con los
# modelos siguientes (importante cuando se deja corriendo de noche).

# Pairs plot con divergencias marcadas (obs. 5 y hallazgo de la raiz unitaria)
try({
  p <- mcmc_pairs(fit_m0$draws(PARS_AR1), np = nuts_params(fit_m0),
                  off_diag_args = list(size = 0.6, alpha = 0.4))
  ggsave(file.path(FIG, "f2cg_M0_pairs.png"), p, width = 9, height = 9, dpi = 150)
  cat("  Guardado: f2cg_M0_pairs.png\n")
})

# Traceplots (Figura 3-9)
try({
  p <- mcmc_trace(fit_m0$draws(PARS_AR1), facet_args = list(ncol = 2))
  ggsave(file.path(FIG, "f2cg_M0_traceplots.png"), p, width = 10, height = 6, dpi = 150)
  ggsave(file.path(FIG, "f2cg_M0_traceplots.pdf"), p, width = 10, height = 6)
  cat("  Guardado: f2cg_M0_traceplots.png/.pdf\n")
})


# =============================================================================
#  PASO 3 - M1: AR(2) + GAMMA   (orden del proceso)                 ~40-60 min
#
#  Para comparar: el AR(2) China-Gamma exploratorio (adapt_delta 0.95,
#  1000 draws por cadena) dio rho1 = 0.766, rho2 = 0.151, sigma = 0.088,
#  phi = 3.245, elpd = -1232.5. Deberia salir parecido.
# =============================================================================

cat("\n=== PASO 3: M1 AR(2) + Gamma ===\n")
fit_m1 <- correr("f2cg_M1_ar2_gamma", STAN_M1, datos_gamma, init_ar2, PARS_AR2)

try({
  p <- mcmc_pairs(fit_m1$draws(PARS_AR2), np = nuts_params(fit_m1),
                  off_diag_args = list(size = 0.6, alpha = 0.4))
  ggsave(file.path(FIG, "f2cg_M1_pairs.png"), p, width = 10, height = 10, dpi = 150)
  cat("  Guardado: f2cg_M1_pairs.png\n")
})


# =============================================================================
#  PASO 4 - M2: RANDOM WALK + GAMMA   (reversion a la media, obs. 4)  ~30-45 min
# =============================================================================

cat("\n=== PASO 4: M2 Random walk + Gamma ===\n")
fit_m2 <- correr("f2cg_M2_rw_gamma", STAN_M2, datos_gamma, init_rw, PARS_RW)


# =============================================================================
#  PASO 5 - M3: AR(1) + EXPONENCIAL   (familia del exogeno)          ~35-50 min
# =============================================================================

cat("\n=== PASO 5: M3 AR(1) + Exponencial ===\n")
fit_m3 <- correr("f2cg_M3_ar1_exp", STAN_M3, datos_exp, init_exp, PARS_AR1)


# =============================================================================
#  PASO 6 - COMPARACION DE LOS CUATRO MODELOS                  ~5 min, sin Stan
#
#  Requiere que los PASOS 2-5 hayan terminado (carga los .rds guardados).
#    6a. Tabla de diagnosticos
#    6b. LOO-CV: loo_compare contra el definitivo, Pareto k
#    6c. Trayectorias de R_t: fechas clave, maximo, cruce bajo 1, diferencias
#        con el definitivo (sobre MEDIANAS, diferencia entre 2 modelos)
#    6d. Parametros de largo plazo del AR(1): media estacionaria y vida media
# =============================================================================

cat("\n=== PASO 6: COMPARACION M0-M3 ===\n")

modelos <- list(
  M0_AR1_Gamma = list(rds = "f2cg_M0_ar1_gamma", pars = PARS_AR1),
  M1_AR2_Gamma = list(rds = "f2cg_M1_ar2_gamma", pars = PARS_AR2),
  M2_RW_Gamma  = list(rds = "f2cg_M2_rw_gamma",  pars = PARS_RW),
  M3_AR1_Exp   = list(rds = "f2cg_M3_ar1_exp",   pars = PARS_AR1))
fits <- lapply(modelos, function(m) readRDS(file.path(OUT, paste0(m$rds, ".rds"))))

# --- 6a. Diagnosticos -----------------------------------------------------------
diag <- do.call(rbind, lapply(names(fits), function(nm) {
  fit <- fits[[nm]]; s <- fit$summary(modelos[[nm]]$pars)
  dg  <- fit$diagnostic_summary(quiet = TRUE)
  data.frame(modelo = nm, rhat_max = max(s$rhat), ess_bulk_min = min(s$ess_bulk),
             ess_tail_min = min(s$ess_tail), divergencias = sum(dg$num_divergent),
             pct_div = 100 * sum(dg$num_divergent) / (N_CHAINS * N_SAMPLING),
             treedepth_max = sum(dg$num_max_treedepth), ebfmi_min = min(dg$ebfmi))
}))
cat("\n--- 6a. Diagnosticos ---\n"); print(diag, digits = 4, row.names = FALSE)
write.csv(diag, file.path(OUT, "f2cg_diagnosticos.csv"), row.names = FALSE)

post <- do.call(rbind, lapply(names(fits), function(nm) {
  s <- fits[[nm]]$summary(modelos[[nm]]$pars, mean, median,
                          ~quantile(.x, probs = c(0.025, 0.975)))
  cbind(modelo = nm, as.data.frame(s))
}))
cat("\n--- Posteriores (media, mediana, IC95) ---\n"); print(post, digits = 3, row.names = FALSE)
write.csv(post, file.path(OUT, "f2cg_posteriores.csv"), row.names = FALSE)

# --- 6b. LOO-CV -----------------------------------------------------------------
cid  <- rep(1:N_CHAINS, each = N_SAMPLING)
loos <- lapply(fits, function(fit) {
  ll <- fit$draws("log_lik", format = "matrix")
  loo(ll, r_eff = relative_eff(exp(ll), chain_id = cid))
})
cmp <- loo_compare(loos)
cat("\n--- 6b. loo_compare (fila de arriba = mayor elpd) ---\n"); print(cmp)

# Diferencias SIEMPRE respecto al definitivo M0, no respecto al mejor
ll_m0 <- loos$M0_AR1_Gamma$pointwise[, "elpd_loo"]
tabla_loo <- do.call(rbind, lapply(names(loos), function(nm) {
  lo <- loos[[nm]]; d <- lo$pointwise[, "elpd_loo"] - ll_m0
  data.frame(modelo = nm,
             elpd = lo$estimates["elpd_loo", "Estimate"], se = lo$estimates["elpd_loo", "SE"],
             p_loo = lo$estimates["p_loo", "Estimate"],
             dif_vs_M0 = sum(d), se_dif_vs_M0 = sqrt(length(d)) * sd(d),
             k_max = max(lo$diagnostics$pareto_k), dias_k_07 = sum(lo$diagnostics$pareto_k > 0.7))
}))
cat("\n--- 6b. LOO respecto al definitivo M0 ---\n"); print(tabla_loo, digits = 4, row.names = FALSE)
write.csv(tabla_loo, file.path(OUT, "f2cg_loo.csv"), row.names = FALSE)
saveRDS(loos, file.path(OUT, "f2cg_loos.rds"))

# --- 6c. Trayectorias de R_t ------------------------------------------------------
Rt_med <- sapply(fits, function(fit) apply(fit$draws("Rt", format = "matrix"), 2, median))
FECHAS_CLAVE <- as.Date(c("2020-03-25", "2020-04-01", "2020-04-15", "2020-04-30",
                          "2020-05-15", "2020-06-15", "2020-07-10", "2020-08-15",
                          "2020-09-09"))
tray <- do.call(rbind, lapply(colnames(Rt_med), function(nm) {
  v <- Rt_med[, nm]; bajo <- which(v < 1)
  dif <- v - Rt_med[, "M0_AR1_Gamma"]
  cbind(data.frame(modelo = nm, Rt_max = max(v), fecha_max = format(fechas[which.max(v)]),
                   primer_dia_bajo_1 = if (length(bajo)) format(fechas[bajo[1]]) else NA,
                   dif_abs_media_vs_M0 = mean(abs(dif)), dif_abs_max_vs_M0 = max(abs(dif)),
                   dias_dif_mayor_02_vs_M0 = sum(abs(dif) > 0.2)),
        t(setNames(v[match(FECHAS_CLAVE, fechas)], paste0("Rt_", format(FECHAS_CLAVE, "%d%b")))))
}))
cat("\n--- 6c. Trayectorias de R_t (medianas posteriores) ---\n"); print(tray, digits = 3, row.names = FALSE)
write.csv(tray, file.path(OUT, "f2cg_trayectorias_resumen.csv"), row.names = FALSE)
write.csv(data.frame(fecha = fechas, Rt_med), file.path(OUT, "f2cg_Rt_medianas_diarias.csv"),
          row.names = FALSE)

colores <- c(M0_AR1_Gamma = "black", M1_AR2_Gamma = "#377EB8",
             M2_RW_Gamma = "#E41A1C", M3_AR1_Exp = "#4DAF4A")
dfg <- do.call(rbind, lapply(colnames(Rt_med), function(nm) {
  M <- fits[[nm]]$draws("Rt", format = "matrix")
  data.frame(fecha = fechas, modelo = nm, med = apply(M, 2, median),
             lo = apply(M, 2, quantile, 0.025), hi = apply(M, 2, quantile, 0.975))
}))
p <- ggplot(dfg, aes(fecha, med, colour = modelo, fill = modelo)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.7) + geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = colores) + scale_fill_manual(values = colores) +
  labs(x = NULL, y = expression(R[t]), colour = NULL, fill = NULL,
       title = "Fase 2 (IS China-Gamma): trayectorias de Rt, mediana e IC95") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "f2cg_trayectorias.png"), p, width = 10, height = 5.5, dpi = 150)
ggsave(file.path(FIG, "f2cg_trayectorias.pdf"), p, width = 10, height = 5.5)
cat("  Guardado: f2cg_trayectorias.png/.pdf\n")

# --- 6d. Largo plazo del AR(1) (M0 y M3) ------------------------------------------
largo <- do.call(rbind, lapply(c("M0_AR1_Gamma", "M3_AR1_Exp"), function(nm) {
  d <- fits[[nm]]$draws(c("mu_rt", "rho1", "sigma_epsilon"), format = "df")
  s2 <- d$sigma_epsilon^2 / (1 - d$rho1^2)
  cuant <- function(v) sprintf("%.3f [%.3f, %.3f]", median(v),
                               quantile(v, 0.025), quantile(v, 0.975))
  data.frame(modelo = nm,
             rho1 = cuant(d$rho1),
             P_rho1_mayor_09 = mean(d$rho1 > 0.9),
             mediana_estacionaria = cuant(exp(d$mu_rt)),
             media_estacionaria = cuant(exp(d$mu_rt + s2 / 2)),
             vida_media_dias = cuant(log(0.5) / log(d$rho1)))
}))
cat("\n--- 6d. Parametros de largo plazo (mediana [IC95]) ---\n"); print(largo, row.names = FALSE)
write.csv(largo, file.path(OUT, "f2cg_largo_plazo.csv"), row.names = FALSE)


# =============================================================================
#  PASO 7 - DEFINITIVO VIEJO (10.5 d) vs NUEVO (China-Gamma)   ~2 min, sin Stan
#
#  Base del paso 3 del plan: que se recicla y que cambia en el documento.
#  Viejo: resultados/corrida1_ar1_ad099.rds y corrida3_rw.rds (script 03).
# =============================================================================

cat("\n=== PASO 7: DEFINITIVO VIEJO vs NUEVO ===\n")

old_ar1 <- readRDS(file.path(OUT_OLD, "corrida1_ar1_ad099.rds"))
old_rw  <- readRDS(file.path(OUT_OLD, "corrida3_rw.rds"))
pares <- list(AR1 = list(viejo = old_ar1, nuevo = fits$M0_AR1_Gamma, pars = PARS_AR1),
              RW  = list(viejo = old_rw,  nuevo = fits$M2_RW_Gamma,  pars = PARS_RW))

vn <- bind_rows(lapply(names(pares), function(nm) {   # bind_rows, no rbind: las filas del RW tienen otras columnas
  pr <- pares[[nm]]
  fila <- function(fit, etiqueta) {
    s  <- fit$summary(pr$pars); dg <- fit$diagnostic_summary(quiet = TRUE)
    ll <- fit$draws("log_lik", format = "matrix")
    lo <- loo(ll, r_eff = relative_eff(exp(ll), chain_id = cid))
    Rt <- apply(fit$draws("Rt", format = "matrix"), 2, median); bajo <- which(Rt < 1)
    cbind(data.frame(proceso = nm, version = etiqueta,
                     divergencias = sum(dg$num_divergent), rhat_max = max(s$rhat),
                     ess_bulk_min = min(s$ess_bulk), ebfmi_min = min(dg$ebfmi),
                     elpd = lo$estimates["elpd_loo", "Estimate"],
                     k_max = max(lo$diagnostics$pareto_k),
                     Rt_max = max(Rt), fecha_max = format(fechas[which.max(Rt)]),
                     primer_dia_bajo_1 = if (length(bajo)) format(fechas[bajo[1]]) else NA),
          t(setNames(s$mean, paste0("media_", s$variable))),
          t(setNames(Rt[match(FECHAS_CLAVE, fechas)], paste0("Rt_", format(FECHAS_CLAVE, "%d%b")))))
  }
  bind_rows(fila(pr$viejo, "viejo_IS_10.5d"), fila(pr$nuevo, "nuevo_China-Gamma"))
}))
print(vn, digits = 3, row.names = FALSE)
write.csv(vn, file.path(OUT, "f2cg_viejo_vs_nuevo.csv"), row.names = FALSE)

Rt_viejo <- apply(old_ar1$draws("Rt", format = "matrix"), 2, median)
Rt_nuevo <- Rt_med[, "M0_AR1_Gamma"]
crec <- Rt_nuevo > 1
cat(sprintf("\n  AR(1): razon media Rt viejo/nuevo cuando Rt > 1: %.3f | cuando Rt < 1: %.3f\n",
            mean(Rt_viejo[crec] / Rt_nuevo[crec]), mean(Rt_viejo[!crec] / Rt_nuevo[!crec])))
cat(sprintf("  AR(1): dias con |Rt viejo - Rt nuevo| > 0.2: %d de %d\n",
            sum(abs(Rt_viejo - Rt_nuevo) > 0.2), T_PERIODO))
cat("  Guardado: f2cg_viejo_vs_nuevo.csv\n")

cat("\n=== SCRIPT 06 COMPLETO ===\n")
cat("Pega en el chat la salida de los PASOS 1, 6 y 7 (y los diagnosticos de 2-5).\n")
