# ==============================================================================
#  ECUACIÓN DE RENOVACIÓN BAYESIANA — COVID-19 BOGOTÁ
#  Exógeno: Gamma(alpha=1.401, beta=0.168) — calibrado por MoM
#  14-Mar-2020 → 09-Sep-2020 (180 días)
#  5 IS × 2 modelos (NB-2 + Poisson) = 10 runs
#
#  CAMBIO RESPECTO A analisis_bogota_ar2.R:
#    - Exógeno Rayleigh(sigma_mu) → Gamma(alpha_mu, beta_mu)
#    - alpha_mu=1.401, beta_mu=0.168 calibrados por MoM sobre datos reales
#    - sigma_mu eliminado del espacio de parámetros
#    - Hiperparámetros pasados como datos Stan (flexibilidad sin recompilar)
#
#  CALIBRACIÓN MoM (datos importados 14-25 mar 2020, n=12):
#    media_obs = 8.33, var_obs = 49.52
#    alpha = media^2 / var = 69.39 / 49.52 = 1.401
#    beta  = media   / var =  8.33 / 49.52 = 0.168
#    Moda  = (alpha-1)/beta = 0.401/0.168  = 2.39 casos/día (> 0)
#    Media implicada = alpha/beta = 1.401/0.168 = 8.33 ✓
# ==============================================================================

# ==============================================================================
# SECCIÓN 1: CONFIGURACIÓN
# ==============================================================================

# Rtools (solo Windows — ajustar versión si es necesario)
if (.Platform$OS.type == "windows") {
  rtools_path <- "C:/rtools44/usr/bin;C:/rtools44/mingw64/bin"
  Sys.setenv(PATH = paste(rtools_path, Sys.getenv("PATH"), sep = ";"))
}

# --- Rutas relativas a la raíz del repositorio ---
library(here)
RUTA_DATOS  <- here("MODELO_STAN_4", "datos_agregados.xlsx")
RUTA_STAN   <- here("MODELO_STAN_4", "EXOGENO_GAMMA_ANALISIS_PRELIMINAR")
RUTA_OUTPUT <- RUTA_STAN   # mismo directorio

# Creamos carpeta de output (si no existe)
if (!dir.exists(RUTA_OUTPUT)) dir.create(RUTA_OUTPUT, recursive = TRUE)
setwd(RUTA_OUTPUT)

# --- Control de ejecución ---
TEST_RUN      <- FALSE   # TRUE = prueba rápida (2 cadenas, 500 iter, solo IS ref)
N_CHAINS      <- 4
N_WARMUP      <- 1000
N_SAMPLING    <- 1000    # iteraciones de muestreo por cadena (total draws = 4000)
ADAPT_DELTA   <- 0.95
MAX_TREEDEPTH <- 12
#ADAPT_DELTA   <- 0.98
#MAX_TREEDEPTH <- 15

# --- Hiperparámetros exógeno Gamma (MoM) ---
ALPHA_MU <- 1.401   # shape
BETA_MU  <- 0.168   # rate

cat("=== CONFIGURACIÓN ===\n")
cat(sprintf("Output:        %s\n", RUTA_OUTPUT))
cat(sprintf("Cadenas:       %d × %d warmup + %d sampling\n",
            N_CHAINS, N_WARMUP, N_SAMPLING))
cat(sprintf("Exógeno:       Gamma(alpha=%.3f, beta=%.3f)\n", ALPHA_MU, BETA_MU))
cat(sprintf("  Media:       %.3f casos/día\n", ALPHA_MU / BETA_MU))
cat(sprintf("  Moda:        %.3f casos/día\n", (ALPHA_MU - 1) / BETA_MU))
cat(sprintf("  Varianza:    %.3f\n", ALPHA_MU / BETA_MU^2))

# --- Paquetes ---
suppressPackageStartupMessages({
  library(readxl)
  library(loo)
  library(posterior)
  library(bayesplot)
  library(gridExtra)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
})

USE_CMDSTANR <- requireNamespace("cmdstanr", quietly = TRUE)
if (USE_CMDSTANR) {
  library(cmdstanr)
  cat("Backend: cmdstanr\n")
} else {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  cat("Backend: rstan\n")
}

if (TEST_RUN) {
  N_CHAINS <- 2; N_WARMUP <- 250; N_SAMPLING <- 250
  cat("MODO TEST: 2 cadenas × 500 iter\n")
}

set.seed(42)

# ==============================================================================
# SECCIÓN 2: DATOS
# ==============================================================================

cat("\n=== CARGANDO DATOS ===\n")

df_raw        <- as.data.frame(read_excel(RUTA_DATOS, sheet = "NO_IMPORTADOS"))
names(df_raw) <- c("fecha", "casos")
df_raw        <- df_raw[!grepl("nan", as.character(df_raw$fecha),
                                ignore.case = TRUE), ]
df_raw        <- df_raw[!is.na(df_raw$fecha), ]
df_raw$fecha  <- as.Date(as.numeric(df_raw$fecha), origin = "1899-12-30")
df_raw$casos  <- as.integer(df_raw$casos)
df_local      <- df_raw[order(df_raw$fecha), ]

FECHA_INICIO <- as.Date("2020-03-14")
FECHA_FIN    <- FECHA_INICIO + 179L
T_PERIODO    <- 180L
N_EXOGENO    <- 12L   # días pre-cuarentena (14-25 mar)

grilla <- seq(FECHA_INICIO, FECHA_FIN, by = "day")
df_periodo        <- merge(data.frame(fecha = grilla),
                           df_local[df_local$fecha <= FECHA_FIN, ],
                           by = "fecha", all.x = TRUE)
df_periodo        <- df_periodo[order(df_periodo$fecha), ]
df_periodo$casos[is.na(df_periodo$casos)] <- 0L
y_obs             <- as.integer(df_periodo$casos)

cat(sprintf("Período:       %s → %s (%d días)\n",
            FECHA_INICIO, FECHA_FIN, T_PERIODO))
cat(sprintf("Casos totales: %s | Media: %.1f/día | Max: %d\n",
            format(sum(y_obs), big.mark=","), mean(y_obs), max(y_obs)))
cat(sprintf("Var/Media:     %.2f (sobredispersión)\n",
            var(y_obs) / mean(y_obs)))
cat(sprintf("Días exógenos: t=1..%d (hasta 25-mar-2020)\n", N_EXOGENO))

# ==============================================================================
# SECCIÓN 3: INTERVALOS SERIALES
# ==============================================================================

cat("\n=== INTERVALOS SERIALES ===\n")

MAX_SI <- 30L

discret_si <- function(pfun, ...) {
  g <- numeric(MAX_SI)
  g[1] <- pfun(1.5, ...) - pfun(0, ...)
  for (s in 2:MAX_SI) g[s] <- pfun(s+0.5,...) - pfun(s-0.5,...)
  g[g < 0] <- 0; g / sum(g)
}

discret_si_normal_trunc <- function(mu_si, sd_si) {
  p_pos <- pnorm(0, mu_si, sd_si, lower.tail = FALSE)
  g <- numeric(MAX_SI)
  g[1] <- (pnorm(1.5,mu_si,sd_si) - pnorm(0,mu_si,sd_si)) / p_pos
  for (s in 2:MAX_SI)
    g[s] <- (pnorm(s+0.5,mu_si,sd_si) - pnorm(s-0.5,mu_si,sd_si)) / p_pos
  g[g < 0] <- 0; g / sum(g)
}

si_lista <- list(
  referencia  = discret_si(pgamma, shape=6.5,  rate=0.62),
  china_norm  = discret_si_normal_trunc(5.0, 5.2),
  china_gamma = discret_si(pgamma, shape=2.29, rate=0.36),
  colombia    = discret_si(pgamma, shape=1.96, rate=0.51),
  burkina     = discret_si(pgamma, shape=1.04, rate=0.18)
)

for (nm in names(si_lista))
  cat(sprintf("  %-14s media = %.2f días\n", nm,
              sum(1:MAX_SI * si_lista[[nm]])))

IS_a_usar <- if (TEST_RUN) si_lista["referencia"] else si_lista

# ==============================================================================
# SECCIÓN 4: FUNCIÓN AUXILIAR DATOS STAN
# ==============================================================================

make_stan_data <- function(g_vec) {
  list(
    T          = T_PERIODO,
    n_exogeno  = N_EXOGENO,
    max_si     = MAX_SI,
    y          = y_obs,
    g          = g_vec,
    alpha_mu   = ALPHA_MU,
    beta_mu    = BETA_MU
  )
}

# ==============================================================================
# SECCIÓN 5: FUNCIONES AUXILIARES EXTRACCIÓN
# ==============================================================================

extr_mat <- function(fit, param) {
  if (inherits(fit, "CmdStanMCMC")) fit$draws(param, format="matrix")
  else {
    ex <- rstan::extract(fit, pars=param)[[param]]
    if (is.matrix(ex)) ex else matrix(ex, ncol=1)
  }
}

extr_ll <- function(fit) {
  if (inherits(fit, "CmdStanMCMC")) fit$draws("log_lik", format="matrix")
  else {
    ll <- loo::extract_log_lik(fit, merge_chains=FALSE)
    n_i <- dim(ll)[1]; n_c <- dim(ll)[2]; n_o <- dim(ll)[3]
    t(matrix(aperm(ll, c(3,1,2)), nrow=n_o, ncol=n_i*n_c))
  }
}

get_rhat_max <- function(fit, es_nb2 = TRUE) {
  pars <- c("mu_rt","rho1","rho2","sigma_epsilon")
  if (es_nb2) pars <- c(pars, "phi")
  tryCatch({
    if (USE_CMDSTANR) max(fit$summary(pars)$rhat, na.rm=TRUE)
    else max(summary(fit, pars=pars)$summary[,"Rhat"], na.rm=TRUE)
  }, error=function(e) NA_real_)
}

get_ess_min <- function(fit, es_nb2 = TRUE) {
  pars <- c("mu_rt","rho1","rho2","sigma_epsilon")
  if (es_nb2) pars <- c(pars, "phi")
  tryCatch({
    if (USE_CMDSTANR) min(fit$summary(pars)$ess_bulk, na.rm=TRUE)
    else min(summary(fit, pars=pars)$summary[,"n_eff"], na.rm=TRUE)
  }, error=function(e) NA_real_)
}

get_divergencias <- function(fit) {
  tryCatch({
    if (USE_CMDSTANR) {
      diags <- fit$diagnostic_summary(quiet=TRUE)
      sum(diags$num_divergent)
    } else {
      sp <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
      sum(sapply(sp, function(x) sum(x[,"divergent__"])))
    }
  }, error=function(e) NA_integer_)
}

# ==============================================================================
# SECCIÓN 6: AJUSTE DE MODELOS
# ==============================================================================

cat("\n=== AJUSTE DE MODELOS (Exógeno Gamma) ===\n")
cat(sprintf("5 IS × 2 modelos = 10 runs | %d cadenas × %d sampling\n",
            N_CHAINS, N_SAMPLING))

modelos_stan <- list(
  NB2     = file.path(RUTA_STAN, "renewal_bogota_nb2_ar2_gamma.stan"),
  Poisson = file.path(RUTA_STAN, "renewal_bogota_poisson_ar2_gamma.stan")
)

# Verificamos que los .stan existen
for (nm in names(modelos_stan)) {
  if (!file.exists(modelos_stan[[nm]]))
    stop(sprintf("No encontrado: %s\nAsegúrate de copiar los .stan a la carpeta destino.",
                 modelos_stan[[nm]]))
}

ajustar_stan <- function(stan_file, stan_data, tag) {
  cat(sprintf("  [%s] Ajustando...\n", tag))
  es_poisson <- grepl("poisson", tolower(stan_file))

  hacer_init <- function() {
    ini <- list(
      mu_rt         = 0.0,
      rho1          = 0.7,
      rho2          = 0.1,
      sigma_epsilon = 0.1,
      epsilon_raw   = rep(0.0, stan_data$T),
      mu_exo        = rep(ALPHA_MU / BETA_MU, stan_data$n_exogeno)  # media Gamma
    )
    if (!es_poisson) ini$phi <- 4.0
    ini
  }

  if (USE_CMDSTANR) {
    mod <- cmdstan_model(stan_file)
    mod$sample(
      data            = stan_data,
      iter_warmup     = N_WARMUP,
      iter_sampling   = N_SAMPLING,
      chains          = N_CHAINS,
      parallel_chains = min(N_CHAINS, parallel::detectCores()),
      adapt_delta     = ADAPT_DELTA,
      max_treedepth   = MAX_TREEDEPTH,
      init            = hacer_init,
      seed            = 42,
      refresh         = 200,
      show_messages   = FALSE
    )
  } else {
    inits_lista <- replicate(N_CHAINS, hacer_init(), simplify=FALSE)
    rstan::stan(
      file    = stan_file,
      data    = stan_data,
      iter    = N_WARMUP + N_SAMPLING,
      warmup  = N_WARMUP,
      chains  = N_CHAINS,
      init    = inits_lista,
      control = list(adapt_delta=ADAPT_DELTA, max_treedepth=MAX_TREEDEPTH),
      seed    = 42,
      refresh = 200
    )
  }
}

resultados <- list()
t_inicio_total <- proc.time()

for (nm_si in names(IS_a_usar)) {
  resultados[[nm_si]] <- list()
  for (nm_mod in names(modelos_stan)) {
    tag      <- paste0(nm_si, "_", nm_mod)
    rds_file <- file.path(RUTA_OUTPUT, paste0("fit_gamma_", tag, ".rds"))

    if (file.exists(rds_file)) {
      cat(sprintf("  [%s] Cargando existente: %s\n", tag, basename(rds_file)))
      resultados[[nm_si]][[nm_mod]] <- readRDS(rds_file)
      next
    }

    sd_i   <- make_stan_data(IS_a_usar[[nm_si]])
    t_ini  <- proc.time()

    tryCatch({
      fit <- ajustar_stan(modelos_stan[[nm_mod]], sd_i, tag)
      resultados[[nm_si]][[nm_mod]] <- fit
      saveRDS(fit, rds_file)
      t_elapsed <- (proc.time() - t_ini)[["elapsed"]]
      cat(sprintf("    ✓ Guardado: %s (%.1f min)\n",
                  basename(rds_file), t_elapsed/60))
    }, error = function(e) {
      cat(sprintf("    ✗ ERROR [%s]: %s\n", tag, conditionMessage(e)))
    })
  }
}

t_total <- (proc.time() - t_inicio_total)[["elapsed"]]
cat(sprintf("\nTiempo total ajuste: %.1f min\n", t_total/60))

# ==============================================================================
# SECCIÓN 7: TABLA DE CONVERGENCIA (diagnósticos rápidos)
# ==============================================================================

cat("\n=== DIAGNÓSTICOS DE CONVERGENCIA ===\n")
cat(sprintf("%-35s %8s %10s %12s %6s\n",
            "Modelo", "Rhat_max", "ESS_min", "Divergencias", "Estado"))
cat(strrep("-", 75), "\n")

tabla_conv <- do.call(rbind, lapply(names(IS_a_usar), function(nm_si) {
  do.call(rbind, lapply(names(modelos_stan), function(nm_mod) {
    fit <- resultados[[nm_si]][[nm_mod]]
    if (is.null(fit)) return(NULL)
    tag    <- paste0(nm_si, "_", nm_mod)
    es_nb2 <- nm_mod == "NB2"
    rhat   <- get_rhat_max(fit, es_nb2)
    ess    <- get_ess_min(fit, es_nb2)
    ndiv   <- get_divergencias(fit)
    ok_rhat <- !is.na(rhat) && rhat < 1.01
    ok_ess  <- !is.na(ess)  && ess  > 400
    ok_div  <- !is.na(ndiv) && ndiv == 0
    estado  <- if (ok_rhat && ok_ess && ok_div) "✓ OK" else "⚠ REVISAR"
    cat(sprintf("  %-33s %8.4f %10.0f %12d %6s\n",
                tag, rhat, ess, ndiv, estado))
    data.frame(IS=nm_si, Modelo=nm_mod, Rhat_max=rhat,
               ESS_min=ess, Divergencias=ndiv, Estado=estado,
               stringsAsFactors=FALSE)
  }))
}))

write.csv(tabla_conv,
          file.path(RUTA_OUTPUT, "diagnosticos_convergencia_gamma.csv"),
          row.names=FALSE)
cat("\n✓ diagnosticos_convergencia_gamma.csv\n")

# ==============================================================================
# SECCIÓN 8: DRAWS ESENCIALES — Guardamos ligero para análisis posteriores
# ==============================================================================
# Guardamos solo lo necesario para no tener que repetir el MCMC:
#   - Parámetros escalares: mu_rt, rho1, rho2, sigma_epsilon, phi (si NB2)
#   - Vectores clave: Rt[T], log_lik[T]  (necesarios para LOO-CV y sensibilidad)
#   - NO guardamos: f[T], mu[T], y_rep[T], epsilon[T] (recalculables o pesados)
# ==============================================================================

cat("\n=== GUARDANDO DRAWS ESENCIALES ===\n")

for (nm_si in names(IS_a_usar)) {
  for (nm_mod in names(modelos_stan)) {
    fit <- resultados[[nm_si]][[nm_mod]]
    if (is.null(fit)) next
    tag       <- paste0(nm_si, "_", nm_mod)
    draws_file <- file.path(RUTA_OUTPUT,
                            paste0("draws_gamma_", tag, ".rds"))
    if (file.exists(draws_file)) {
      cat(sprintf("  [%s] draws ya existen\n", tag))
      next
    }
    tryCatch({
      es_nb2 <- nm_mod == "NB2"
      pars_escalares <- c("mu_rt","rho1","rho2","sigma_epsilon")
      if (es_nb2) pars_escalares <- c(pars_escalares, "phi")

      draws_list <- list(
        escalares = lapply(pars_escalares, function(p)
          as.vector(extr_mat(fit, p))),
        Rt        = extr_mat(fit, "Rt"),
        log_lik   = extr_ll(fit),
        log_Rt    = extr_mat(fit, "log_Rt"),
        y_pred    = extr_mat(fit, "y_pred"),
        meta      = list(IS=nm_si, modelo=nm_mod,
                         n_draws = N_CHAINS * N_SAMPLING,
                         alpha_mu=ALPHA_MU, beta_mu=BETA_MU,
                         fecha_ajuste=Sys.time())
      )
      names(draws_list$escalares) <- pars_escalares
      saveRDS(draws_list, draws_file)
      sz <- round(file.size(draws_file)/1e6, 1)
      cat(sprintf("  ✓ [%s] draws guardados (%.1f MB)\n", tag, sz))
    }, error=function(e)
      cat(sprintf("  ✗ ERROR draws [%s]: %s\n", tag, conditionMessage(e))))
  }
}

# ==============================================================================
# SECCIÓN 9: LOO-CV (NB-2 vs Poisson, IS referencia)
# ==============================================================================

cat("\n=== LOO-CV: NB-2 vs Poisson (IS: referencia) ===\n")

loo_res <- list()
for (nm_mod in names(modelos_stan)) {
  fit <- resultados[["referencia"]][[nm_mod]]
  if (is.null(fit)) next
  ll <- tryCatch(extr_ll(fit), error=function(e) NULL)
  if (is.null(ll)) next
  #r_eff <- relative_eff(exp(ll))
  r_eff <- relative_eff(exp(ll), 
                        chain_id = rep(1:N_CHAINS, each = N_SAMPLING))
  loo_res[[nm_mod]] <- tryCatch(loo(ll, r_eff=r_eff), error=function(e) NULL)
  if (!is.null(loo_res[[nm_mod]]))
    cat(sprintf("  %-8s elpd_loo = %.1f (SE=%.1f)\n", nm_mod,
                loo_res[[nm_mod]]$estimates["elpd_loo","Estimate"],
                loo_res[[nm_mod]]$estimates["elpd_loo","SE"]))
}

if (length(loo_res) == 2) {
  comp <- loo_compare(loo_res)
  cat("\nComparación LOO-CV:\n"); print(comp)
  write.csv(as.data.frame(comp),
            file.path(RUTA_OUTPUT, "loo_comparison_gamma.csv"))
  delta <- comp[2,"elpd_diff"]
  se    <- comp[2,"se_diff"]
  cat(sprintf("\nNB-2 vs Poisson: Δelpd=%.1f (%.1f SE) → %s\n",
              abs(delta), abs(delta/se),
              if (abs(delta/se)>2) "NB-2 SUPERIOR" else "diferencia moderada"))
}

###### Para revisión de por qué Poisson pareciera tener mejor resultado
# Ver diagnóstico k de ambos modelos
print(loo_res[["NB2"]])
print(loo_res[["Poisson"]])

# Cuántas observaciones tienen k > 0.7
cat("NB2 - k > 0.7:", 
    sum(loo_res[["NB2"]]$diagnostics$pareto_k > 0.7), "\n")
cat("Poisson - k > 0.7:", 
    sum(loo_res[["Poisson"]]$diagnostics$pareto_k > 0.7), "\n")

# ==============================================================================
# SECCIÓN 10: POSTERIOR PARÁMETROS CLAVE (IS referencia, NB-2)
# ==============================================================================

cat("\n=== POSTERIOR PARÁMETROS (IS: referencia, NB-2) ===\n")

fit_ref <- resultados[["referencia"]][["NB2"]]
if (!is.null(fit_ref)) {
  fmt <- function(v)
    sprintf("%.3f [%.3f, %.3f]", mean(v), quantile(v,.025), quantile(v,.975))

  phi_v  <- as.vector(extr_mat(fit_ref, "phi"))
  rho1_v <- as.vector(extr_mat(fit_ref, "rho1"))
  rho2_v <- as.vector(extr_mat(fit_ref, "rho2"))
  sig_v  <- as.vector(extr_mat(fit_ref, "sigma_epsilon"))
  mu_rt_v <- as.vector(extr_mat(fit_ref, "mu_rt"))

  cat(sprintf("  phi           = %s\n", fmt(phi_v)))
  cat(sprintf("  rho1          = %s  (prior: N(0.7,0.15))\n", fmt(rho1_v)))
  cat(sprintf("  rho2          = %s  (prior: N(0.1,0.10))\n", fmt(rho2_v)))
  cat(sprintf("  sigma_epsilon = %s\n", fmt(sig_v)))
  cat(sprintf("  mu_rt         = %s\n", fmt(mu_rt_v)))
  cat(sprintf("  Rt medio      = %.3f\n", mean(exp(mu_rt_v))))

  # Guardar tabla posterior
  post_df <- data.frame(
    parametro = c("phi","rho1","rho2","sigma_epsilon","mu_rt"),
    media     = c(mean(phi_v), mean(rho1_v), mean(rho2_v),
                  mean(sig_v), mean(mu_rt_v)),
    q025      = c(quantile(phi_v,.025), quantile(rho1_v,.025),
                  quantile(rho2_v,.025), quantile(sig_v,.025),
                  quantile(mu_rt_v,.025)),
    q975      = c(quantile(phi_v,.975), quantile(rho1_v,.975),
                  quantile(rho2_v,.975), quantile(sig_v,.975),
                  quantile(mu_rt_v,.975))
  )
  write.csv(post_df,
            file.path(RUTA_OUTPUT, "posterior_params_gamma.csv"),
            row.names=FALSE)
  cat("✓ posterior_params_gamma.csv\n")
}

# ==============================================================================
# SECCIÓN 11: SENSIBILIDAD AL IS
# ==============================================================================

cat("\n=== SENSIBILIDAD AL IS (NB-2) ===\n")

sens_df <- do.call(rbind, lapply(names(IS_a_usar), function(nm_si) {
  fit <- resultados[[nm_si]][["NB2"]]
  if (is.null(fit)) return(NULL)
  Rt_m <- extr_mat(fit, "Rt")
  data.frame(fecha  = grilla,
             media  = colMeans(Rt_m),
             q025   = apply(Rt_m, 2, quantile, .025),
             q975   = apply(Rt_m, 2, quantile, .975),
             IS     = nm_si,
             stringsAsFactors=FALSE)
}))

if (!is.null(sens_df) && nrow(sens_df) > 0) {
  rng <- tapply(sens_df$media, sens_df$fecha, function(x) max(x) - min(x))
  cat(sprintf("  Rango medio Rt entre IS:  %.3f\n", mean(rng)))
  cat(sprintf("  Días con rango > 0.2:     %d de %d\n",
              sum(rng > 0.2), length(rng)))
  write.csv(sens_df,
            file.path(RUTA_OUTPUT, "sensibilidad_IS_gamma.csv"),
            row.names=FALSE)
  cat("✓ sensibilidad_IS_gamma.csv\n")
}

# ==============================================================================
# SECCIÓN 12: TABLA RESUMEN FINAL
# ==============================================================================

cat("\n=== TABLA RESUMEN FINAL ===\n")

tabla_final <- do.call(rbind, lapply(names(IS_a_usar), function(nm_si) {
  do.call(rbind, lapply(names(modelos_stan), function(nm_mod) {
    fit <- resultados[[nm_si]][[nm_mod]]
    if (is.null(fit)) return(NULL)
    es_nb2 <- nm_mod == "NB2"

    fmt3 <- function(p) tryCatch({
      v <- as.vector(extr_mat(fit, p))
      sprintf("%.3f [%.3f,%.3f]", mean(v), quantile(v,.025), quantile(v,.975))
    }, error=function(e) NA_character_)

    ll <- tryCatch(extr_ll(fit), error=function(e) NULL)
    loo_val <- NA_real_; loo_se <- NA_real_
    if (!is.null(ll)) {
      r_eff   <- tryCatch(relative_eff(exp(ll)), error=function(e) NULL)
      loo_obj <- tryCatch(loo(ll, r_eff=r_eff), error=function(e) NULL)
      if (!is.null(loo_obj)) {
        loo_val <- round(loo_obj$estimates["elpd_loo","Estimate"], 1)
        loo_se  <- round(loo_obj$estimates["elpd_loo","SE"], 1)
      }
    }

    data.frame(
      IS           = nm_si,
      Modelo       = nm_mod,
      elpd_loo     = loo_val,
      SE_loo       = loo_se,
      Rhat_max     = round(get_rhat_max(fit, es_nb2), 4),
      ESS_min      = round(get_ess_min(fit, es_nb2), 0),
      Divergencias = get_divergencias(fit),
      phi          = if (es_nb2) fmt3("phi") else NA_character_,
      rho1         = fmt3("rho1"),
      rho2         = fmt3("rho2"),
      sigma_eps    = fmt3("sigma_epsilon"),
      stringsAsFactors=FALSE
    )
  }))
}))

print(tabla_final, row.names=FALSE)
write.csv(tabla_final,
          file.path(RUTA_OUTPUT, "tabla_resumen_gamma.csv"),
          row.names=FALSE)

# ==============================================================================
# SECCIÓN 13: FIGURAS
# ==============================================================================

cat("\n=== GENERANDO FIGURAS ===\n")

COLS <- c(referencia="#E41A1C", china_norm="#377EB8",
          china_gamma="#4DAF4A", colombia="#984EA3", burkina="#FF7F00")

EVENTOS <- data.frame(
  fecha = as.Date(c("2020-03-25","2020-06-01","2020-09-01")),
  label = c("Cuarentena","Reapertura\npiloto","Reapertura\ngradual"),
  stringsAsFactors=FALSE)

# --- Fig 1: traceplots parámetros clave (modelo referencia NB2) ---
fit_ref <- resultados[["referencia"]][["NB2"]]
if (!is.null(fit_ref)) {
  tryCatch({
    if (USE_CMDSTANR) {
      draws_arr <- fit_ref$draws(
        c("mu_rt","rho1","rho2","sigma_epsilon","phi"),
        format="array")
    } else {
      draws_arr <- as.array(fit_ref,
        pars=c("mu_rt","rho1","rho2","sigma_epsilon","phi"))
    }
    p_trace <- mcmc_trace(draws_arr,
      pars=c("mu_rt","rho1","rho2","sigma_epsilon","phi")) +
      labs(title="Traceplots — IS referencia, NB-2, Exógeno Gamma") +
      theme_minimal(base_size=10)
    ggsave(file.path(RUTA_OUTPUT, "fig1_traceplots_gamma.pdf"),
           p_trace, width=14, height=10)
    cat("  ✓ fig1_traceplots_gamma.pdf\n")
  }, error=function(e) cat("  ✗ traceplots:", conditionMessage(e), "\n"))
}

# --- Fig 2: Rhat y ESS por modelo ---
if (!is.null(tabla_conv) && nrow(tabla_conv) > 0) {
  tryCatch({
    tc_plot <- tabla_conv
    tc_plot$tag <- paste0(tc_plot$IS, "\n", tc_plot$Modelo)

    p_rhat <- ggplot(tc_plot, aes(x=tag, y=Rhat_max, fill=Estado)) +
      geom_col(width=0.6) +
      geom_hline(yintercept=1.01, linetype="dashed", color="red") +
      scale_fill_manual(values=c("✓ OK"="#4DAF4A","⚠ REVISAR"="#E41A1C")) +
      labs(title="Rhat máximo por modelo (umbral=1.01)",
           x=NULL, y="Rhat_max") +
      theme_minimal(base_size=9) +
      theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none")

    p_ess <- ggplot(tc_plot, aes(x=tag, y=ESS_min, fill=Estado)) +
      geom_col(width=0.6) +
      geom_hline(yintercept=400, linetype="dashed", color="red") +
      scale_fill_manual(values=c("✓ OK"="#4DAF4A","⚠ REVISAR"="#E41A1C")) +
      labs(title="ESS mínimo por modelo (umbral=400)",
           x=NULL, y="ESS_min") +
      theme_minimal(base_size=9) +
      theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none")

    pdf(file.path(RUTA_OUTPUT, "fig2_diagnosticos_gamma.pdf"),
        width=14, height=6)
    grid.arrange(p_rhat, p_ess, ncol=2,
      top="Diagnósticos de convergencia — Exógeno Gamma")
    dev.off()
    cat("  ✓ fig2_diagnosticos_gamma.pdf\n")
  }, error=function(e) cat("  ✗ fig diagnósticos:", conditionMessage(e), "\n"))
}

# --- Fig 3: Ajuste principal (referencia NB2) ---
fit_ref <- resultados[["referencia"]][["NB2"]]
if (!is.null(fit_ref)) {
  tryCatch({
    f_m  <- extr_mat(fit_ref, "f")
    Rt_m <- extr_mat(fit_ref, "Rt")
    mu_m <- extr_mat(fit_ref, "mu")

    rf  <- data.frame(fecha=grilla, obs=y_obs,
                      pred=colMeans(f_m),
                      q025=apply(f_m,2,quantile,.025),
                      q975=apply(f_m,2,quantile,.975))
    rRt <- data.frame(fecha=grilla,
                      media=colMeans(Rt_m),
                      q025=apply(Rt_m,2,quantile,.025),
                      q975=apply(Rt_m,2,quantile,.975))
    rmu <- data.frame(fecha=grilla[1:N_EXOGENO],
                      media=colMeans(mu_m[,1:N_EXOGENO,drop=FALSE]),
                      q025=apply(mu_m[,1:N_EXOGENO,drop=FALSE],
                                 2,quantile,.025),
                      q975=apply(mu_m[,1:N_EXOGENO,drop=FALSE],
                                 2,quantile,.975))

    pa <- ggplot(rf, aes(fecha)) +
      geom_col(aes(y=obs), fill="salmon", alpha=.7, width=1) +
      geom_ribbon(aes(ymin=q025,ymax=q975), fill="#2166AC", alpha=.25) +
      geom_line(aes(y=pred), color="#2166AC", linewidth=.9) +
      geom_vline(data=EVENTOS, aes(xintercept=fecha),
                 linetype="dashed", color="gray40") +
      scale_y_continuous(labels=comma) +
      scale_x_date(date_breaks="3 weeks", date_labels="%d\n%b") +
      labs(title="Casos observados vs predichos — AR(2), NB-2, Exógeno Gamma",
           x=NULL, y="Casos/día") + theme_minimal(base_size=12)

    pb <- ggplot(rRt, aes(fecha)) +
      geom_hline(yintercept=1, linetype="dotted", linewidth=.8) +
      geom_ribbon(aes(ymin=q025,ymax=q975), fill="#1B7837", alpha=.25) +
      geom_line(aes(y=media), color="#1B7837", linewidth=.9) +
      geom_vline(data=EVENTOS, aes(xintercept=fecha),
                 linetype="dashed", color="gray40") +
      scale_x_date(date_breaks="3 weeks", date_labels="%d\n%b") +
      scale_y_log10(breaks=c(.1,.5,1,2,5,10)) +
      labs(title="Número reproductivo efectivo Rt — AR(2)",
           x=NULL, y="Rt (escala log)") + theme_minimal(base_size=12)

    pc <- ggplot(rmu, aes(fecha)) +
      geom_ribbon(aes(ymin=q025,ymax=q975), fill="#D95F02", alpha=.3) +
      geom_line(aes(y=media), color="#D95F02", linewidth=1) +
      geom_point(data=data.frame(fecha=grilla[1:N_EXOGENO],
                                  obs=y_obs[1:N_EXOGENO]),
                 aes(y=obs), color="black", size=1.5, alpha=.7) +
      scale_x_date(date_breaks="3 days", date_labels="%d %b") +
      labs(title="Componente exógeno µ(t) — prior Gamma(1.401, 0.168)",
           subtitle="Línea: media posterior | Puntos: importados observados",
           x=NULL, y="µ(t) casos/día") + theme_minimal(base_size=12)

    pdf(file.path(RUTA_OUTPUT, "fig3_ajuste_gamma.pdf"), width=12, height=12)
    grid.arrange(pa, pb, pc, ncol=1, heights=c(3,2.5,1.8))
    dev.off()
    cat("  ✓ fig3_ajuste_gamma.pdf\n")
  }, error=function(e) cat("  ✗ fig ajuste:", conditionMessage(e), "\n"))
}

# --- Fig 4: PPC ---
if (!is.null(fit_ref)) {
  tryCatch({
    yrep_m <- extr_mat(fit_ref, "y_rep")
    idx_s  <- sample(seq_len(nrow(yrep_m)), min(200L, nrow(yrep_m)))
    pdf(file.path(RUTA_OUTPUT, "fig4_ppc_gamma.pdf"), width=10, height=5)
    print(ppc_dens_overlay(y=y_obs, yrep=yrep_m[idx_s,],
                           size=.3, alpha=.2) +
      xlim(0, quantile(y_obs,.99)*1.5) +
      labs(title="Posterior Predictive Check — AR(2), NB-2, Exógeno Gamma") +
      theme_minimal())
    dev.off()
    cat("  ✓ fig4_ppc_gamma.pdf\n")
  }, error=function(e) cat("  ✗ PPC:", conditionMessage(e), "\n"))
}

# --- Fig 5: Sensibilidad IS ---
if (!is.null(sens_df) && nrow(sens_df) > 0) {
  tryCatch({
    p5 <- ggplot(sens_df, aes(fecha)) +
      geom_ribbon(aes(ymin=q025,ymax=q975,fill=IS), alpha=.15) +
      geom_line(aes(y=media,color=IS), linewidth=.75) +
      geom_hline(yintercept=1, linetype="dotted") +
      scale_color_manual(values=COLS) +
      scale_fill_manual(values=COLS) +
      scale_x_date(date_breaks="3 weeks", date_labels="%d\n%b") +
      scale_y_log10(breaks=c(.1,.5,1,2,5)) +
      labs(title="Sensibilidad de Rt al IS — AR(2), NB-2, Exógeno Gamma",
           subtitle="5 especificaciones de IS",
           x=NULL, y="Rt", color=NULL, fill=NULL) +
      theme_minimal(base_size=12) + theme(legend.position="bottom")
    ggsave(file.path(RUTA_OUTPUT, "fig5_sensibilidad_IS_gamma.pdf"),
           p5, width=12, height=5)
    cat("  ✓ fig5_sensibilidad_IS_gamma.pdf\n")
  }, error=function(e) cat("  ✗ fig sensibilidad:", conditionMessage(e), "\n"))
}

# --- Fig 6: Posterior parámetros AR(2) ---
if (!is.null(fit_ref)) {
  tryCatch({
    rho1_v  <- as.vector(extr_mat(fit_ref, "rho1"))
    rho2_v  <- as.vector(extr_mat(fit_ref, "rho2"))
    phi_v   <- as.vector(extr_mat(fit_ref, "phi"))
    sig_v   <- as.vector(extr_mat(fit_ref, "sigma_epsilon"))
    mu_rt_v <- as.vector(extr_mat(fit_ref, "mu_rt"))

    pdf(file.path(RUTA_OUTPUT, "fig6_posterior_params_gamma.pdf"),
        width=15, height=5)
    par(mfrow=c(1,5))
    for (v_list in list(
      list(v=rho1_v,  pr=0.7, nm="ρ₁",   xl=c(0,1)),
      list(v=rho2_v,  pr=0.1, nm="ρ₂",   xl=c(0,1)),
      list(v=phi_v,   pr=6.0, nm="φ",    xl=c(0,20)),
      list(v=sig_v,   pr=0.1, nm="σ_ε",  xl=c(0,.5)),
      list(v=mu_rt_v, pr=0.0, nm="μ_Rt", xl=c(-1,1))
    )) {
      hist(v_list$v, breaks=50, col="#8DA0CB", border="white",
           main=sprintf("%s  media=%.3f", v_list$nm, mean(v_list$v)),
           xlab=v_list$nm, xlim=v_list$xl)
      abline(v=mean(v_list$v),          col="red",  lwd=2)
      abline(v=quantile(v_list$v,c(.025,.975)), col="red", lwd=1.5, lty=2)
      abline(v=v_list$pr,               col="blue", lwd=1.5, lty=3)
    }
    legend("topright", c("Post. media","IC95%","Prior media"),
           col=c("red","red","blue"), lwd=c(2,1.5,1.5), lty=c(1,2,3))
    dev.off()
    cat("  ✓ fig6_posterior_params_gamma.pdf\n")
  }, error=function(e) cat("  ✗ fig posteriors:", conditionMessage(e), "\n"))
}

# ==============================================================================
# SECCIÓN 14: RESUMEN FINAL
# ==============================================================================

cat("\n", strrep("=",60), "\n")
cat("ANÁLISIS COMPLETADO — EXÓGENO GAMMA\n")
cat(strrep("=",60), "\n")
cat("Archivos generados en:\n", RUTA_OUTPUT, "\n\n")
cat("FITS (.rds, completos):\n")
for (nm_si in names(IS_a_usar))
  for (nm_mod in names(modelos_stan))
    cat(sprintf("  fit_gamma_%s_%s.rds\n", nm_si, nm_mod))
cat("\nDRAWS ESENCIALES (.rds, ligeros):\n")
for (nm_si in names(IS_a_usar))
  for (nm_mod in names(modelos_stan))
    cat(sprintf("  draws_gamma_%s_%s.rds\n", nm_si, nm_mod))
cat("\nTABLAS (.csv):\n")
cat("  diagnosticos_convergencia_gamma.csv\n")
cat("  loo_comparison_gamma.csv\n")
cat("  posterior_params_gamma.csv\n")
cat("  sensibilidad_IS_gamma.csv\n")
cat("  tabla_resumen_gamma.csv\n")
cat("\nFIGURAS (.pdf):\n")
cat("  fig1_traceplots_gamma.pdf\n")
cat("  fig2_diagnosticos_gamma.pdf\n")
cat("  fig3_ajuste_gamma.pdf\n")
cat("  fig4_ppc_gamma.pdf\n")
cat("  fig5_sensibilidad_IS_gamma.pdf\n")
cat("  fig6_posterior_params_gamma.pdf\n")

