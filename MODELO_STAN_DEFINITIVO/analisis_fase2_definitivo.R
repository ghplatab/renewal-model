# ==============================================================================
#  FASE 2 + MODELO DEFINITIVO — COVID-19 BOGOTÁ
#  Sensibilidad al proceso latente (AR1 vs AR2) y al exógeno (Gamma vs Exp)
#
#  Modelos a correr:
#    M1: AR(1) + Gamma exógeno       → sensibilidad proceso latente
#    M4: AR(2) + Exponencial exógeno → sensibilidad prior exógeno
#    MD: modelo definitivo (ganador de M1/M4) con 1500+1500 iteraciones
#
#  IS: referencia (Bi et al. 2020) para todos los modelos
#  Comparación base: referencia_NB2 de análisis preliminar (AR2+Gamma)
#
#  Carpeta de resultados:
#  C:\Users\GeorgeHarrisonPlataB\OneDrive - SQDM\Documentos George\
#  MODELO_STAN\MODELO_STAN_DEFINITIVO
# ==============================================================================

# ==============================================================================
# SECCIÓN 1: CONFIGURACIÓN
# ==============================================================================

# Rtools (solo Windows — ajustar versión si es necesario)
if (.Platform$OS.type == "windows") {
  rtools_path <- "C:/rtools44/usr/bin;C:/rtools44/mingw64/bin"
  Sys.setenv(PATH = paste(rtools_path, Sys.getenv("PATH"), sep = ";"))
}

library(here)
RUTA_DATOS  <- here("MODELO_STAN_4", "datos_agregados.xlsx")
RUTA_STAN   <- here("MODELO_STAN_DEFINITIVO")
RUTA_OUTPUT <- RUTA_STAN

# Ruta al fit de referencia del análisis preliminar (AR2+Gamma, IS referencia)
RUTA_PRELIM <- here("MODELO_STAN_4", "EXOGENO_GAMMA_ANALISIS_PRELIMINAR")
FIT_BASE    <- here("MODELO_STAN_4", "EXOGENO_GAMMA_ANALISIS_PRELIMINAR",
                    "fit_gamma_referencia_NB2.rds")

if (!dir.exists(RUTA_OUTPUT)) dir.create(RUTA_OUTPUT, recursive = TRUE)
setwd(RUTA_OUTPUT)

# --- Hiperparámetros exógeno ---
ALPHA_MU  <- 1.401   # Gamma shape (MoM)
BETA_MU   <- 0.168   # Gamma rate  (MoM)
LAMBDA_MU <- 0.120   # Exponencial rate = 1/8.33

# --- Control de ejecución ---
# Fase 2: configuración igual al análisis preliminar para comparabilidad
N_CHAINS_F2   <- 4
N_WARMUP_F2   <- 1000
N_SAMPLING_F2 <- 1000
# Modelo definitivo: muestreo extendido
N_CHAINS_MD   <- 4
N_WARMUP_MD   <- 1500
N_SAMPLING_MD <- 1500

ADAPT_DELTA   <- 0.95
MAX_TREEDEPTH <- 12

cat("=== FASE 2 + MODELO DEFINITIVO ===\n")
cat(sprintf("Fase 2:    %d cadenas × %d warmup + %d sampling\n",
            N_CHAINS_F2, N_WARMUP_F2, N_SAMPLING_F2))
cat(sprintf("Definitivo:%d cadenas × %d warmup + %d sampling\n",
            N_CHAINS_MD, N_WARMUP_MD, N_SAMPLING_MD))
cat(sprintf("Exógeno Gamma:  alpha=%.3f, beta=%.3f\n", ALPHA_MU, BETA_MU))
cat(sprintf("Exógeno Exp:    lambda=%.3f\n", LAMBDA_MU))

# --- Paquetes ---
suppressPackageStartupMessages({
  library(readxl)
  library(loo)
  library(posterior)
  library(bayesplot)
  library(cmdstanr)
  library(ggplot2)
  library(gridExtra)
  library(scales)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

# ==============================================================================
# SECCIÓN 2: DATOS
# ==============================================================================

cat("\n=== CARGANDO DATOS ===\n")

df_raw        <- as.data.frame(read_excel(RUTA_DATOS, sheet="NO_IMPORTADOS"))
names(df_raw) <- c("fecha", "casos")
df_raw        <- df_raw[!grepl("nan", as.character(df_raw$fecha),
                                ignore.case=TRUE), ]
df_raw        <- df_raw[!is.na(df_raw$fecha), ]
df_raw$fecha  <- as.Date(as.numeric(df_raw$fecha), origin="1899-12-30")
df_raw$casos  <- as.integer(df_raw$casos)
df_local      <- df_raw[order(df_raw$fecha), ]

FECHA_INICIO <- as.Date("2020-03-14")
FECHA_FIN    <- FECHA_INICIO + 179L
T_PERIODO    <- 180L
N_EXOGENO    <- 12L
MAX_SI       <- 30L

grilla <- seq(FECHA_INICIO, FECHA_FIN, by="day")
df_periodo        <- merge(data.frame(fecha=grilla),
                           df_local[df_local$fecha <= FECHA_FIN,],
                           by="fecha", all.x=TRUE)
df_periodo        <- df_periodo[order(df_periodo$fecha),]
df_periodo$casos[is.na(df_periodo$casos)] <- 0L
y_obs             <- as.integer(df_periodo$casos)

cat(sprintf("Período: %s → %s (%d días)\n",
            FECHA_INICIO, FECHA_FIN, T_PERIODO))
cat(sprintf("Casos totales: %s | Media: %.1f/día\n",
            format(sum(y_obs), big.mark=","), mean(y_obs)))

# ==============================================================================
# SECCIÓN 3: INTERVALO SERIAL (referencia: Bi et al. 2020)
# ==============================================================================

discret_si <- function(pfun, ...) {
  g <- numeric(MAX_SI)
  g[1] <- pfun(1.5,...) - pfun(0,...)
  for (s in 2:MAX_SI) g[s] <- pfun(s+0.5,...) - pfun(s-0.5,...)
  g[g < 0] <- 0; g / sum(g)
}

g_ref <- discret_si(pgamma, shape=6.5, rate=0.62)
cat(sprintf("IS referencia: media=%.2f días\n", sum(1:MAX_SI * g_ref)))

# ==============================================================================
# SECCIÓN 4: FUNCIONES AUXILIARES
# ==============================================================================

extr_mat <- function(fit, param)
  fit$draws(param, format="matrix")

extr_ll <- function(fit)
  fit$draws("log_lik", format="matrix")

get_rhat_max <- function(fit, tiene_rho2=TRUE) {
  pars <- c("mu_rt","rho1","sigma_epsilon","phi")
  if (tiene_rho2) pars <- c(pars, "rho2")
  tryCatch(max(fit$summary(pars)$rhat, na.rm=TRUE), error=function(e) NA_real_)
}

get_ess_min <- function(fit, tiene_rho2=TRUE) {
  pars <- c("mu_rt","rho1","sigma_epsilon","phi")
  if (tiene_rho2) pars <- c(pars, "rho2")
  tryCatch(min(fit$summary(pars)$ess_bulk, na.rm=TRUE), error=function(e) NA_real_)
}

get_divergencias <- function(fit) {
  tryCatch(sum(fit$diagnostic_summary(quiet=TRUE)$num_divergent),
           error=function(e) NA_integer_)
}

calc_loo <- function(fit, n_chains, n_sampling) {
  ll    <- extr_ll(fit)
  r_eff <- relative_eff(exp(ll),
                         chain_id=rep(1:n_chains, each=n_sampling))
  loo(ll, r_eff=r_eff)
}

fmt_post <- function(fit, param) {
  v <- as.vector(extr_mat(fit, param))
  sprintf("%.3f [%.3f, %.3f]", mean(v), quantile(v,.025), quantile(v,.975))
}

# ==============================================================================
# SECCIÓN 5: MODELOS FASE 2
# ==============================================================================

cat("\n=== FASE 2: SENSIBILIDAD DE COMPONENTES ===\n")

modelos_f2 <- list(
  M1_AR1_Gamma = list(
    stan    = file.path(RUTA_STAN, "renewal_bogota_nb2_ar1_gamma.stan"),
    data    = list(T=T_PERIODO, n_exogeno=N_EXOGENO, max_si=MAX_SI,
                   y=y_obs, g=g_ref, alpha_mu=ALPHA_MU, beta_mu=BETA_MU),
    init_fn = function() list(mu_rt=0.0, rho1=0.8,
                               sigma_epsilon=0.1,
                               epsilon_raw=rep(0.0, T_PERIODO),
                               mu_exo=rep(ALPHA_MU/BETA_MU, N_EXOGENO),
                               phi=4.0),
    rho2    = FALSE,
    desc    = "AR(1) + Gamma exógeno"
  ),
  M4_AR2_Exp = list(
    stan    = file.path(RUTA_STAN, "renewal_bogota_nb2_ar2_exponencial.stan"),
    data    = list(T=T_PERIODO, n_exogeno=N_EXOGENO, max_si=MAX_SI,
                   y=y_obs, g=g_ref, lambda_mu=LAMBDA_MU),
    init_fn = function() list(mu_rt=0.0, rho1=0.7, rho2=0.1,
                               sigma_epsilon=0.1,
                               epsilon_raw=rep(0.0, T_PERIODO),
                               mu_exo=rep(1/LAMBDA_MU, N_EXOGENO),
                               phi=4.0),
    rho2    = TRUE,
    desc    = "AR(2) + Exponencial exógeno"
  )
)

fits_f2 <- list()
for (nm in names(modelos_f2)) {
  rds_file <- file.path(RUTA_OUTPUT, paste0("fit_fase2_", nm, ".rds"))
  m        <- modelos_f2[[nm]]

  if (file.exists(rds_file)) {
    cat(sprintf("  [%s] Cargando existente\n", nm))
    fits_f2[[nm]] <- readRDS(rds_file)
    next
  }

  cat(sprintf("  [%s] Ajustando: %s ...\n", nm, m$desc))
  t_ini <- proc.time()

  tryCatch({
    mod <- cmdstan_model(m$stan)
    fit <- mod$sample(
      data            = m$data,
      iter_warmup     = N_WARMUP_F2,
      iter_sampling   = N_SAMPLING_F2,
      chains          = N_CHAINS_F2,
      parallel_chains = min(N_CHAINS_F2, parallel::detectCores()),
      adapt_delta     = ADAPT_DELTA,
      max_treedepth   = MAX_TREEDEPTH,
      init            = m$init_fn,
      seed            = 42,
      refresh         = 200,
      show_messages   = FALSE
    )
    fits_f2[[nm]] <- fit
    saveRDS(fit, rds_file)
    cat(sprintf("    ✓ Guardado (%.1f min)\n",
                (proc.time()-t_ini)[["elapsed"]]/60))
  }, error=function(e)
    cat(sprintf("    ✗ ERROR: %s\n", conditionMessage(e))))
}

# ==============================================================================
# SECCIÓN 6: DIAGNÓSTICOS FASE 2
# ==============================================================================

cat("\n=== DIAGNÓSTICOS FASE 2 ===\n")
cat(sprintf("%-25s %8s %10s %12s\n",
            "Modelo", "Rhat_max", "ESS_min", "Divergencias"))
cat(strrep("-", 58), "\n")

diag_f2 <- do.call(rbind, lapply(names(fits_f2), function(nm) {
  fit    <- fits_f2[[nm]]
  rho2   <- modelos_f2[[nm]]$rho2
  rhat   <- get_rhat_max(fit, rho2)
  ess    <- get_ess_min(fit, rho2)
  ndiv   <- get_divergencias(fit)
  estado <- if(!is.na(rhat) && rhat<1.01 && !is.na(ess) && ess>400) "✓" else "⚠"
  cat(sprintf("  %s %-23s %8.4f %10.0f %12d\n",
              estado, nm, rhat, ess, ndiv))
  data.frame(Modelo=nm, Rhat_max=rhat, ESS_min=ess,
             Divergencias=ndiv, stringsAsFactors=FALSE)
}))

# ==============================================================================
# SECCIÓN 7: LOO-CV COMPARACIÓN FASE 2
# ==============================================================================

cat("\n=== LOO-CV: COMPARACIÓN FASE 2 vs BASE ===\n")

# Cargar modelo base (AR2+Gamma del análisis preliminar)
loo_lista <- list()

if (file.exists(FIT_BASE)) {
  fit_base <- readRDS(FIT_BASE)
  loo_lista[["AR2_Gamma_base"]] <- calc_loo(fit_base, 4, 1000)
  cat(sprintf("  %-25s elpd=%.1f (SE=%.1f)\n", "AR2_Gamma_base",
              loo_lista[["AR2_Gamma_base"]]$estimates["elpd_loo","Estimate"],
              loo_lista[["AR2_Gamma_base"]]$estimates["elpd_loo","SE"]))
} else {
  cat("  AVISO: fit base no encontrado en:", FIT_BASE, "\n")
}

for (nm in names(fits_f2)) {
  fit  <- fits_f2[[nm]]
  if (is.null(fit)) next
  n_s  <- if (grepl("AR1", nm)) N_SAMPLING_F2 else N_SAMPLING_F2
  loo_lista[[nm]] <- tryCatch(
    calc_loo(fit, N_CHAINS_F2, n_s),
    error=function(e) NULL)
  if (!is.null(loo_lista[[nm]]))
    cat(sprintf("  %-25s elpd=%.1f (SE=%.1f)\n", nm,
                loo_lista[[nm]]$estimates["elpd_loo","Estimate"],
                loo_lista[[nm]]$estimates["elpd_loo","SE"]))
}

if (length(loo_lista) >= 2) {
  comp_f2 <- loo_compare(loo_lista)
  cat("\nComparación LOO-CV:\n")
  print(comp_f2)
  write.csv(as.data.frame(comp_f2),
            file.path(RUTA_OUTPUT, "loo_fase2_comparacion.csv"))
}

# ==============================================================================
# SECCIÓN 8: POSTERIORS FASE 2
# ==============================================================================

cat("\n=== POSTERIORS FASE 2 ===\n")

for (nm in names(fits_f2)) {
  fit  <- fits_f2[[nm]]
  if (is.null(fit)) next
  rho2 <- modelos_f2[[nm]]$rho2
  cat(sprintf("\n--- %s ---\n", nm))
  cat(sprintf("  phi           = %s\n", fmt_post(fit, "phi")))
  cat(sprintf("  rho1          = %s\n", fmt_post(fit, "rho1")))
  if (rho2) cat(sprintf("  rho2          = %s\n", fmt_post(fit, "rho2")))
  cat(sprintf("  sigma_epsilon = %s\n", fmt_post(fit, "sigma_epsilon")))
  cat(sprintf("  mu_rt         = %s\n", fmt_post(fit, "mu_rt")))
}

# ==============================================================================
# SECCIÓN 9: SELECCIÓN MODELO DEFINITIVO
# ==============================================================================

cat("\n=== SELECCIÓN MODELO DEFINITIVO ===\n")
cat("Criterios:\n")
cat("  1. Convergencia (Rhat < 1.01, ESS > 400, divergencias mínimas)\n")
cat("  2. Capacidad predictiva (LOO-CV)\n")
cat("  3. Parsimonia (AR1 preferido si diferencia LOO < 2 unidades)\n")
cat("\nRevisa los resultados anteriores y edita MD_STAN y MD_DATA\n")
cat("antes de correr la Sección 10.\n\n")

# ==============================================================================
# SECCIÓN 10: MODELO DEFINITIVO (editar según resultados Fase 2)
# ==============================================================================
# INSTRUCCIÓN: Después de revisar los resultados de las secciones 6-8,
# editar las variables MD_STAN y MD_DATA según el modelo ganador:
#
#   Si AR(1) ≈ AR(2) en LOO-CV → usar AR(1) por parsimonia
#     MD_STAN <- file.path(RUTA_STAN, "renewal_bogota_nb2_ar1_gamma.stan")
#     MD_DATA <- modelos_f2[["M1_AR1_Gamma"]]$data
#     MD_INIT <- modelos_f2[["M1_AR1_Gamma"]]$init_fn
#
#   Si AR(2) >> AR(1) → usar AR(2)
#     MD_STAN <- file.path(RUTA_STAN, "renewal_bogota_nb2_ar2_gamma.stan")
#     MD_DATA <- modelos_f2[["M1_AR1_Gamma"]]$data  # mismo data con alpha/beta
#     MD_INIT <- modelos_f2[["M1_AR1_Gamma"]]$init_fn
#
#   Para el exógeno: si Exp ≈ Gamma → usar Exp (sigue Mishra)
#                    si Gamma >> Exp → usar Gamma



# --- EDITAR AQUÍ según resultados ---
MD_STAN <- file.path(RUTA_STAN, 
                     "renewal_bogota_nb2_ar1_gamma.stan")
MD_DATA <- list(T=T_PERIODO, n_exogeno=N_EXOGENO, max_si=MAX_SI,
                y=y_obs, g=g_ref, 
                alpha_mu=ALPHA_MU, beta_mu=BETA_MU)
MD_INIT <- function() list(
  mu_rt         = 0.0,
  rho1          = 0.7,
  sigma_epsilon = 0.1,
  epsilon_raw   = rep(0.0, T_PERIODO),
  mu_exo        = rep(ALPHA_MU/BETA_MU, N_EXOGENO),
  phi           = 4.0
)
MD_NOMBRE <- "AR1_Gamma_definitivo"
# ------------------------------------

rds_md <- file.path(RUTA_OUTPUT, paste0("fit_", MD_NOMBRE, ".rds"))

if (file.exists(rds_md)) {
  cat(sprintf("Cargando modelo definitivo existente: %s\n", MD_NOMBRE))
  fit_md <- readRDS(rds_md)
} else {
  cat(sprintf("Ajustando modelo definitivo: %s\n", MD_NOMBRE))
  cat(sprintf("Configuración: %d cadenas × %d warmup + %d sampling\n",
              N_CHAINS_MD, N_WARMUP_MD, N_SAMPLING_MD))
  t_ini <- proc.time()
  mod_md <- cmdstan_model(MD_STAN)
  fit_md <- mod_md$sample(
    data            = MD_DATA,
    iter_warmup     = N_WARMUP_MD,
    iter_sampling   = N_SAMPLING_MD,
    chains          = N_CHAINS_MD,
    parallel_chains = min(N_CHAINS_MD, parallel::detectCores()),
    adapt_delta     = ADAPT_DELTA,
    max_treedepth   = MAX_TREEDEPTH,
    init            = MD_INIT,
    seed            = 42,
    refresh         = 300,
    show_messages   = FALSE
  )
  saveRDS(fit_md, rds_md)
  cat(sprintf("✓ Modelo definitivo guardado (%.1f min)\n",
              (proc.time()-t_ini)[["elapsed"]]/60))
}

# ==============================================================================
# SECCIÓN 11: DIAGNÓSTICOS MODELO DEFINITIVO
# ==============================================================================

cat("\n=== DIAGNÓSTICOS MODELO DEFINITIVO ===\n")

rhat_md <- get_rhat_max(fit_md, grepl("AR2", MD_NOMBRE))
ess_md  <- get_ess_min(fit_md,  grepl("AR2", MD_NOMBRE))
div_md  <- get_divergencias(fit_md)

cat(sprintf("  Rhat_max:     %.4f %s\n", rhat_md,
            if(!is.na(rhat_md) && rhat_md<1.01) "✓" else "⚠"))
cat(sprintf("  ESS_min:      %.0f %s\n", ess_md,
            if(!is.na(ess_md) && ess_md>400) "✓" else "⚠"))
cat(sprintf("  Divergencias: %d\n", div_md))

# Diagnóstico detallado
diags_md <- fit_md$diagnostic_summary(quiet=TRUE)
cat(sprintf("  Treedepth hits: %d (%.1f%%)\n",
            sum(diags_md$num_max_treedepth),
            100*sum(diags_md$num_max_treedepth) /
              (N_CHAINS_MD * N_SAMPLING_MD)))


# ==============================================================================
# Checkeo de donde se concentran las Divergencias


# Verificar dimensiones
cat("Dims diag_draws:", dim(diag_draws), "\n")
cat("Dims draws_params:", dim(draws_params), "\n")

# Solución: usar índices directamente
div_idx    <- which(diag_draws[, "divergent__"] == 1)
nodiv_idx  <- which(diag_draws[, "divergent__"] == 0)

cat(sprintf("N divergentes: %d | N no divergentes: %d\n",
            length(div_idx), length(nodiv_idx)))

cat("\nMedias en draws divergentes vs no divergentes:\n")
for (p in c("mu_rt","rho1","sigma_epsilon","phi")) {
  v       <- draws_params[[p]]
  m_div   <- mean(v[div_idx])
  m_nodiv <- mean(v[nodiv_idx])
  cat(sprintf("  %-15s div=%.3f  no_div=%.3f  diff=%.3f\n",
              p, m_div, m_nodiv, abs(m_div - m_nodiv)))
}

# ==============================================================================

# ==============================================================================
# SECCIÓN 12: POSTERIORS MODELO DEFINITIVO
# ==============================================================================

cat("\n=== POSTERIORS MODELO DEFINITIVO ===\n")

tiene_rho2 <- grepl("AR2", MD_NOMBRE)
cat(sprintf("  phi           = %s\n", fmt_post(fit_md, "phi")))
cat(sprintf("  rho1          = %s\n", fmt_post(fit_md, "rho1")))
if (tiene_rho2)
  cat(sprintf("  rho2          = %s\n", fmt_post(fit_md, "rho2")))
cat(sprintf("  sigma_epsilon = %s\n", fmt_post(fit_md, "sigma_epsilon")))
cat(sprintf("  mu_rt         = %s\n", fmt_post(fit_md, "mu_rt")))
cat(sprintf("  Rt medio      = %.3f\n",
            mean(exp(as.vector(extr_mat(fit_md, "mu_rt"))))))

# Guardar tabla posterior
pars_md <- c("phi","rho1","sigma_epsilon","mu_rt")
if (tiene_rho2) pars_md <- append(pars_md, "rho2", after=2)
post_md <- do.call(rbind, lapply(pars_md, function(p) {
  v <- as.vector(extr_mat(fit_md, p))
  data.frame(parametro=p, media=mean(v),
             q025=quantile(v,.025), q975=quantile(v,.975))
}))
write.csv(post_md,
          file.path(RUTA_OUTPUT, "posterior_modelo_definitivo.csv"),
          row.names=FALSE)
cat("✓ posterior_modelo_definitivo.csv\n")

# ==============================================================================
# SECCIÓN 13: DRAWS ESENCIALES MODELO DEFINITIVO
# ==============================================================================

cat("\n=== GUARDANDO DRAWS ESENCIALES ===\n")

draws_md_file <- file.path(RUTA_OUTPUT,
                            paste0("draws_", MD_NOMBRE, ".rds"))
if (!file.exists(draws_md_file)) {
  pars_esc <- c("mu_rt","rho1","sigma_epsilon","phi")
  if (tiene_rho2) pars_esc <- append(pars_esc, "rho2", after=2)
  draws_md <- list(
    escalares = setNames(
      lapply(pars_esc, function(p) as.vector(extr_mat(fit_md, p))),
      pars_esc),
    Rt       = extr_mat(fit_md, "Rt"),
    log_lik  = extr_ll(fit_md),
    log_Rt   = extr_mat(fit_md, "log_Rt"),
    y_pred   = extr_mat(fit_md, "y_pred"),
    f        = extr_mat(fit_md, "f"),
    meta     = list(modelo=MD_NOMBRE, n_draws=N_CHAINS_MD*N_SAMPLING_MD,
                    fecha=Sys.time())
  )
  saveRDS(draws_md, draws_md_file)
  cat(sprintf("✓ draws guardados (%.1f MB)\n",
              file.size(draws_md_file)/1e6))
}

# ==============================================================================
# SECCIÓN 14: FIGURAS
# ==============================================================================

cat("\n=== GENERANDO FIGURAS ===\n")

EVENTOS <- data.frame(
  fecha = as.Date(c("2020-03-25","2020-06-01","2020-09-01")),
  label = c("Cuarentena","Reapertura\npiloto","Reapertura\ngradual"))

# --- Fig 1: Traceplots modelo definitivo ---
tryCatch({
  pars_trace <- c("mu_rt","rho1","sigma_epsilon","phi")
  if (tiene_rho2) pars_trace <- append(pars_trace, "rho2", after=2)
  draws_arr <- fit_md$draws(pars_trace, format="array")
  p_trace <- mcmc_trace(draws_arr, pars=pars_trace) +
    labs(title=sprintf("Traceplots — %s", MD_NOMBRE)) +
    theme_minimal(base_size=10)
  ggsave(file.path(RUTA_OUTPUT, "fig1_traceplots_definitivo.pdf"),
         p_trace, width=14, height=10)
  cat("  ✓ fig1_traceplots_definitivo.pdf\n")
}, error=function(e) cat("  ✗ traceplots:", conditionMessage(e), "\n"))

# --- Fig 2: Comparación Rt — AR1 vs AR2 vs base ---
tryCatch({
  modelos_rt <- list()
  if (file.exists(FIT_BASE)) {
    fb <- readRDS(FIT_BASE)
    Rt_b <- extr_mat(fb, "Rt")
    modelos_rt[["AR2+Gamma (base)"]] <- data.frame(
      fecha=grilla, media=colMeans(Rt_b),
      q025=apply(Rt_b,2,quantile,.025),
      q975=apply(Rt_b,2,quantile,.975))
  }
  for (nm in names(fits_f2)) {
    fit_i <- fits_f2[[nm]]
    if (is.null(fit_i)) next
    Rt_i  <- extr_mat(fit_i, "Rt")
    modelos_rt[[nm]] <- data.frame(
      fecha=grilla, media=colMeans(Rt_i),
      q025=apply(Rt_i,2,quantile,.025),
      q975=apply(Rt_i,2,quantile,.975))
  }
  sens_df <- do.call(rbind, lapply(names(modelos_rt), function(nm) {
    cbind(modelos_rt[[nm]], Modelo=nm)
  }))

  COLS <- c("AR2+Gamma (base)"="#1F77B4",
            "M1_AR1_Gamma"="#D62728",
            "M4_AR2_Exp"="#2CA02C")

  p_sens <- ggplot(sens_df, aes(fecha)) +
    geom_ribbon(aes(ymin=q025, ymax=q975, fill=Modelo), alpha=.15) +
    geom_line(aes(y=media, color=Modelo), linewidth=.8) +
    geom_hline(yintercept=1, linetype="dotted", linewidth=.8) +
    geom_vline(data=EVENTOS, aes(xintercept=fecha),
               linetype="dashed", color="gray50", linewidth=.6) +
    scale_color_manual(values=COLS) +
    scale_fill_manual(values=COLS) +
    scale_x_date(date_breaks="3 weeks", date_labels="%d\n%b") +
    scale_y_log10(breaks=c(.5,1,2,5)) +
    labs(title="Sensibilidad de Rt: AR(1) vs AR(2) y Gamma vs Exponencial",
         subtitle="IS referencia (Bi et al. 2020) | NB-2",
         x=NULL, y="Rt (escala log)", color=NULL, fill=NULL) +
    theme_minimal(base_size=12) +
    theme(legend.position="bottom")
  ggsave(file.path(RUTA_OUTPUT, "fig2_sensibilidad_componentes.pdf"),
         p_sens, width=12, height=6)
  cat("  ✓ fig2_sensibilidad_componentes.pdf\n")
}, error=function(e) cat("  ✗ fig sensibilidad:", conditionMessage(e), "\n"))

# --- Fig 3: Ajuste modelo definitivo ---
tryCatch({
  f_m  <- extr_mat(fit_md, "f")
  Rt_m <- extr_mat(fit_md, "Rt")

  rf <- data.frame(fecha=grilla, obs=y_obs,
                   pred=colMeans(f_m),
                   q025=apply(f_m,2,quantile,.025),
                   q975=apply(f_m,2,quantile,.975))
  rRt <- data.frame(fecha=grilla,
                    media=colMeans(Rt_m),
                    q025=apply(Rt_m,2,quantile,.025),
                    q975=apply(Rt_m,2,quantile,.975))

  pa <- ggplot(rf, aes(fecha)) +
    geom_col(aes(y=obs), fill="salmon", alpha=.7, width=1) +
    geom_ribbon(aes(ymin=q025,ymax=q975), fill="#2166AC", alpha=.25) +
    geom_line(aes(y=pred), color="#2166AC", linewidth=.9) +
    geom_vline(data=EVENTOS, aes(xintercept=fecha),
               linetype="dashed", color="gray40") +
    scale_y_continuous(labels=comma) +
    scale_x_date(date_breaks="3 weeks", date_labels="%d\n%b") +
    labs(title=sprintf("Casos observados vs predichos — %s", MD_NOMBRE),
         x=NULL, y="Casos/día") +
    theme_minimal(base_size=12)

  pb <- ggplot(rRt, aes(fecha)) +
    geom_hline(yintercept=1, linetype="dotted", linewidth=.8) +
    geom_ribbon(aes(ymin=q025,ymax=q975), fill="#1B7837", alpha=.25) +
    geom_line(aes(y=media), color="#1B7837", linewidth=.9) +
    geom_vline(data=EVENTOS, aes(xintercept=fecha),
               linetype="dashed", color="gray40") +
    scale_x_date(date_breaks="3 weeks", date_labels="%d\n%b") +
    scale_y_log10(breaks=c(.1,.5,1,2,5,10)) +
    labs(title="Número reproductivo efectivo Rt — modelo definitivo",
         x=NULL, y="Rt (escala log)") +
    theme_minimal(base_size=12)

  pdf(file.path(RUTA_OUTPUT, "fig3_ajuste_definitivo.pdf"),
      width=12, height=9)
  grid.arrange(pa, pb, ncol=1, heights=c(3,2.5))
  dev.off()
  cat("  ✓ fig3_ajuste_definitivo.pdf\n")
}, error=function(e) cat("  ✗ fig ajuste:", conditionMessage(e), "\n"))

# --- Fig 4: PPC modelo definitivo ---
tryCatch({
  yrep_m <- extr_mat(fit_md, "y_rep")
  idx_s  <- sample(seq_len(nrow(yrep_m)), min(200L, nrow(yrep_m)))
  pdf(file.path(RUTA_OUTPUT, "fig4_ppc_definitivo.pdf"), width=10, height=5)
  print(ppc_dens_overlay(y=y_obs, yrep=yrep_m[idx_s,],
                         size=.3, alpha=.2) +
    xlim(0, quantile(y_obs,.99)*1.5) +
    labs(title=sprintf("PPC — %s", MD_NOMBRE)) +
    theme_minimal())
  dev.off()
  cat("  ✓ fig4_ppc_definitivo.pdf\n")
}, error=function(e) cat("  ✗ PPC:", conditionMessage(e), "\n"))



### Figura Ajuste
# Generar fig2_sensibilidad_componentes.pdf
library(ggplot2)
library(dplyr)

RUTA_OUTPUT <- "C:/Users/GeorgeHarrisonPlataB/OneDrive - SQDM/Documentos George/MODELO_STAN/MODELO_STAN_DEFINITIVO"
RUTA_PRELIM <- "C:/Users/GeorgeHarrisonPlataB/OneDrive - SQDM/Documentos George/MODELO_STAN/MODELO_STAN_4/EXOGENO_GAMMA_ANALISIS_PRELIMINAR"

grilla <- seq(as.Date("2020-03-14"), as.Date("2020-03-14") + 179L, by="day")

EVENTOS <- data.frame(
  fecha = as.Date(c("2020-03-25", "2020-04-27", "2020-06-01")),
  label = c("Cuarentena\nobligatoria", "Cuarentena\nlocalidades", "Reapertura\npiloto")
)

# Cargar Rt de los tres modelos
modelos_rt <- list()

# Base: AR2+Gamma
fit_base <- readRDS(file.path(RUTA_PRELIM, "fit_gamma_referencia_NB2.rds"))
Rt_b <- fit_base$draws("Rt", format="matrix")
modelos_rt[["AR(2) + Gamma (base)"]] <- data.frame(
  fecha = grilla,
  media = colMeans(Rt_b),
  q025  = apply(Rt_b, 2, quantile, .025),
  q975  = apply(Rt_b, 2, quantile, .975))

# M1: AR1+Gamma
fit_m1 <- readRDS(file.path(RUTA_OUTPUT, "fit_fase2_M1_AR1_Gamma.rds"))
Rt_m1  <- fit_m1$draws("Rt", format="matrix")
modelos_rt[["AR(1) + Gamma (M1)"]] <- data.frame(
  fecha = grilla,
  media = colMeans(Rt_m1),
  q025  = apply(Rt_m1, 2, quantile, .025),
  q975  = apply(Rt_m1, 2, quantile, .975))

# Definitivo: AR1+Gamma extendido
fit_md <- readRDS(file.path(RUTA_OUTPUT, "fit_AR1_Gamma_definitivo.rds"))
Rt_md  <- fit_md$draws("Rt", format="matrix")
modelos_rt[["AR(1) + Gamma (definitivo)"]] <- data.frame(
  fecha = grilla,
  media = colMeans(Rt_md),
  q025  = apply(Rt_md, 2, quantile, .025),
  q975  = apply(Rt_md, 2, quantile, .975))

sens_df <- do.call(rbind, lapply(names(modelos_rt), function(nm)
  cbind(modelos_rt[[nm]], Modelo=nm)))

COLS <- c(
  "AR(2) + Gamma (base)"        = "#AEC6CF",
  "AR(1) + Gamma (M1)"          = "#FF7F0E",
  "AR(1) + Gamma (definitivo)"  = "#1F77B4")

p <- ggplot(sens_df, aes(fecha)) +
  geom_ribbon(aes(ymin=q025, ymax=q975, fill=Modelo), alpha=.15) +
  geom_line(aes(y=media, color=Modelo), linewidth=.85) +
  geom_hline(yintercept=1, linetype="dotted", linewidth=.8, color="gray30") +
  geom_vline(data=EVENTOS, aes(xintercept=fecha),
             linetype="dashed", color="gray50", linewidth=.6) +
  geom_text(data=EVENTOS, aes(x=fecha, label=label),
            y=0.45, size=2.8, hjust=0, nudge_x=1, color="gray40") +
  scale_color_manual(values=COLS) +
  scale_fill_manual(values=COLS) +
  scale_x_date(date_breaks="3 weeks", date_labels="%d\n%b") +
  scale_y_log10(breaks=c(.5, 1, 2, 5),
                labels=c("0.5","1","2","5")) +
  labs(
    title    = "Sensibilidad de Rt: AR(1) vs AR(2) y modelo definitivo",
    subtitle = "IS referencia (Bi et al. 2020) | NB-2 | Exógeno Gamma",
    x        = NULL,
    y        = expression(R[t]~"(escala log)"),
    color    = NULL, fill = NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(RUTA_OUTPUT, "fig2_sensibilidad_componentes.pdf"),
       p, width=12, height=6)
ggsave(file.path(RUTA_OUTPUT, "fig2_sensibilidad_componentes.png"),
       p, width=12, height=6, dpi=300)
cat("✓ fig2_sensibilidad_componentes generada\n")

# ==============================================================================
# SECCIÓN 15: RESUMEN FINAL
# ==============================================================================

cat("\n", strrep("=",60), "\n")
cat("FASE 2 + MODELO DEFINITIVO — COMPLETADO\n")
cat(strrep("=",60), "\n")
cat("Archivos en:", RUTA_OUTPUT, "\n\n")
cat("Fase 2:\n")
cat("  fit_fase2_M1_AR1_Gamma.rds\n")
cat("  fit_fase2_M4_AR2_Exp.rds\n")
cat("  loo_fase2_comparacion.csv\n\n")
cat("Modelo definitivo:\n")
cat(sprintf("  fit_%s.rds\n", MD_NOMBRE))
cat(sprintf("  draws_%s.rds\n", MD_NOMBRE))
cat("  posterior_modelo_definitivo.csv\n\n")
cat("Figuras:\n")
cat("  fig1_traceplots_definitivo.pdf\n")
cat("  fig2_sensibilidad_componentes.pdf\n")
cat("  fig3_ajuste_definitivo.pdf\n")
cat("  fig4_ppc_definitivo.pdf\n")
cat("\nPROXIMO PASO:\n")
cat("  Revisar LOO-CV (Sección 7) y diagnosticos (Sección 6)\n")
cat("  Editar MD_STAN y MD_DATA en Sección 10 si el modelo\n")
cat("  ganador no es AR2+Gamma\n")
