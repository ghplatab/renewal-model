# =============================================================================
#  AJUSTE DE LAS CORRIDAS AR(1) Y PASEO ALEATORIO
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#    AR(1) con adapt_delta = 0.99
#      adapt_delta es un parametro del muestreador, no de la especificacion:
#      la posterior es la misma y solo cambia como se explora. Esta corrida
#      deja ademas los CSV con energy__ y divergent__, necesarios para el
#      E-BFMI y los graficos de pares.
#
#    Paseo aleatorio sobre log R(t)
#      Especificacion alternativa sin nivel de largo plazo, para evaluar si
#      la persistencia estimada justifica un proceso estacionario.
#
#    Comparacion
#      loo_compare entre ambas, por capacidad predictiva fuera de muestra.
#
#  Ejecutar primero la configuracion de rutas, luego cada corrida por
#  separado. Los CSV crudos se escriben en C:/stan_out: las rutas con tildes
#  corrompen la salida de CmdStan en Windows.
# =============================================================================


# =============================================================================
#  PASO 0 - CONFIGURACION
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(cmdstanr)
  library(posterior); library(loo); library(bayesplot); library(ggplot2)
})

# --- CmdStan -----------------------------------------------------------------
# .Renviron solo se lee cuando R ARRANCA. Si la sesion se abrio antes de que
# ese archivo existiera, la variable CMDSTAN no esta y cmdstanr no encuentra
# nada. Esto lo resuelve sin reiniciar la sesion.
CMDSTAN_DIR <- Sys.getenv("CMDSTAN", unset = "C:/cmdstan/cmdstan-2.39.0")
if (!dir.exists(CMDSTAN_DIR))
  stop("No encuentro CmdStan en: ", CMDSTAN_DIR,
       "\nRevisa la ruta o corre install_cmdstan(dir = 'C:/cmdstan').")
set_cmdstan_path(CMDSTAN_DIR)

BASE      <- file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Escritorio",
                       "Tesis Maestría", "AJUSTES_2026", "AJUSTES_PROFE_JCS")
STAN_ROOT <- file.path(BASE, "Documentos George", "MODELO_STAN")
SEP       <- file.path(BASE, "AJUSTES_SEP_2026")

OUT  <- file.path(SEP, "resultados")
FIG  <- file.path(SEP, "figuras")
STAN <- file.path(SEP, "stan")
for (p in c(OUT, FIG, STAN)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# Rutas ASCII y fuera de OneDrive: CmdStan corrompe CSV en rutas con acento
BUILD <- "C:/stan_build"; CSVOUT <- "C:/stan_out"
dir.create(BUILD,  showWarnings = FALSE, recursive = TRUE)
dir.create(CSVOUT, showWarnings = FALSE, recursive = TRUE)

# Copiamos el .stan del AR(1) sin modificarlo, para no compilar dentro de
# tus carpetas. La copia es identica byte a byte: misma especificacion.
AR1_ORIG <- file.path(STAN_ROOT, "MODELO_STAN_DEFINITIVO",
                      "renewal_bogota_nb2_ar1_gamma.stan")
AR1_STAN <- file.path(STAN, "renewal_bogota_nb2_ar1_gamma.stan")
if (!file.exists(AR1_STAN)) file.copy(AR1_ORIG, AR1_STAN)
RW_STAN  <- file.path(STAN, "renewal_bogota_nb2_rw_gamma.stan")

# --- Datos -------------------------------------------------------------------
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
ALPHA_MU  <- 1.401; BETA_MU <- 0.168

discret_si <- function(pfun, ...) {
  g <- numeric(MAX_SI); g[1] <- pfun(1.5, ...) - pfun(0, ...)
  for (s in 2:MAX_SI) g[s] <- pfun(s + 0.5, ...) - pfun(s - 0.5, ...)
  g[g < 0] <- 0; g / sum(g)
}
g_ref <- discret_si(pgamma, shape = 6.5, rate = 0.62)

datos_stan <- list(T = T_PERIODO, n_exogeno = N_EXOGENO, max_si = MAX_SI,
                   y = y_obs, g = g_ref, alpha_mu = ALPHA_MU, beta_mu = BETA_MU)

# Configuracion IDENTICA a la del documento, salvo adapt_delta
N_CHAINS <- 4L; N_WARMUP <- 1500L; N_SAMPLING <- 1500L; MAX_TREEDEPTH <- 12L

cat("=== PASO 0 OK ===\n")
cat(sprintf("  R %s | CmdStan %s\n", getRversion(), cmdstan_version()))
cat(sprintf("  Datos: %d dias, %s casos\n", length(y_obs),
            format(sum(y_obs), big.mark = ",")))
cat(sprintf("  IS referencia: media %.2f dias\n", sum(1:MAX_SI * g_ref)))
cat(sprintf("  Stan AR(1): %s\n", file.exists(AR1_STAN)))
cat(sprintf("  Stan RW   : %s\n", file.exists(RW_STAN)))


# =============================================================================
#  CORRIDA 1 - AR(1) CON adapt_delta = 0.99
#
#  Linea base a superar (lo que reporta el documento, adapt_delta = 0.95):
#      120 divergencias de 6000 (2.0%)   rhat 1.0036   ESS 1718
#
#  Tiempo estimado: 35-50 min. Con adapt_delta mas alto el paso es menor,
#  las trayectorias mas largas, y cada iteracion cuesta mas.
#  Dejalo corriendo y vete a hacer otra cosa.
# =============================================================================

cat("\n=== CORRIDA 1: AR(1) con adapt_delta = 0.99 ===\n")

init_ar1 <- function() list(mu_rt = 0.0, rho1 = 0.8, sigma_epsilon = 0.1,
                            epsilon_raw = rep(0.0, T_PERIODO),
                            mu_exo = rep(ALPHA_MU / BETA_MU, N_EXOGENO),
                            phi = 4.0)

mod_ar1 <- cmdstan_model(AR1_STAN, dir = BUILD)
t0 <- proc.time()
fit_ar1 <- mod_ar1$sample(
  data = datos_stan, iter_warmup = N_WARMUP, iter_sampling = N_SAMPLING,
  chains = N_CHAINS, parallel_chains = N_CHAINS,
  adapt_delta = 0.99, max_treedepth = MAX_TREEDEPTH,
  init = init_ar1, seed = 42, refresh = 250,
  show_messages = FALSE, output_dir = CSVOUT)
t_ar1 <- (proc.time() - t0)[["elapsed"]] / 60
cat(sprintf("\n  Tiempo: %.1f min\n", t_ar1))

pars_ar1 <- c("mu_rt", "rho1", "sigma_epsilon", "phi")
s1 <- fit_ar1$summary(pars_ar1)
print(s1)

dg1 <- fit_ar1$diagnostic_summary(quiet = TRUE)
cat(sprintf("\n  Divergencias   : %d de %d  (%.2f%%)   [documento: 120, 2.00%%]\n",
            sum(dg1$num_divergent), N_CHAINS * N_SAMPLING,
            100 * sum(dg1$num_divergent) / (N_CHAINS * N_SAMPLING)))
cat(sprintf("  Treedepth max  : %d\n", sum(dg1$num_max_treedepth)))
cat(sprintf("  Rhat max       : %.4f   [documento: 1.0036]\n", max(s1$rhat)))
cat(sprintf("  ESS bulk min   : %.0f     [documento: 1718]\n", min(s1$ess_bulk)))
if (!is.null(dg1$ebfmi))
  cat(sprintf("  E-BFMI por cadena: %s  (< 0.3 es problematico)\n",
              paste(sprintf("%.2f", dg1$ebfmi), collapse = "  ")))

# Guardar el ajuste completo (save_object conserva los draws, saveRDS no)
fit_ar1$save_object(file.path(OUT, "corrida1_ar1_ad099.rds"))
cat("  Guardado: corrida1_ar1_ad099.rds\n")

# --- Pairs plot con las divergencias marcadas --------------
np1 <- nuts_params(fit_ar1)
p1  <- mcmc_pairs(fit_ar1$draws(pars_ar1), np = np1,
                  off_diag_args = list(size = 0.6, alpha = 0.4))
ggsave(file.path(FIG, "corrida1_pairs_ar1.png"), p1,
       width = 9, height = 9, dpi = 150)

# --- Comparacion de posteriores: la clave del "riesgo cero" -----------------
cat("\n  Posterior nueva vs la del documento:\n")
ref <- c(mu_rt = 0.203, rho1 = 0.962, sigma_epsilon = 0.099, phi = 3.263)
for (p in pars_ar1)
  cat(sprintf("    %-14s %7.4f   documento %6.3f   dif %+.4f\n",
              p, s1$mean[s1$variable == p], ref[[p]],
              s1$mean[s1$variable == p] - ref[[p]]))
cat("  (diferencias pequenas = solo error de Monte Carlo, el modelo no cambio)\n")


# =============================================================================
#  CORRIDA 3 - RANDOM WALK
#
#  La pregunta: si dejamos que log Rt sea un paseo aleatorio, sin media a la
#  que revertir, se ajusta peor que el AR(1)?
#
#    Si el AR(1) resulta claramente superior -> la reversion a la media queda
#    respaldada por evidencia predictiva.
#    Si empatan -> mu_R no esta identificado por los datos y la conclusion
#    sobre el nivel de largo plazo debe reencuadrarse.
#
#  Nota: este modelo NO tiene mu_rt. En un random walk el nivel global no es
#  identificable por separado del estado inicial. Esa ausencia es el punto.
# =============================================================================

cat("\n=== CORRIDA 3: RANDOM WALK ===\n")

init_rw <- function() list(log_Rt1 = 0.0, sigma_epsilon = 0.1,
                           epsilon_raw = rep(0.0, T_PERIODO - 1),
                           mu_exo = rep(ALPHA_MU / BETA_MU, N_EXOGENO),
                           phi = 4.0)

mod_rw <- cmdstan_model(RW_STAN, dir = BUILD)
t0 <- proc.time()
fit_rw <- mod_rw$sample(
  data = datos_stan, iter_warmup = N_WARMUP, iter_sampling = N_SAMPLING,
  chains = N_CHAINS, parallel_chains = N_CHAINS,
  adapt_delta = 0.99, max_treedepth = MAX_TREEDEPTH,
  init = init_rw, seed = 42, refresh = 250,
  show_messages = FALSE, output_dir = CSVOUT)
t_rw <- (proc.time() - t0)[["elapsed"]] / 60
cat(sprintf("\n  Tiempo: %.1f min\n", t_rw))

pars_rw <- c("log_Rt1", "sigma_epsilon", "phi")
s3 <- fit_rw$summary(pars_rw)
print(s3)

dg3 <- fit_rw$diagnostic_summary(quiet = TRUE)
cat(sprintf("\n  Divergencias  : %d de %d (%.2f%%)\n",
            sum(dg3$num_divergent), N_CHAINS * N_SAMPLING,
            100 * sum(dg3$num_divergent) / (N_CHAINS * N_SAMPLING)))
cat(sprintf("  Treedepth max : %d\n", sum(dg3$num_max_treedepth)))
cat(sprintf("  Rhat max      : %.4f\n", max(s3$rhat)))
cat(sprintf("  ESS bulk min  : %.0f\n", min(s3$ess_bulk)))
if (!is.null(dg3$ebfmi))
  cat(sprintf("  E-BFMI        : %s\n", paste(sprintf("%.2f", dg3$ebfmi), collapse = "  ")))

fit_rw$save_object(file.path(OUT, "corrida3_rw.rds"))
cat("  Guardado: corrida3_rw.rds\n")


# =============================================================================
#  COMPARACION AR(1) vs RANDOM WALK
#
#  Por capacidad PREDICTIVA (elpd), no por convergencia: el comportamiento
#  computacional informa sobre la geometria de la posterior, no sobre la
#  adecuacion del modelo a los datos.
#
#  Como leer loo_compare: la fila de arriba es el mejor modelo. Lo que
#  importa es elpd_diff frente a se_diff. Regla de trabajo habitual:
#    |elpd_diff| > 4 * se_diff  -> diferencia clara
#    |elpd_diff| < 2 * se_diff  -> indistinguibles
# =============================================================================

cat("\n=== COMPARACION AR(1) vs RANDOM WALK ===\n")

ll_ar1 <- fit_ar1$draws("log_lik", format = "matrix")
ll_rw  <- fit_rw$draws("log_lik",  format = "matrix")
cid    <- rep(1:N_CHAINS, each = N_SAMPLING)

loo_ar1 <- loo(ll_ar1, r_eff = relative_eff(exp(ll_ar1), chain_id = cid))
loo_rw  <- loo(ll_rw,  r_eff = relative_eff(exp(ll_rw),  chain_id = cid))

cat("\n--- AR(1) ---\n"); print(loo_ar1)
cat("\n--- Random walk ---\n"); print(loo_rw)

cmp <- loo_compare(list(AR1 = loo_ar1, RandomWalk = loo_rw))
cat("\n--- loo_compare ---\n"); print(cmp)

ed <- abs(cmp[2, "elpd_diff"]); se <- cmp[2, "se_diff"]
cat(sprintf("\n  elpd_diff = %.1f   se_diff = %.1f   razon = %.1f\n", ed, se, ed / se))
cat(sprintf("  Gana: %s\n", rownames(cmp)[1]))

# Nota: escrito como UNA sola expresion a proposito. Un `else` que empieza
# una linea en el nivel superior es error de sintaxis cuando se ejecuta linea
# por linea en la consola, aunque funcione al hacer source() del script entero.
veredicto <- if (ed > 4 * se) {
  "  -> Diferencia CLARA entre ambas especificaciones.\n"
} else if (ed < 2 * se) {
  "  -> INDISTINGUIBLES. Hay que reencuadrar la conclusion.\n"
} else {
  "  -> Zona gris. Reportar ambos y ser prudente.\n"
}
cat(veredicto)

# k de Pareto de los dos, para no comparar sobre aproximaciones malas
cat(sprintf("\n  Pareto k max  AR(1) = %.3f | RW = %.3f\n",
            max(loo_ar1$diagnostics$pareto_k), max(loo_rw$diagnostics$pareto_k)))

saveRDS(list(loo_ar1 = loo_ar1, loo_rw = loo_rw, comparacion = cmp,
             tiempos = c(ar1_min = t_ar1, rw_min = t_rw)),
        file.path(OUT, "corrida_comparacion_ar1_vs_rw.rds"))
write.csv(as.data.frame(cmp), file.path(OUT, "corrida_loo_ar1_vs_rw.csv"))
cat("\nGuardado: corrida_loo_ar1_vs_rw.csv\n")

cat("\n=== CORRIDAS AR(1) Y PASEO ALEATORIO COMPLETAS ===\n")
