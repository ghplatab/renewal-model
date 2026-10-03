# =============================================================================
#  DIAGNOSTICOS SOBRE LA POSTERIOR DEL MODELO DEFINITIVO
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  Autocontenido: no requiere haber ejecutado el script anterior.
#  Ejecutar primero la configuracion de rutas; cada bloque es independiente.
#
#    - Diagnostico k de Pareto para la aproximacion PSIS-LOO
#    - Verificacion predictiva de la dispersion: Poisson frente a NB-2
#    - Residuos de Pearson e identificacion de dias atipicos
#    - Media y mediana estacionarias de R(t), calculadas draw a draw
#    - Escala de las innovaciones en las diez especificaciones preliminares
#    - Verificacion predictiva a priori
#
#  Ninguno ajusta modelos: todos operan sobre draws ya guardados.
# =============================================================================


# =============================================================================
#  PASO 0 - CONFIGURACION  (correlo siempre primero)
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(loo); library(ggplot2)
})

BASE <- file.path(Sys.getenv("USERPROFILE"),
                  "OneDrive", "Escritorio", "Tesis Maestría",
                  "AJUSTES_2026", "AJUSTES_PROFE_JCS")
STAN_ROOT <- file.path(BASE, "Documentos George", "MODELO_STAN")

# La copia canonica que identificamos en el PASO 1 del script anterior
CANONICA <- file.path(STAN_ROOT, "MODELO_REPO_GITHUB", "MODELO_STAN_DEFINITIVO",
                      "draws_AR1_Gamma_definitivo.rds")
PRELIM   <- file.path(STAN_ROOT, "MODELO_REPO_GITHUB", "MODELO_STAN_4",
                      "EXOGENO_GAMMA_ANALISIS_PRELIMINAR")
DATOS    <- file.path(STAN_ROOT, "MODELO_STAN_4", "datos_agregados.xlsx")

OUT <- file.path(BASE, "AJUSTES_SEP_2026", "resultados")
FIG <- file.path(BASE, "AJUSTES_SEP_2026", "figuras")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# --- Datos observados --------------------------------------------------------
df_raw <- as.data.frame(read_excel(DATOS, sheet = "NO_IMPORTADOS"))
names(df_raw) <- c("fecha", "casos")
df_raw <- df_raw[!grepl("nan", as.character(df_raw$fecha), ignore.case = TRUE), ]
df_raw <- df_raw[!is.na(df_raw$fecha), ]
df_raw$fecha <- as.Date(as.numeric(df_raw$fecha), origin = "1899-12-30")

FECHA_INICIO <- as.Date("2020-03-14")
FECHA_FIN    <- as.Date("2020-09-09")
fechas <- seq(FECHA_INICIO, FECHA_FIN, by = "day")
df     <- left_join(data.frame(fecha = fechas), df_raw, by = "fecha")
df$casos[is.na(df$casos)] <- 0L
y_obs  <- as.integer(df$casos)

# Configuracion del muestreo que produjo estos draws
N_CHAINS <- 4L; N_SAMPLING <- 1500L        # 4 x 1500 = 6000 draws
T_PERIODO <- 180L; N_EXOGENO <- 12L; MAX_SI <- 30L

cat("=== PASO 0 OK ===\n")
cat(sprintf("  Datos: %d dias | %s casos | media %.1f/dia\n",
            length(y_obs), format(sum(y_obs), big.mark = ","), mean(y_obs)))
cat(sprintf("  Draws: %s\n", basename(CANONICA)))
cat(sprintf("  Salidas: %s\n", OUT))


# =============================================================================
#  PASO 3 - TABLA DE PARETO k
#
#  Que es: LOO-CV pregunta que tan bien predice el modelo el dia t si no
#  hubiera visto el dia t. Reajustar 180 veces es inviable, asi que PSIS-LOO
#  lo aproxima reponderando los draws que ya existen. Pareto k mide si esa
#  aproximacion es de fiar para cada dia.
#
#    k < 0.5   bien
#    0.5-0.7   aceptable
#    k > 0.7   la aproximacion esta fallando en ese dia
#    k > 1     inservible
#
#  Por que importa: con un proceso latente autorregresivo fuerte, quitar un
#  dia casi no cambia la posterior porque los vecinos cargan la informacion.
#  Eso produce k altos. Si muchos dias superan 0.7, las comparaciones de elpd
#  entre los cinco intervalos seriales no son fiables.
# =============================================================================

cat("\n=== PASO 3: PARETO k ===\n")

d  <- readRDS(CANONICA)
ll <- as.matrix(d$log_lik)          # 6000 draws x 180 dias
rm(d); gc(verbose = FALSE)

r_eff <- relative_eff(exp(ll), chain_id = rep(1:N_CHAINS, each = N_SAMPLING))
loo_def <- loo(ll, r_eff = r_eff)

print(loo_def)

cat("\n--- Distribucion de los k ---\n")
print(pareto_k_table(loo_def))

k <- loo_def$diagnostics$pareto_k
cat(sprintf("\n  k maximo : %.3f\n", max(k)))
cat(sprintf("  k mediano: %.3f\n", median(k)))
cat(sprintf("  dias con k > 0.7 : %d de %d\n", sum(k > 0.7), length(k)))
cat(sprintf("  dias con k > 1.0 : %d de %d\n", sum(k > 1.0), length(k)))

# Que dias son los problematicos, y cuantos casos tenian
if (any(k > 0.7)) {
  idx <- which(k > 0.7)
  cat("\n  Dias problematicos (los 15 peores):\n")
  ord <- idx[order(k[idx], decreasing = TRUE)][1:min(15, length(idx))]
  for (i in ord)
    cat(sprintf("    dia %3d  %s  k=%.2f  casos=%d\n",
                i, format(fechas[i]), k[i], y_obs[i]))
}

saveRDS(loo_def, file.path(OUT, "paso3_loo_definitivo.rds"))
write.csv(data.frame(dia = 1:length(k), fecha = fechas,
                     casos = y_obs, pareto_k = k),
          file.path(OUT, "paso3_pareto_k.csv"), row.names = FALSE)
cat(sprintf("\nGuardado: paso3_pareto_k.csv y paso3_loo_definitivo.rds\n"))
rm(ll, r_eff); gc(verbose = FALSE)


# =============================================================================
#  PASO 4 - PPC DE DISPERSION
#
#  Idea: el modelo ajustado sabe generar datos falsos. Generamos muchos y
#  medimos si los datos REALES se parecen a un conjunto falso tipico, usando
#  un estadistico que capture dispersion:
#
#      T(y, f) = suma_t (y_t - f_t)^2 / f_t        (discrepancia de Pearson)
#
#  La clave del procedimiento: simulamos tambien datos POISSON
#  desde las MISMAS medias f ya ajustadas. Si el Poisson no alcanza la
#  dispersion observada, tenemos evidencia estadistica contra el Poisson
#  SIN necesitar un ajuste Poisson convergido.
#
#  Como leerlo: la proporcion es un p-valor bayesiano.
#    cerca de 0.5 -> el modelo reproduce bien el estadistico
#    cerca de 0 o 1 -> el modelo no puede generar lo que se observa
# =============================================================================

cat("\n=== PASO 4: PPC DE DISPERSION ===\n")

d       <- readRDS(CANONICA)
f_draws <- as.matrix(d$y_pred)   # f(t) por draw: la media condicional
y_nb2   <- as.matrix(d$y_rep)    # replicas bajo NB-2
phi_dr  <- as.vector(d$escalares$phi)
rm(d); gc(verbose = FALSE)

S <- nrow(f_draws)
set.seed(42)
# Replicas Poisson desde las mismas medias f
y_pois <- matrix(rpois(length(f_draws), as.vector(f_draws)), nrow = S)

Tp <- function(y, f) sum((y - f)^2 / f)

T_obs  <- vapply(1:S, function(i) Tp(y_obs,       f_draws[i, ]), numeric(1))
T_nb2  <- vapply(1:S, function(i) Tp(y_nb2[i, ],  f_draws[i, ]), numeric(1))
T_pois <- vapply(1:S, function(i) Tp(y_pois[i, ], f_draws[i, ]), numeric(1))

p_nb2  <- mean(T_nb2  >= T_obs)
p_pois <- mean(T_pois >= T_obs)

cat("\n  Discrepancia de Pearson (media sobre draws):\n")
cat(sprintf("    datos observados : %10.0f\n", mean(T_obs)))
cat(sprintf("    replicas NB-2    : %10.0f\n", mean(T_nb2)))
cat(sprintf("    replicas Poisson : %10.0f\n", mean(T_pois)))
cat("\n  p-valor bayesiano (proporcion de replicas que superan lo observado):\n")
cat(sprintf("    NB-2    : %.3f   <- cerca de 0.5 = reproduce bien\n", p_nb2))
cat(sprintf("    Poisson : %.3f   <- cerca de 0   = no alcanza la dispersion\n", p_pois))

# --- Estadisticos complementarios que sugiere el plan ------------------------
cat("\n  Estadisticos complementarios (p-valores bayesianos):\n")
comp <- function(fn, etiqueta) {
  o <- fn(y_obs)
  a <- vapply(1:S, function(i) fn(y_nb2[i, ]),  numeric(1))
  b <- vapply(1:S, function(i) fn(y_pois[i, ]), numeric(1))
  cat(sprintf("    %-12s observado=%9.1f | NB-2 p=%.3f | Poisson p=%.3f\n",
              etiqueta, o, mean(a >= o), mean(b >= o)))
}
comp(max, "max(y)")
comp(sd,  "sd(y)")
comp(function(v) sum(diff(v) > 0), "dias al alza")

res_ppc <- data.frame(modelo = c("NB-2", "Poisson"),
                      p_valor_dispersion = c(p_nb2, p_pois),
                      T_medio = c(mean(T_nb2), mean(T_pois)),
                      T_observado = mean(T_obs))
write.csv(res_ppc, file.path(OUT, "paso4_ppc_dispersion.csv"), row.names = FALSE)
cat("\nGuardado: paso4_ppc_dispersion.csv\n")

rm(y_pois, T_obs, T_nb2, T_pois); gc(verbose = FALSE)


# =============================================================================
#  PASO 5 - RESIDUOS Y DIAS ATIPICOS
#
#  Motivo: en el PASO 2 vimos cobertura EXACTA al 95% pero SOBRECOBERTURA en
#  los niveles centrales (72% donde deberia haber 50%). Hipotesis: unos pocos
#  dias atipicos inflan la sobredispersion, lo que ensancha toda la banda
#  predictiva; el resto de los dias quedan mucho mas cerca de f(t).
#
#  Este paso la pone a prueba de tres formas:
#    (a) residuos de Pearson contra el tiempo -> se ven los atipicos
#    (b) residuos por dia de la semana        -> descarta efecto de reporte
#    (c) recalcula la cobertura EXCLUYENDO los atipicos
#
#  Si al quitar los atipicos la sobrecobertura baja, la hipotesis se sostiene.
# =============================================================================

cat("\n=== PASO 5: RESIDUOS Y DIAS ATIPICOS ===\n")

f_med  <- apply(f_draws, 2, mean)          # media posterior de f(t)
phi_m  <- mean(phi_dr)
var_nb <- f_med + f_med^2 / phi_m          # varianza NB-2
resid  <- (y_obs - f_med) / sqrt(var_nb)   # residuo de Pearson estandarizado

cat(sprintf("  Residuos de Pearson: media=%.3f  sd=%.3f\n",
            mean(resid), sd(resid)))
cat("  (si el modelo calibra, la sd deberia rondar 1)\n")

atip <- which(abs(resid) > 2)
cat(sprintf("\n  Dias con |residuo| > 2 : %d de %d\n", length(atip), length(y_obs)))
if (length(atip)) {
  for (i in atip)
    cat(sprintf("    dia %3d  %s  obs=%6d  f=%8.0f  resid=%6.2f\n",
                i, format(fechas[i]), y_obs[i], f_med[i], resid[i]))
}

# --- (b) Por dia de la semana ------------------------------------------------
dsem <- weekdays(fechas)
orden <- c("lunes","martes","miércoles","jueves","viernes","sábado","domingo")
if (!all(unique(dsem) %in% orden))            # por si el locale esta en ingles
  orden <- c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")
cat("\n  Residuo medio por dia de la semana:\n")
for (dd in orden) {
  s <- dsem == dd
  if (any(s)) cat(sprintf("    %-10s n=%2d  media=%6.3f\n", dd, sum(s), mean(resid[s])))
}
cat("  (si un dia se desvia sistematicamente hay efecto de reporte)\n")

# --- (c) Autocorrelacion -----------------------------------------------------
ac <- acf(resid, plot = FALSE, lag.max = 14)
cat("\n  ACF de los residuos (lags 1 a 7):\n")
cat(sprintf("    %s\n", paste(sprintf("lag%d=%.2f", 1:7, ac$acf[2:8]), collapse = "  ")))
cat("  (valores altos = el modelo deja estructura temporal sin capturar)\n")

# --- (d) Cobertura excluyendo los atipicos -----------------------------------
d2    <- readRDS(CANONICA)
M_rep <- as.matrix(d2$y_rep)
rm(d2); gc(verbose = FALSE)

NIVELES <- c(0.50, 0.80, 0.95)
cob_sub <- function(M, y, keep) {
  sapply(NIVELES, function(p) {
    a  <- (1 - p) / 2
    lo <- apply(M[, keep, drop = FALSE], 2, quantile, probs = a,     names = FALSE)
    hi <- apply(M[, keep, drop = FALSE], 2, quantile, probs = 1 - a, names = FALSE)
    mean(y[keep] >= lo & y[keep] <= hi)
  })
}
todos <- seq_along(y_obs)
sin_at <- setdiff(todos, atip)

cat("\n  Cobertura con y sin los dias atipicos:\n")
c_all <- cob_sub(M_rep, y_obs, todos)
c_sin <- cob_sub(M_rep, y_obs, sin_at)
cat("    Nivel   todos    sin atipicos\n")
for (i in seq_along(NIVELES))
  cat(sprintf("    %3.0f%%   %5.1f%%      %5.1f%%\n",
              NIVELES[i]*100, c_all[i]*100, c_sin[i]*100))

# --- Figura de residuos ------------------------------------------------------
df_res <- data.frame(fecha = fechas, resid = resid, casos = y_obs)
p <- ggplot(df_res, aes(fecha, resid)) +
  geom_hline(yintercept = c(-2, 0, 2), linetype = c(2, 1, 2),
             color = c("grey50", "black", "grey50")) +
  geom_point(size = 1.3, alpha = .8) +
  labs(title = "Residuos de Pearson - modelo definitivo AR(1)+Gamma+NB-2",
       subtitle = "Lineas punteadas: +/- 2 desviaciones",
       x = NULL, y = "Residuo de Pearson") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG, "paso5_residuos.png"), p, width = 9, height = 4.5, dpi = 200)

write.csv(df_res, file.path(OUT, "paso5_residuos.csv"), row.names = FALSE)
cat("\nGuardado: paso5_residuos.csv y figuras/paso5_residuos.png\n")
rm(M_rep, f_draws, y_nb2); gc(verbose = FALSE)


# =============================================================================
#  PASO 6 - MEDIA Y MEDIANA ESTACIONARIA DRAW A DRAW
#
#  Tres precisiones sobre las cantidades estacionarias:
#    1. exp(mu_R) es la MEDIANA estacionaria de R(t), no la media.
#    2. La media estacionaria es  exp(mu_R + sigma^2 / (2(1-rho1^2)))
#    3. Ambas deben calcularse draw a draw, no sobre las medias posteriores.
#
#  Y el plan anade: calcularlo draw a draw, no enchufando las medias
#  posteriores. Motivo: desigualdad de Jensen. Para una funcion convexa como
#  la exponencial, la funcion de las medias NO es la media de la funcion, y
#  ademas el enchufado no da intervalo de credibilidad.
#
#  Ojo con la sensibilidad: el termino 1/(1-rho1^2) explota cuando rho1 -> 1.
#  Este paso cuantifica cuanto.
# =============================================================================

cat("\n=== PASO 6: MEDIA Y MEDIANA ESTACIONARIA ===\n")

d   <- readRDS(CANONICA)
esc <- d$escalares
mu  <- as.vector(esc$mu_rt)
rho <- as.vector(esc$rho1)
sig <- as.vector(esc$sigma_epsilon)
rm(d); gc(verbose = FALSE)

mediana_dr <- exp(mu)                              # mediana estacionaria
media_dr   <- exp(mu + sig^2 / (2 * (1 - rho^2)))  # media estacionaria
vida_media <- log(0.5) / log(rho)                  # vida media de las perturbaciones

resumen <- function(v, etiqueta) {
  q <- quantile(v, c(.025, .5, .975), names = FALSE)
  cat(sprintf("  %-22s media=%7.3f  mediana=%7.3f  IC95%%=[%6.3f, %7.3f]\n",
              etiqueta, mean(v), q[2], q[1], q[3]))
  c(media = mean(v), q2.5 = q[1], q50 = q[2], q97.5 = q[3])
}

cat("\n  Calculado draw a draw (correcto):\n")
r_med <- resumen(mediana_dr, "Mediana estacionaria")
r_mea <- resumen(media_dr,   "Media estacionaria")
r_vid <- resumen(vida_media, "Vida media (dias)")

cat("\n  Enchufando las medias posteriores (lo incorrecto, para comparar):\n")
mu_m <- mean(mu); rho_m <- mean(rho); sig_m <- mean(sig)
cat(sprintf("    exp(mu_R)                          = %.4f\n", exp(mu_m)))
cat(sprintf("    exp(mu_R + sig^2/(2(1-rho^2)))     = %.4f\n",
            exp(mu_m + sig_m^2 / (2 * (1 - rho_m^2)))))
cat(sprintf("    documento dice                     = 1.25\n"))

cat("\n  Sensibilidad a rho1 (por que las obs. 4 y 11 van juntas):\n")
for (rr in c(0.90, 0.962, 0.98, 0.99, 0.997))
  cat(sprintf("    rho1=%.3f -> media estacionaria = %8.3f\n",
              rr, exp(mu_m + sig_m^2 / (2 * (1 - rr^2)))))
cat("  (si rho1 se acerca a 1 la media estacionaria deja de tener sentido)\n")

write.csv(rbind(mediana = r_med, media = r_mea, vida_media = r_vid),
          file.path(OUT, "paso6_estacionarias.csv"))
cat("\nGuardado: paso6_estacionarias.csv\n")


# =============================================================================
#  PASO 7 - sigma_epsilon EN LOS 10 AJUSTES PRELIMINARES
#
#  Hipotesis a probar (seccion 7.5 del plan): la verosimilitud Poisson, al no
#  tener parametro de dispersion, obliga al proceso latente a absorber toda la
#  variabilidad. El muestreador infla sigma_epsilon contra una previa
#  N+(0, 0.2) que se resiste, y eso genera la geometria mala.
#
#  Si sigma_epsilon sale sistematicamente MAYOR bajo Poisson -> hipotesis
#  sostenida, y se puede escribir el mecanismo en el manuscrito.
#  Si sale IGUAL -> la explicacion es falsa y hay que omitirla, reportando
#  solo el patron de convergencia.
#
#  Ojo: estos son modelos AR(2), su lista de escalares incluye rho2.
# =============================================================================

cat("\n=== PASO 7: sigma_epsilon POISSON vs NB-2 ===\n")

archivos <- list.files(PRELIM, pattern = "^draws_gamma_.*\\.rds$", full.names = TRUE)
cat(sprintf("  Encontrados %d archivos\n\n", length(archivos)))

filas <- list()
for (a in archivos) {
  nm  <- sub("^draws_gamma_", "", tools::file_path_sans_ext(basename(a)))
  dd  <- readRDS(a)
  e   <- dd$escalares
  s   <- if (!is.null(e$sigma_epsilon)) as.vector(e$sigma_epsilon) else NA_real_
  div <- if (!is.null(dd$diagnostics)) dd$diagnostics$n_div else NA
  rh  <- if (!is.null(dd$diagnostics)) dd$diagnostics$rhat  else NA
  cat(sprintf("  %-26s sigma=%.4f (sd %.4f)  rhat=%s  div=%s\n",
              nm, mean(s), sd(s), signif(rh, 5), div))
  filas[[nm]] <- data.frame(modelo = nm,
                            verosimilitud = ifelse(grepl("Poisson", nm), "Poisson", "NB-2"),
                            sigma_media = mean(s), sigma_sd = sd(s),
                            rhat = as.numeric(rh)[1], n_div = as.numeric(div)[1])
  rm(dd); gc(verbose = FALSE)
}

tab7 <- do.call(rbind, filas)
cat("\n  Promedio por verosimilitud:\n")
print(aggregate(sigma_media ~ verosimilitud, tab7, mean))

cat("\n  LECTURA: si Poisson tiene sigma claramente mayor, el mecanismo del\n")
cat("           plan se sostiene. Si son parecidos, hay que omitirlo.\n")

write.csv(tab7, file.path(OUT, "paso7_sigma_poisson_vs_nb2.csv"), row.names = FALSE)
cat("\nGuardado: paso7_sigma_poisson_vs_nb2.csv\n")


# =============================================================================
#  PASO 8 - PRIOR PREDICTIVE CHECK
#
#  No usa draws, no usa Stan, no usa tus datos. Solo simula de las PREVIAS
#  tal como estan escritas en el modelo, y mira que datos implican.
#
#  Pregunta: mis previas, antes de ver nada, implican epidemias plausibles
#  para Bogota? Si implican R(t) de 15 o brotes de un millon de casos al dia,
#  son demasiado laxas. Si todo sale pegado a R(t)=1, son demasiado
#  informativas. El caso de N(0.7, 0.15) sobre rho1 merece esa comprobacion.
# =============================================================================

cat("\n=== PASO 8: PRIOR PREDICTIVE CHECK ===\n")

# Intervalo serial de referencia (Bi et al. 2020), igual que en el modelo
discret_si <- function(pfun, ...) {
  g <- numeric(MAX_SI); g[1] <- pfun(1.5, ...) - pfun(0, ...)
  for (s in 2:MAX_SI) g[s] <- pfun(s + 0.5, ...) - pfun(s - 0.5, ...)
  g[g < 0] <- 0; g / sum(g)
}
g_ref <- discret_si(pgamma, shape = 6.5, rate = 0.62)

# Sorteo de una normal truncada por rechazo
rnorm_trunc <- function(mean, sd, lo, hi) {
  repeat { v <- rnorm(1, mean, sd); if (v > lo && v < hi) return(v) }
}

set.seed(2026)
N_SIM <- 1000
sim <- data.frame(Rt_ini = NA_real_, Rt_med = NA_real_, Rt_fin = NA_real_,
                  Rt_max = NA_real_, casos_tot = NA_real_, casos_max = NA_real_)[0, ]
n_explotan <- 0

for (s in 1:N_SIM) {
  # --- Sorteo de las PREVIAS del modelo ---
  mu_rt  <- rnorm(1, 0, 0.3)                          # mu_rt ~ N(0, 0.3)
  rho1   <- rnorm_trunc(0.7, 0.15, 0, 1)              # rho1 ~ N(0.7,0.15) T[0,1]
  sigma  <- abs(rnorm(1, 0, 0.2))                     # sigma ~ N+(0, 0.2)
  phi    <- abs(rnorm(1, 0, 5))                       # phi ~ N+(0, 5)
  mu_exo <- rgamma(N_EXOGENO, shape = 1.401, rate = 0.168)

  # --- Proceso AR(1) con inicio estacionario ---
  eps <- numeric(T_PERIODO)
  eps[1] <- (sigma / sqrt(1 - rho1^2)) * rnorm(1)
  for (t in 2:T_PERIODO) eps[t] <- rho1 * eps[t-1] + sigma * rnorm(1)
  Rt <- exp(mu_rt + eps)

  # --- Ecuacion de renovacion ---
  f <- numeric(T_PERIODO)
  for (t in 1:T_PERIODO) {
    mut  <- if (t <= N_EXOGENO) mu_exo[t] else 0
    lagm <- min(t - 1, MAX_SI)
    conv <- if (lagm >= 1) sum(f[t - (1:lagm)] * g_ref[1:lagm]) else 0
    f[t] <- max(mut + Rt[t] * conv, 1e-6)
  }
  if (!all(is.finite(f)) || max(f) > 1e9) { n_explotan <- n_explotan + 1; next }

  y <- rnbinom(T_PERIODO, mu = f, size = max(phi, 1e-3))
  if (!all(is.finite(y))) { n_explotan <- n_explotan + 1; next }

  sim <- rbind(sim, data.frame(
    Rt_ini = Rt[1], Rt_med = Rt[90], Rt_fin = Rt[180], Rt_max = max(Rt),
    casos_tot = sum(y), casos_max = max(y)))
}

cat(sprintf("  Simulaciones validas: %d de %d (%d explotaron)\n",
            nrow(sim), N_SIM, n_explotan))

q <- function(v) quantile(v, c(.025, .25, .5, .75, .975), names = FALSE)
cat("\n  Que implican tus previas (cuantiles 2.5 / 25 / 50 / 75 / 97.5):\n")
for (cl in names(sim))
  cat(sprintf("    %-10s %s\n", cl,
              paste(formatC(q(sim[[cl]]), format = "f", digits = 1, big.mark = ","),
                    collapse = "  ")))

cat(sprintf("\n  REALIDAD para comparar:\n"))
cat(sprintf("    casos_tot observado = %s\n", format(sum(y_obs), big.mark = ",")))
cat(sprintf("    casos_max observado = %s\n", format(max(y_obs), big.mark = ",")))

cat("\n  LECTURA:\n")
cat("    - Si el observado cae DENTRO del rango previo -> previas razonables.\n")
cat("    - Si las previas permiten epidemias absurdas   -> demasiado laxas.\n")
cat("    - Si el observado queda FUERA del rango previo -> previas en conflicto\n")
cat("      con los datos, que es lo mas delicado de reportar.\n")

write.csv(sim, file.path(OUT, "paso8_prior_predictive.csv"), row.names = FALSE)
cat("\nGuardado: paso8_prior_predictive.csv\n")

cat("\n=== DIAGNOSTICOS POSTERIORES COMPLETOS ===\n")
