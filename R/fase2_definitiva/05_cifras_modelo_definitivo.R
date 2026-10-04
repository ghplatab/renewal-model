# =============================================================================
#  CIFRAS DEL MODELO DEFINITIVO
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  Calcula, sobre los ajustes producidos por el script 04, el conjunto de
#  cantidades derivadas que se reportan en el documento:
#
#    - Cobertura de los intervalos predictivos
#    - LOO-CV y diagnostico k de Pareto
#    - Verificacion predictiva de la dispersion, NB-2 frente a Poisson
#    - Residuos de Pearson
#    - Mediana, media estacionaria y vida media del proceso
#    - Figura de la verificacion predictiva posterior, en dos paneles
#    - Dias fuera de la banda del 95 % y represas de notificacion
#    - Comparacion de trayectorias AR(1) frente a paseo aleatorio
#
#  Entrada: resultados/fase2_chinagamma/f2cg_M0_ar1_gamma.rds
#           resultados/fase2_chinagamma/f2cg_M2_rw_gamma.rds
#  Salida : resultados/fase2_chinagamma/def_*.csv y las figuras asociadas
#
#  No ajusta modelos: opera sobre draws ya guardados.
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(loo)
  library(readxl); library(dplyr); library(ggplot2)
  library(scales); library(gridExtra); library(grid); library(here)
})

# --- Rutas relativas a la raíz del repositorio -------------------------------
OLD   <- here("resultados")
OUT   <- here("resultados", "fase2_chinagamma")
FIG   <- here("figuras",    "fase2_chinagamma")
DATOS <- here("data", "datos_agregados.xlsx")
for (p in c(OUT, FIG)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

df_raw <- as.data.frame(read_excel(DATOS, sheet = "NO_IMPORTADOS"))
names(df_raw) <- c("fecha", "casos")
df_raw <- df_raw[!grepl("nan", as.character(df_raw$fecha), ignore.case = TRUE), ]
df_raw <- df_raw[!is.na(df_raw$fecha), ]
df_raw$fecha <- as.Date(as.numeric(df_raw$fecha), origin = "1899-12-30")
fechas <- seq(as.Date("2020-03-14"), as.Date("2020-09-09"), by = "day")
dfp    <- left_join(data.frame(fecha = fechas), df_raw, by = "fecha")
dfp$casos[is.na(dfp$casos)] <- 0L
y_obs  <- as.integer(dfp$casos); T_P <- length(y_obs)
N_CHAINS <- 4L; N_SAMPLING <- 1500L

cat("=== CARGANDO EL DEFINITIVO NUEVO (AR(1) + Gamma, IS China-Gamma) ===\n")
fit    <- readRDS(file.path(OUT, "f2cg_M0_ar1_gamma.rds"))
M_pred <- fit$draws("y_pred",  format = "matrix")
M_rep  <- fit$draws("y_rep",   format = "matrix")
LL     <- fit$draws("log_lik", format = "matrix")
esc    <- fit$draws(c("mu_rt", "rho1", "sigma_epsilon", "phi"), format = "df")
Rt_ar1 <- fit$draws("Rt", format = "matrix")
rm(fit); gc(verbose = FALSE)
stopifnot(ncol(M_rep) == T_P, ncol(M_pred) == T_P, ncol(LL) == T_P)
cat(sprintf("  %d draws x %d dias | %s casos observados\n",
            nrow(M_rep), T_P, format(sum(y_obs), big.mark = ",")))

nuevo <- list()


# =============================================================================
#  PASO 1 - COBERTURA
# =============================================================================

cat("\n=== PASO 1: COBERTURA ===\n")
NIVELES <- c(0.50, 0.80, 0.95)
bandas <- function(M, p) {
  a <- (1 - p) / 2
  list(lo = apply(M, 2, quantile, a,     names = FALSE),
       hi = apply(M, 2, quantile, 1 - a, names = FALSE))
}
cob <- function(M, y, keep = seq_along(y)) sapply(NIVELES, function(p) {
  b <- bandas(M[, keep, drop = FALSE], p)
  mean(y[keep] >= b$lo & y[keep] <= b$hi)
})
cob_rep  <- cob(M_rep,  y_obs)
cob_pred <- cob(M_pred, y_obs)
b95_rep  <- bandas(M_rep, 0.95); b95_pred <- bandas(M_pred, 0.95)
for (i in seq_along(NIVELES))
  cat(sprintf("  %2.0f%%   con y_pred %5.1f%%   con y_rep %5.1f%%\n",
              NIVELES[i] * 100, cob_pred[i] * 100, cob_rep[i] * 100))
ancho_pred <- mean(b95_pred$hi - b95_pred$lo)
ancho_rep  <- mean(b95_rep$hi  - b95_rep$lo)
fuera_rep  <- sum(y_obs < b95_rep$lo | y_obs > b95_rep$hi)
cat(sprintf("  Ancho banda 95%%: %.1f -> %.1f casos (factor %.2f)\n",
            ancho_pred, ancho_rep, ancho_rep / ancho_pred))
cat(sprintf("  Dias fuera de la banda 95%%: %d de %d\n", fuera_rep, T_P))
nuevo$cob50 <- cob_rep[1] * 100; nuevo$cob80 <- cob_rep[2] * 100
nuevo$cob95 <- cob_rep[3] * 100
nuevo$ancho_pred <- ancho_pred; nuevo$ancho_rep <- ancho_rep
nuevo$dias_fuera <- fuera_rep


# =============================================================================
#  PASO 2 - LOO Y PARETO k
# =============================================================================

cat("\n=== PASO 2: LOO Y PARETO k ===\n")
lo_obj <- loo(LL, r_eff = relative_eff(exp(LL),
                                       chain_id = rep(1:N_CHAINS, each = N_SAMPLING)))
print(lo_obj)
k <- lo_obj$diagnostics$pareto_k
cat(sprintf("  k max %.3f | k mediano %.3f | dias k>0.7: %d\n",
            max(k), median(k), sum(k > 0.7)))
nuevo$elpd    <- lo_obj$estimates["elpd_loo", "Estimate"]
nuevo$elpd_se <- lo_obj$estimates["elpd_loo", "SE"]
nuevo$p_loo   <- lo_obj$estimates["p_loo", "Estimate"]
nuevo$k_max   <- max(k); nuevo$k_med <- median(k); nuevo$k_sobre07 <- sum(k > 0.7)
saveRDS(lo_obj, file.path(OUT, "def_loo.rds"))
write.csv(data.frame(dia = 1:T_P, fecha = fechas, casos = y_obs, pareto_k = k),
          file.path(OUT, "def_pareto_k.csv"), row.names = FALSE)


# =============================================================================
#  PASO 3 - PPC DE DISPERSION
# =============================================================================

cat("\n=== PASO 3: PPC DE DISPERSION ===\n")
S <- nrow(M_pred); set.seed(42)
M_pois <- matrix(rpois(length(M_pred), as.vector(M_pred)), nrow = S)
Tp <- function(y, f) sum((y - f)^2 / f)
T_obs  <- vapply(1:S, function(i) Tp(y_obs,       M_pred[i, ]), numeric(1))
T_nb2  <- vapply(1:S, function(i) Tp(M_rep[i, ],  M_pred[i, ]), numeric(1))
T_pois <- vapply(1:S, function(i) Tp(M_pois[i, ], M_pred[i, ]), numeric(1))
cat(sprintf("  T observado %.0f | NB-2 %.0f | Poisson %.0f\n",
            mean(T_obs), mean(T_nb2), mean(T_pois)))
cat(sprintf("  Razon observado/Poisson: %.0f veces\n", mean(T_obs) / mean(T_pois)))
cat(sprintf("  p-valor NB-2 %.3f | Poisson %.3f\n",
            mean(T_nb2 >= T_obs), mean(T_pois >= T_obs)))
comp <- function(fn) {
  o <- fn(y_obs)
  c(nb2  = mean(vapply(1:S, function(i) fn(M_rep[i, ]),  numeric(1)) >= o),
    pois = mean(vapply(1:S, function(i) fn(M_pois[i, ]), numeric(1)) >= o))
}
p_max <- comp(max); p_sd <- comp(sd)
cat(sprintf("  max(y): NB-2 p=%.3f Poisson p=%.3f\n", p_max["nb2"], p_max["pois"]))
cat(sprintf("  sd(y) : NB-2 p=%.3f Poisson p=%.3f\n", p_sd["nb2"],  p_sd["pois"]))
nuevo$T_obs <- mean(T_obs); nuevo$T_nb2 <- mean(T_nb2); nuevo$T_pois <- mean(T_pois)
nuevo$razon_pois <- mean(T_obs) / mean(T_pois)
nuevo$p_nb2 <- mean(T_nb2 >= T_obs); nuevo$p_pois <- mean(T_pois >= T_obs)
nuevo$p_max_nb2 <- p_max[["nb2"]]; nuevo$p_sd_nb2 <- p_sd[["nb2"]]
rm(M_pois, T_obs, T_nb2, T_pois); gc(verbose = FALSE)


# =============================================================================
#  PASO 4 - RESIDUOS DE PEARSON
# =============================================================================

cat("\n=== PASO 4: RESIDUOS ===\n")
f_med <- colMeans(M_pred); phi_m <- mean(esc$phi)
resid <- (y_obs - f_med) / sqrt(f_med + f_med^2 / phi_m)
tercios <- sapply(list(1:60, 61:120, 121:180), function(ix) sd(resid[ix]))
atip <- which(abs(resid) > 2)
ac   <- acf(resid, plot = FALSE, lag.max = 7)$acf[2:8]
cat(sprintf("  sd %.3f | media %.3f\n", sd(resid), mean(resid)))
cat(sprintf("  sd por tercios: %.3f / %.3f / %.3f\n",
            tercios[1], tercios[2], tercios[3]))
cat(sprintf("  ACF max |lag 1-7|: %.3f (banda +/- %.3f)\n",
            max(abs(ac)), 1.96 / sqrt(T_P)))
cat(sprintf("  Dias con |resid| > 2: %d\n", length(atip)))
for (i in atip)
  cat(sprintf("    %s  obs=%5d  f=%6.0f  resid=%5.2f\n",
              format(fechas[i]), y_obs[i], f_med[i], resid[i]))
cat(sprintf("  Dias sin IVA: 19-jun resid=%.3f | 03-jul resid=%.3f\n",
            resid[fechas == as.Date("2020-06-19")],
            resid[fechas == as.Date("2020-07-03")]))
cob_sin <- cob(M_rep, y_obs, setdiff(seq_along(y_obs), atip))
cat(sprintf("  Cobertura 50%% sin atipicos: %.1f%%\n", cob_sin[1] * 100))
nuevo$resid_sd <- sd(resid)
nuevo$tercio1 <- tercios[1]; nuevo$tercio2 <- tercios[2]; nuevo$tercio3 <- tercios[3]
nuevo$n_atip <- length(atip); nuevo$cob50_sin <- cob_sin[1] * 100
write.csv(data.frame(fecha = fechas, obs = y_obs, f = f_med, resid = resid),
          file.path(OUT, "def_residuos.csv"), row.names = FALSE)


# =============================================================================
#  PASO 5 - CANTIDADES ESTACIONARIAS
# =============================================================================

cat("\n=== PASO 5: CANTIDADES ESTACIONARIAS ===\n")
mediana_dr <- exp(esc$mu_rt)
media_dr   <- exp(esc$mu_rt + esc$sigma_epsilon^2 / (2 * (1 - esc$rho1^2)))
vida_dr    <- log(0.5) / log(esc$rho1)
resumen <- function(v, nm) {
  q <- quantile(v, c(.025, .5, .975), names = FALSE)
  cat(sprintf("  %-22s media %.3f  mediana %.3f  IC95 [%.3f, %.3f]\n",
              nm, mean(v), q[2], q[1], q[3]))
  q
}
q_med <- resumen(mediana_dr, "Mediana estacionaria")
q_mea <- resumen(media_dr,   "Media estacionaria")
q_vid <- resumen(vida_dr,    "Vida media (dias)")
cat(sprintf("  P(rho1 > 0.9) = %.3f | P(media estacionaria > 1) = %.3f\n",
            mean(esc$rho1 > 0.9), mean(media_dr > 1)))
nuevo$med_est <- mean(mediana_dr); nuevo$med_lo <- q_med[1]; nuevo$med_hi <- q_med[3]
nuevo$mea_est <- mean(media_dr);   nuevo$mea_lo <- q_mea[1]; nuevo$mea_hi <- q_mea[3]
nuevo$vid_med <- q_vid[2];         nuevo$vid_lo <- q_vid[1]; nuevo$vid_hi <- q_vid[3]


# =============================================================================
#  PASO 6 - FIGURA DEL PPC
# =============================================================================

cat("\n=== PASO 6: FIGURA DEL PPC ===\n")
b80 <- bandas(M_rep, 0.80)
df_ppc <- data.frame(fecha = fechas, y_obs = y_obs,
                     mediana = apply(M_rep, 2, median),
                     lo95 = b95_rep$lo, hi95 = b95_rep$hi,
                     lo80 = b80$lo, hi80 = b80$hi)
p_a <- ggplot(df_ppc, aes(x = fecha)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#AEC6CF", alpha = 0.55) +
  geom_ribbon(aes(ymin = lo80, ymax = hi80), fill = "#1F4E79", alpha = 0.30) +
  geom_line(aes(y = mediana), color = "#1F4E79", linewidth = 0.8) +
  geom_point(aes(y = y_obs), size = 0.8, alpha = 0.85, color = "black") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  scale_y_continuous(labels = label_comma(big.mark = ".")) +
  labs(title = "(a) Verificación predictiva posterior temporal",
       subtitle = "Banda clara: IC95% | Banda oscura: IC80% | Línea: mediana | Puntos: observados",
       x = NULL, y = "Casos diarios") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8.5, color = "gray40"))
df_cal <- data.frame(
  nivel = factor(paste0("IC", NIVELES * 100, "%"),
                 levels = paste0("IC", NIVELES * 100, "%")),
  esperado = NIVELES, observado = cob_rep)
p_b <- ggplot(df_cal, aes(x = nivel)) +
  geom_col(aes(y = observado), fill = "#AEC6CF", width = 0.5) +
  geom_point(aes(y = esperado), color = "#D62728", size = 4.5, shape = 18) +
  geom_text(aes(y = observado, label = sprintf("%.1f%%", observado * 100)),
            vjust = -0.6, size = 3.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  labs(title = "(b) Calibración de los intervalos predictivos",
       subtitle = "Barras: cobertura empírica | Rombos: cobertura nominal",
       x = "Intervalo predictivo", y = "Proporción de días cubiertos") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8.5, color = "gray40"))
fig <- arrangeGrob(p_a, p_b, ncol = 2, widths = c(1.65, 1))
ggsave(file.path(FIG, "fig_ppc_chinagamma.png"), fig, width = 13, height = 5.5, dpi = 300)
ok_pdf <- tryCatch({
  ggsave(file.path(FIG, "fig_ppc_chinagamma.pdf"), fig,
         width = 13, height = 5.5, device = cairo_pdf); TRUE
}, error = function(e) FALSE)
cat(sprintf("  Guardada: fig_ppc_chinagamma.png%s\n",
            if (ok_pdf) " y .pdf" else " (sin PDF)"))


# =============================================================================
#  PASO 7 - DIAS FUERA DE LA BANDA 95% Y REPRESAS DE REPORTE
# =============================================================================

cat("\n=== PASO 7: DIAS FUERA DE LA BANDA 95% ===\n")
pit <- sapply(1:T_P, function(i) mean(M_rep[, i] < y_obs[i]))
ix  <- which(y_obs < b95_rep$lo | y_obs > b95_rep$hi)
tab_fuera <- data.frame(
  fecha = fechas[ix], dia_sem = weekdays(fechas[ix]), obs = y_obs[ix],
  f = round(f_med[ix]), banda_lo = round(b95_rep$lo[ix]),
  banda_hi = round(b95_rep$hi[ix]),
  lado = ifelse(y_obs[ix] > b95_rep$hi[ix], "ARRIBA", "ABAJO"),
  resid = round(resid[ix], 2), pct_predictivo = round(100 * pit[ix], 2))
print(tab_fuera, row.names = FALSE)
cat(sprintf("\n  Arriba: %d | Abajo: %d | Total: %d\n",
            sum(tab_fuera$lado == "ARRIBA"), sum(tab_fuera$lado == "ABAJO"),
            nrow(tab_fuera)))
cat("\n  Dia siguiente a cada dia por debajo (represa de reporte):\n")
for (i in ix[y_obs[ix] < b95_rep$lo[ix]])
  if (i < T_P)
    cat(sprintf("    %s obs=%5d  ->  %s obs=%5d  (f=%5.0f)\n",
                format(fechas[i]), y_obs[i], format(fechas[i + 1]),
                y_obs[i + 1], f_med[i + 1]))
write.csv(tab_fuera, file.path(OUT, "def_dias_fuera_banda95.csv"), row.names = FALSE)


# =============================================================================
#  PASO 8 - TRAYECTORIAS AR(1) vs RANDOM WALK
# =============================================================================

cat("\n=== PASO 8: TRAYECTORIAS AR(1) vs RANDOM WALK ===\n")
fit_rw <- readRDS(file.path(OUT, "f2cg_M2_rw_gamma.rds"))
Rt_rw  <- fit_rw$draws("Rt", format = "matrix"); rm(fit_rw); gc(verbose = FALSE)
resumir <- function(M, etiqueta) data.frame(
  fecha = fechas, modelo = etiqueta,
  mediana = apply(M, 2, median),
  lo = apply(M, 2, quantile, .05, names = FALSE),
  hi = apply(M, 2, quantile, .95, names = FALSE),
  p_menor1 = apply(M, 2, function(v) mean(v < 1)))
s_ar1 <- resumir(Rt_ar1, "AR(1)"); s_rw <- resumir(Rt_rw, "Random walk")
dif <- s_rw$mediana - s_ar1$mediana
cat(sprintf("  Correlacion entre medianas : %.4f\n", cor(s_ar1$mediana, s_rw$mediana)))
cat(sprintf("  Diferencia absoluta media  : %.4f | maxima %.4f (%s)\n",
            mean(abs(dif)), max(abs(dif)), format(fechas[which.max(abs(dif))])))
for (u in c(0.05, 0.10, 0.20))
  cat(sprintf("    dias con |dif| > %.2f : %3d de %d (%.1f%%)\n",
              u, sum(abs(dif) > u), T_P, 100 * mean(abs(dif) > u)))
w_ar1 <- s_ar1$hi - s_ar1$lo; w_rw <- s_rw$hi - s_rw$lo
cat(sprintf("  Ancho medio IC90: AR(1) %.3f | RW %.3f | razon %.2f\n",
            mean(w_ar1), mean(w_rw), mean(w_rw) / mean(w_ar1)))
for (ix2 in list(1:60, 61:120, 121:180))
  cat(sprintf("    dias %3d-%3d : AR(1) %.3f | RW %.3f | razon %.2f\n",
              min(ix2), max(ix2), mean(w_ar1[ix2]), mean(w_rw[ix2]),
              mean(w_rw[ix2]) / mean(w_ar1[ix2])))
cat(sprintf("  Mediana AR(1) dentro de la banda del RW: %.1f%% de los dias\n",
            100 * mean(s_ar1$mediana >= s_rw$lo & s_ar1$mediana <= s_rw$hi)))
write.csv(bind_rows(s_ar1, s_rw),
          file.path(OUT, "def_trayectorias_ar1_vs_rw.csv"), row.names = FALSE)
p <- ggplot(bind_rows(s_ar1, s_rw), aes(fecha, mediana, colour = modelo, fill = modelo)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = c("AR(1)" = "black", "Random walk" = "#E41A1C")) +
  scale_fill_manual(values = c("AR(1)" = "black", "Random walk" = "#E41A1C")) +
  labs(x = NULL, y = expression(R[t]), colour = NULL, fill = NULL,
       title = "Modelo definitivo (IS China-Gamma): AR(1) frente a random walk",
       subtitle = "Mediana posterior y banda del 90%") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "def_trayectorias_ar1_vs_rw.png"), p,
       width = 10, height = 5.5, dpi = 150)


# =============================================================================
#  PASO 9 - SENSIBILIDAD AL INTERVALO SERIAL: 10,5 DIAS vs CHINA-GAMMA
#
#  La columna "nuevo" de consistencia_viejo_vs_nuevo.csv son las cifras del
#  definitivo con el nucleo de 10.5 dias: las que hoy estan redactadas.
# =============================================================================

cat("\n=== PASO 9: SENSIBILIDAD AL INTERVALO SERIAL ===\n")
prev <- read.csv(file.path(OLD, "consistencia_viejo_vs_nuevo.csv"))
etiq <- c(
  cob50 = "Cobertura IC50 (%)", cob80 = "Cobertura IC80 (%)", cob95 = "Cobertura IC95 (%)",
  ancho_pred = "Ancho banda 95 con y_pred", ancho_rep = "Ancho banda 95 con y_rep",
  dias_fuera = "Dias fuera banda 95",
  elpd = "elpd_loo", elpd_se = "SE elpd_loo", p_loo = "p_loo",
  k_max = "Pareto k maximo", k_med = "Pareto k mediano", k_sobre07 = "Dias con k > 0.7",
  T_obs = "Pearson T observado", T_nb2 = "Pearson T NB-2", T_pois = "Pearson T Poisson",
  razon_pois = "Razon observado/Poisson", p_nb2 = "p-valor dispersion NB-2",
  p_pois = "p-valor dispersion Poisson", p_max_nb2 = "p-valor max(y) NB-2",
  p_sd_nb2 = "p-valor sd(y) NB-2",
  resid_sd = "sd residuos Pearson", tercio1 = "sd residuos dias 1-60",
  tercio2 = "sd residuos dias 61-120", tercio3 = "sd residuos dias 121-180",
  n_atip = "Dias con |resid| > 2", cob50_sin = "Cobertura IC50 sin atipicos (%)",
  med_est = "Mediana estacionaria", mea_est = "Media estacionaria",
  vid_med = "Vida media mediana (dias)")
nom <- names(etiq)
viejo_val <- sapply(nom, function(nm) {
  i <- match(etiq[[nm]], prev$cantidad)
  if (is.na(i)) NA_real_ else prev$nuevo[i]
})
tabla <- data.frame(cantidad = unname(etiq[nom]),
                    nucleo_10_5d = unname(viejo_val),
                    china_gamma = unlist(nuevo[nom]), row.names = NULL)
tabla$diferencia <- tabla$china_gamma - tabla$nucleo_10_5d
tabla$revisar <- ifelse(!is.na(tabla$diferencia) &
  abs(tabla$diferencia) >= pmax(0.05 * abs(tabla$nucleo_10_5d), 0.01), "<-- REVISAR", "")
cat(sprintf("\n  %-34s %12s %12s\n", "Cantidad", "IS 10.5 d", "China-Gamma"))
cat("  ", strrep("-", 72), "\n", sep = "")
for (i in seq_len(nrow(tabla)))
  cat(sprintf("  %-34s %12.3f %12.3f  %s\n", tabla$cantidad[i],
              tabla$nucleo_10_5d[i], tabla$china_gamma[i], tabla$revisar[i]))
cat("\n  Intervalos nuevos (IC95):\n")
cat(sprintf("    Mediana estacionaria [%.3f, %.3f] | Media [%.3f, %.3f] | Vida media [%.1f, %.1f]\n",
            nuevo$med_lo, nuevo$med_hi, nuevo$mea_lo, nuevo$mea_hi,
            nuevo$vid_lo, nuevo$vid_hi))
write.csv(tabla, file.path(OUT, "def_viejo_vs_nuevo_cifras.csv"), row.names = FALSE)
cat("\n  Guardado: def_viejo_vs_nuevo_cifras.csv\n")
cat("\n=== CIFRAS DEL MODELO DEFINITIVO COMPLETAS ===\n")
