# =============================================================================
#  CONSISTENCIA DE LAS CIFRAS SOBRE LA CORRIDA FINAL
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  La metodologia declara adapt_delta = 0.99, de modo que todas las cifras
#  reportadas deben provenir de esa corrida. Este script recalcula sobre ella
#  el conjunto completo de cantidades derivadas:
#
#    - Cobertura de los intervalos predictivos
#    - LOO-CV y diagnostico k de Pareto
#    - Verificacion predictiva de la dispersion (NB-2 frente a Poisson)
#    - Residuos de Pearson
#    - Mediana, media estacionaria y vida media, calculadas draw a draw
#    - Figura de la verificacion predictiva posterior, construida con y_rep
#
#  No ajusta modelos: opera sobre draws ya guardados.
#  Puede ejecutarse completo o bloque por bloque.
# =============================================================================


# =============================================================================
#  CONFIGURACION Y CARGA DE DATOS
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(loo)
  library(readxl); library(dplyr); library(ggplot2)
  library(scales); library(gridExtra); library(grid); library(here)
})

# --- Rutas relativas a la raíz del repositorio -------------------------------
OUT   <- here("resultados")
FIG   <- here("figuras")
DATOS <- here("data", "datos_agregados.xlsx")
for (p in c(OUT, FIG)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# --- Datos observados --------------------------------------------------------
df_raw <- as.data.frame(read_excel(DATOS, sheet = "NO_IMPORTADOS"))
names(df_raw) <- c("fecha", "casos")
df_raw <- df_raw[!grepl("nan", as.character(df_raw$fecha), ignore.case = TRUE), ]
df_raw <- df_raw[!is.na(df_raw$fecha), ]
df_raw$fecha <- as.Date(as.numeric(df_raw$fecha), origin = "1899-12-30")

fechas <- seq(as.Date("2020-03-14"), as.Date("2020-09-09"), by = "day")
dfp    <- left_join(data.frame(fecha = fechas), df_raw, by = "fecha")
dfp$casos[is.na(dfp$casos)] <- 0L
y_obs  <- as.integer(dfp$casos)
T_P    <- length(y_obs)

N_CHAINS <- 4L; N_SAMPLING <- 1500L

# --- Corrida final -----------------------------------------------------------
cat("=== CARGANDO CORRIDA FINAL (adapt_delta = 0.99) ===\n")
fit <- readRDS(file.path(OUT, "corrida1_ar1_ad099.rds"))

M_pred <- fit$draws("y_pred",  format = "matrix")   # f(t): la media
M_rep  <- fit$draws("y_rep",   format = "matrix")   # predictivo posterior
LL     <- fit$draws("log_lik", format = "matrix")
esc    <- fit$draws(c("mu_rt", "rho1", "sigma_epsilon", "phi"), format = "df")
rm(fit); gc(verbose = FALSE)

stopifnot(ncol(M_rep) == T_P, ncol(M_pred) == T_P, ncol(LL) == T_P)
cat(sprintf("  %d draws x %d dias | %s casos observados\n",
            nrow(M_rep), T_P, format(sum(y_obs), big.mark = ",")))

nuevo <- list()   # aqui se acumulan las cifras nuevas para la tabla final


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
              NIVELES[i]*100, cob_pred[i]*100, cob_rep[i]*100))
ancho_pred <- mean(b95_pred$hi - b95_pred$lo)
ancho_rep  <- mean(b95_rep$hi  - b95_rep$lo)
fuera_rep  <- sum(y_obs < b95_rep$lo | y_obs > b95_rep$hi)
cat(sprintf("  Ancho banda 95%%: %.1f -> %.1f casos (factor %.2f)\n",
            ancho_pred, ancho_rep, ancho_rep / ancho_pred))
cat(sprintf("  Dias fuera de la banda 95%%: %d de %d\n", fuera_rep, T_P))

nuevo$cob50 <- cob_rep[1]*100; nuevo$cob80 <- cob_rep[2]*100
nuevo$cob95 <- cob_rep[3]*100
nuevo$ancho_pred <- ancho_pred; nuevo$ancho_rep <- ancho_rep
nuevo$dias_fuera <- fuera_rep


# =============================================================================
#  PASO 2 - LOO Y PARETO k
# =============================================================================

cat("\n=== PASO 2: LOO Y PARETO k ===\n")

r_eff <- relative_eff(exp(LL), chain_id = rep(1:N_CHAINS, each = N_SAMPLING))
lo    <- loo(LL, r_eff = r_eff)
print(lo)
k <- lo$diagnostics$pareto_k
cat(sprintf("  k max %.3f | k mediano %.3f | dias k>0.7: %d\n",
            max(k), median(k), sum(k > 0.7)))

nuevo$elpd <- lo$estimates["elpd_loo", "Estimate"]
nuevo$elpd_se <- lo$estimates["elpd_loo", "SE"]
nuevo$p_loo <- lo$estimates["p_loo", "Estimate"]
nuevo$k_max <- max(k); nuevo$k_med <- median(k); nuevo$k_sobre07 <- sum(k > 0.7)
saveRDS(lo, file.path(OUT, "final_loo.rds"))


# =============================================================================
#  PASO 3 - PPC DE DISPERSION
# =============================================================================

cat("\n=== PASO 3: PPC DE DISPERSION ===\n")

S <- nrow(M_pred)
set.seed(42)
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

f_med  <- colMeans(M_pred)
phi_m  <- mean(esc$phi)
resid  <- (y_obs - f_med) / sqrt(f_med + f_med^2 / phi_m)

tercios <- sapply(list(1:60, 61:120, 121:180), function(ix) sd(resid[ix]))
atip    <- which(abs(resid) > 2)
ac      <- acf(resid, plot = FALSE, lag.max = 7)$acf[2:8]

cat(sprintf("  sd %.3f | media %.3f\n", sd(resid), mean(resid)))
cat(sprintf("  sd por tercios: %.3f / %.3f / %.3f\n", tercios[1], tercios[2], tercios[3]))
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
cat(sprintf("  Cobertura 50%% sin atipicos: %.1f%%\n", cob_sin[1]*100))

nuevo$resid_sd <- sd(resid)
nuevo$tercio1 <- tercios[1]; nuevo$tercio2 <- tercios[2]; nuevo$tercio3 <- tercios[3]
nuevo$n_atip <- length(atip); nuevo$cob50_sin <- cob_sin[1]*100


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

nuevo$med_est <- mean(mediana_dr); nuevo$med_lo <- q_med[1]; nuevo$med_hi <- q_med[3]
nuevo$mea_est <- mean(media_dr);   nuevo$mea_lo <- q_mea[1]; nuevo$mea_hi <- q_mea[3]
nuevo$vid_med <- q_vid[2];         nuevo$vid_lo <- q_vid[1]; nuevo$vid_hi <- q_vid[3]


# =============================================================================
#  PASO 6 - FIGURA DEL PPC CORREGIDA
# =============================================================================

cat("\n=== PASO 6: FIGURA DEL PPC ===\n")

b80 <- bandas(M_rep, 0.80)
df_ppc <- data.frame(
  fecha = fechas, y_obs = y_obs,
  mediana = apply(M_rep, 2, median),
  lo95 = b95_rep$lo, hi95 = b95_rep$hi,
  lo80 = b80$lo,     hi80 = b80$hi)

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

ggsave(file.path(FIG, "fig_ppc_corregido.png"), fig,
       width = 13, height = 5.5, dpi = 300)
ok_pdf <- tryCatch({
  ggsave(file.path(FIG, "fig_ppc_corregido.pdf"), fig,
         width = 13, height = 5.5, device = cairo_pdf); TRUE
}, error = function(e) FALSE)
cat(sprintf("  Guardada: fig_ppc_corregido.png%s\n",
            if (ok_pdf) " y .pdf" else " (PDF no disponible, usa el PNG)"))


# =============================================================================
#  COMPARACION ENTRE LAS DOS CONFIGURACIONES DEL MUESTREADOR
#
#  Contrasta cada cantidad derivada entre la corrida con adapt_delta = 0.95 y
#  la corrida con adapt_delta = 0.99. Ambas estiman la misma posterior, de
#  modo que las diferencias deben ser atribuibles a error de Monte Carlo; una
#  discrepancia apreciable indicaria exploracion insuficiente en la primera.
# =============================================================================

cat("\n=== COMPARACION ENTRE CONFIGURACIONES DEL MUESTREADOR ===\n")

viejo <- list(
  cob50 = 72.2, cob80 = 88.3, cob95 = 95.0,
  ancho_pred = 852.5, ancho_rep = 2946.3, dias_fuera = 9,
  elpd = -1233.5, elpd_se = 29.1, p_loo = 18.1,
  k_max = 0.602, k_med = 0.218, k_sobre07 = 0,
  T_obs = 65227, T_nb2 = 73367, T_pois = 180, razon_pois = 362,
  p_nb2 = 0.626, p_pois = 0.000, p_max_nb2 = 0.628, p_sd_nb2 = 0.497,
  resid_sd = 0.761, tercio1 = 0.754, tercio2 = 0.618, tercio3 = 0.896,
  n_atip = 4, cob50_sin = 73.9,
  med_est = 1.249, med_lo = 0.773, med_hi = 1.737,
  mea_est = 1.440, mea_lo = 0.933, mea_hi = 2.205,
  vid_med = 24.6, vid_lo = 5.0, vid_hi = 259)

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
  med_est = "Mediana estacionaria", med_lo = "  IC95 inferior", med_hi = "  IC95 superior",
  mea_est = "Media estacionaria", mea_lo = "  IC95 inferior", mea_hi = "  IC95 superior",
  vid_med = "Vida media mediana (dias)", vid_lo = "  IC95 inferior", vid_hi = "  IC95 superior")

tabla <- data.frame(
  cantidad = unname(etiq[names(viejo)]),
  viejo = unlist(viejo),
  nuevo = unlist(nuevo[names(viejo)]),
  row.names = NULL)
tabla$diferencia <- tabla$nuevo - tabla$viejo
# Marca las que cambian lo suficiente como para notarse al redondear en el texto
tabla$revisar <- ifelse(abs(tabla$diferencia) >= pmax(0.05 * abs(tabla$viejo), 0.01),
                        "<-- REVISAR", "")

cat("\n")
cat(sprintf("  %-34s %11s %11s\n", "Cantidad", "Vieja", "Nueva"))
cat("  ", strrep("-", 66), "\n", sep = "")
for (i in seq_len(nrow(tabla)))
  cat(sprintf("  %-34s %11.3f %11.3f  %s\n",
              tabla$cantidad[i], tabla$viejo[i], tabla$nuevo[i], tabla$revisar[i]))

write.csv(tabla, file.path(OUT, "consistencia_viejo_vs_nuevo.csv"), row.names = FALSE)
cat("\n  Guardado: consistencia_viejo_vs_nuevo.csv\n")
cat("  Las filas marcadas REVISAR cambian mas de un 5% y deben actualizarse.\n")
cat("  Las demas se pueden conservar redondeadas como estan.\n")
cat("\n=== LISTO ===\n")
