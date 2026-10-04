# =============================================================================
#  COMPARACION DE TRAYECTORIAS DE R(t): AR(1) FRENTE A PASEO ALEATORIO
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  Las dos especificaciones son predictivamente indistinguibles, de modo que
#  mu_R no queda identificado por los datos. Esa es, sin embargo, una
#  afirmacion sobre un parametro. Este script evalua la pregunta mas amplia:
#
#      cuanto difieren las TRAYECTORIAS estimadas de R(t) entre ambas?
#
#  Compara las medianas diarias, la diferencia media y maxima entre ambas
#  trayectorias, el valor maximo de R(t) y la fecha en que cruza el umbral
#  de uno. De ello depende que conclusiones del capitulo son robustas a la
#  eleccion del proceso latente y cuales dependen de ella.
#
#  No ajusta modelos: opera sobre los dos ajustes ya guardados.
# =============================================================================


# =============================================================================
#  CONFIGURACION Y CARGA DE DATOS
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(ggplot2)
  library(dplyr); library(tidyr); library(readxl); library(here)
})

# --- Rutas relativas a la raíz del repositorio -------------------------------
OUT <- here("resultados"); FIG <- here("figuras")
for (p in c(OUT, FIG)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

fechas <- seq(as.Date("2020-03-14"), as.Date("2020-09-09"), by = "day")
T_P    <- length(fechas)

cat("=== CARGANDO AJUSTES ===\n")
fit_ar1 <- readRDS(file.path(OUT, "corrida1_ar1_ad099.rds"))
Rt_ar1  <- fit_ar1$draws("Rt", format = "matrix")   # 6000 x 180
rm(fit_ar1); gc(verbose = FALSE)

fit_rw <- readRDS(file.path(OUT, "corrida3_rw.rds"))
Rt_rw  <- fit_rw$draws("Rt", format = "matrix")
rm(fit_rw); gc(verbose = FALSE)

cat(sprintf("  AR(1): %d draws x %d dias\n", nrow(Rt_ar1), ncol(Rt_ar1)))
cat(sprintf("  RW   : %d draws x %d dias\n", nrow(Rt_rw),  ncol(Rt_rw)))
stopifnot(ncol(Rt_ar1) == T_P, ncol(Rt_rw) == T_P)

# Resumen por dia: mediana y banda del 90%
resumir <- function(M, etiqueta) {
  data.frame(
    fecha   = fechas,
    modelo  = etiqueta,
    mediana = apply(M, 2, median),
    lo      = apply(M, 2, quantile, 0.05, names = FALSE),
    hi      = apply(M, 2, quantile, 0.95, names = FALSE),
    p_menor1 = apply(M, 2, function(v) mean(v < 1))   # P(R(t) < 1)
  )
}
s_ar1 <- resumir(Rt_ar1, "AR(1)")
s_rw  <- resumir(Rt_rw,  "Random walk")


# =============================================================================
#  PASO 1 - ACUERDO ENTRE LAS TRAYECTORIAS
# =============================================================================

cat("\n=== PASO 1: ACUERDO ENTRE TRAYECTORIAS ===\n")

dif <- s_rw$mediana - s_ar1$mediana

cat(sprintf("  Correlacion entre medianas    : %.4f\n",
            cor(s_ar1$mediana, s_rw$mediana)))
cat(sprintf("  Diferencia absoluta media     : %.4f\n", mean(abs(dif))))
cat(sprintf("  Diferencia absoluta mediana   : %.4f\n", median(abs(dif))))
cat(sprintf("  Diferencia absoluta maxima    : %.4f  (%s)\n",
            max(abs(dif)), format(fechas[which.max(abs(dif))])))
cat(sprintf("  Sesgo medio (RW - AR1)        : %+.4f\n", mean(dif)))

cat("\n  Dias con diferencia mayor que:\n")
for (u in c(0.05, 0.10, 0.20, 0.43))
  cat(sprintf("    %.2f : %3d de %d  (%.1f%%)\n",
              u, sum(abs(dif) > u), T_P, 100 * mean(abs(dif) > u)))
cat("  (0.43 es el rango promedio diario de R(t) entre especificaciones del\n")
cat("   intervalo serial, incluido como referencia de magnitud)\n")


# =============================================================================
#  PASO 2 - ANCHURA DE LOS INTERVALOS
#
#  El random walk no tiene ancla: nada lo devuelve a un nivel. Cabe esperar
#  que sus intervalos sean mas anchos, sobre todo al final de la serie.
#  Si lo son mucho, es una diferencia sustantiva aunque las medianas
#  coincidan.
# =============================================================================

cat("\n=== PASO 2: ANCHURA DE LOS INTERVALOS (90%) ===\n")

w_ar1 <- s_ar1$hi - s_ar1$lo
w_rw  <- s_rw$hi  - s_rw$lo

cat(sprintf("  Ancho medio  AR(1) = %.3f | RW = %.3f | razon = %.2f\n",
            mean(w_ar1), mean(w_rw), mean(w_rw) / mean(w_ar1)))

cat("\n  Ancho medio por tercio del periodo:\n")
for (ix in list(1:60, 61:120, 121:180))
  cat(sprintf("    dias %3d-%3d : AR(1) = %.3f | RW = %.3f | razon = %.2f\n",
              min(ix), max(ix), mean(w_ar1[ix]), mean(w_rw[ix]),
              mean(w_rw[ix]) / mean(w_ar1[ix])))
cat("  (si la razon crece hacia el final, el RW se desancla con el tiempo)\n")

# Solapamiento: en cuantos dias la mediana de un modelo cae dentro de la
# banda del otro. Es la prueba practica de si son compatibles.
dentro_ar1 <- mean(s_ar1$mediana >= s_rw$lo  & s_ar1$mediana <= s_rw$hi)
dentro_rw  <- mean(s_rw$mediana  >= s_ar1$lo & s_rw$mediana  <= s_ar1$hi)
cat(sprintf("\n  Mediana AR(1) dentro de la banda del RW : %.1f%% de los dias\n",
            100 * dentro_ar1))
cat(sprintf("  Mediana RW dentro de la banda del AR(1) : %.1f%% de los dias\n",
            100 * dentro_rw))


# =============================================================================
#  PASO 3 - CANTIDADES EPIDEMIOLOGICAS EN FECHAS DE REFERENCIA
#
#  El 25 de marzo corresponde a la entrada en vigor de la cuarentena
#  nacional obligatoria (Decreto 457 de 2020); las restantes son marcas
#  mensuales de referencia. La coincidencia de R(t) en estas fechas entre
#  ambas especificaciones indica en que medida las conclusiones sobre la
#  trayectoria dependen del proceso latente elegido.
# =============================================================================

cat("\n=== R(t) EN FECHAS DE REFERENCIA ===\n")

hitos <- as.Date(c("2020-03-14", "2020-03-25", "2020-04-15", "2020-05-15",
                   "2020-06-15", "2020-07-15", "2020-08-15", "2020-09-09"))
etiquetas <- c("inicio serie", "Decreto 457", "abril", "mayo",
               "junio", "julio", "agosto", "fin serie")

cat("\n  Fecha        Hito            AR(1)  [IC90]         RW  [IC90]\n")
cat("  ", strrep("-", 68), "\n", sep = "")
for (i in seq_along(hitos)) {
  j <- which(fechas == hitos[i])
  if (!length(j)) next
  cat(sprintf("  %s  %-14s  %.2f [%.2f,%.2f]   %.2f [%.2f,%.2f]\n",
              format(hitos[i]), etiquetas[i],
              s_ar1$mediana[j], s_ar1$lo[j], s_ar1$hi[j],
              s_rw$mediana[j],  s_rw$lo[j],  s_rw$hi[j]))
}

# Primer cruce por debajo de 1: el hito epidemiologico principal
cruce <- function(s, etiqueta) {
  j <- which(s$mediana < 1)
  if (!length(j)) { cat(sprintf("  %-12s nunca cruza 1\n", etiqueta)); return(NA) }
  cat(sprintf("  %-12s primer dia con mediana < 1: %s (dia %d)\n",
              etiqueta, format(s$fecha[j[1]]), j[1]))
  j[1]
}
cat("\n  Cruce de R(t) por debajo de 1:\n")
c1 <- cruce(s_ar1, "AR(1)"); c2 <- cruce(s_rw, "Random walk")
if (!is.na(c1) && !is.na(c2))
  cat(sprintf("  Diferencia entre modelos: %d dias\n", abs(c2 - c1)))

# Probabilidad de que R(t) < 1 al final del periodo
cat(sprintf("\n  P(R(t) < 1) el ultimo dia: AR(1) = %.3f | RW = %.3f\n",
            s_ar1$p_menor1[T_P], s_rw$p_menor1[T_P]))


# =============================================================================
#  PASO 4 - FIGURA
# =============================================================================

cat("\n=== PASO 4: FIGURA ===\n")

comb <- rbind(s_ar1, s_rw)

p_tray <- ggplot(comb, aes(fecha, mediana, color = modelo, fill = modelo)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = c("AR(1)" = "#1F4E79", "Random walk" = "#C1440E")) +
  scale_fill_manual(values  = c("AR(1)" = "#1F4E79", "Random walk" = "#C1440E")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(title = "Trayectorias de R(t): AR(1) estacionario vs random walk",
       subtitle = "Linea: mediana posterior | Banda: intervalo de credibilidad del 90%",
       x = NULL, y = "R(t)", color = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"))

p_dif <- ggplot(data.frame(fecha = fechas, dif = dif), aes(fecha, dif)) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_line(linewidth = 0.6, color = "#444444") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(title = "Diferencia entre medianas (Random walk menos AR(1))",
       x = NULL, y = "Diferencia en R(t)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave(file.path(FIG, "trayectorias_ar1_vs_rw.png"), p_tray,
       width = 9, height = 4.5, dpi = 200)
ggsave(file.path(FIG, "trayectorias_diferencia.png"), p_dif,
       width = 9, height = 3, dpi = 200)
cat("  Guardadas: trayectorias_ar1_vs_rw.png y trayectorias_diferencia.png\n")


# =============================================================================
#  PASO 5 - GUARDAR
# =============================================================================

tab <- data.frame(
  fecha        = fechas,
  ar1_mediana  = s_ar1$mediana, ar1_lo = s_ar1$lo, ar1_hi = s_ar1$hi,
  ar1_p_menor1 = s_ar1$p_menor1,
  rw_mediana   = s_rw$mediana,  rw_lo  = s_rw$lo,  rw_hi  = s_rw$hi,
  rw_p_menor1  = s_rw$p_menor1,
  diferencia   = dif)
write.csv(tab, file.path(OUT, "trayectorias_ar1_vs_rw.csv"), row.names = FALSE)
cat("  Guardado: trayectorias_ar1_vs_rw.csv (las 180 filas, por si tu\n")
cat("            Tabla 3-14 usa fechas distintas a las de arriba)\n")

cat("\n=== COMO LEER ESTO ===\n")
cat("  Diferencia media < 0.05 y correlacion > 0.99\n")
cat("     -> las trayectorias son la misma. La Seccion 3.2.9 sobrevive.\n")
cat("  Diferencia media > 0.20, o el cruce de 1 se mueve mas de una semana\n")
cat("     -> hay que replantear el analisis de intervenciones.\n")
