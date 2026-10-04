# =============================================================================
#  ANALISIS EXPLORATORIO DE INTERVENCIONES (MATERIAL DE RESPALDO)
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  Evalua el cambio en R(t) asociado a cada intervencion no farmaceutica,
#  comparando el promedio de R(t) en una ventana de 14 dias antes y otra
#  despues de cada medida, con un rezago intermedio excluido y la diferencia
#  calculada muestra a muestra. El calculo se repite con rezagos de 7, 10 y
#  14 dias, y sobre dos especificaciones del nucleo.
#
#  Los resultados no se reportan entre los hallazgos del trabajo, por cinco
#  razones que este script documenta y que se detallan en
#  resultados/respaldo_intervenciones/HALLAZGOS.md:
#
#    1. Ninguna intervencion produce un cambio distinguible de cero.
#    2. El signo del cambio depende del rezago elegido.
#    3. La ventana previa a la cuarentena cae en el periodo semilla, donde
#       R(t) apenas esta informado por los datos.
#    4. El modelo no incorpora las intervenciones como covariables y suaviza
#       los cambios bruscos por construccion.
#    5. Ninguno de los objetivos del trabajo requiere esta estimacion.
#
#  El analisis se conserva aqui como respaldo verificable de esa decision.
#
#  Salida: resultados/respaldo_intervenciones/*.csv
#  No ajusta modelos: opera sobre draws ya guardados.
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(here)
})
set_cmdstan_path(Sys.getenv("CMDSTAN", unset = "C:/cmdstan/cmdstan-2.39.0"))

# --- Rutas relativas a la raíz del repositorio -------------------------------
RES <- here("resultados")
OUT <- here("resultados", "respaldo_intervenciones")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

fechas <- seq(as.Date("2020-03-14"), by = "day", length.out = 180)

# Ventanas exactas de la Tabla 3-14 de la tesis (rezago de 10 dias)
VENTANAS <- data.frame(
  intervencion  = c("Cuarentena obligatoria", "Apertura sectorial", "Reapertura piloto"),
  fecha         = as.Date(c("2020-03-25", "2020-05-04", "2020-06-01")),
  antes_ini     = as.Date(c("2020-03-14", "2020-04-20", "2020-05-18")),
  antes_fin     = as.Date(c("2020-03-24", "2020-05-03", "2020-05-31")),
  despues_ini   = as.Date(c("2020-04-04", "2020-05-14", "2020-06-11")),
  despues_fin   = as.Date(c("2020-04-17", "2020-05-27", "2020-06-24")))

prom <- function(R, a, b) rowMeans(R[, fechas >= a & fechas <= b, drop = FALSE])
resumen <- function(D, antes, despues) data.frame(
  Rt_antes = mean(antes), Rt_despues = mean(despues),
  cambio = mean(D), ic_inf = quantile(D, .025, names = FALSE),
  ic_sup = quantile(D, .975, names = FALSE), P_reduccion = mean(D < 0))

analizar <- function(R, modelo) {
  tesis <- do.call(rbind, lapply(seq_len(nrow(VENTANAS)), function(i) {
    v <- VENTANAS[i, ]
    a <- prom(R, v$antes_ini, v$antes_fin); d <- prom(R, v$despues_ini, v$despues_fin)
    cbind(modelo = modelo, intervencion = v$intervencion, rezago = 10,
          ventana_antes = paste(format(v$antes_ini), format(v$antes_fin)),
          ventana_despues = paste(format(v$despues_ini), format(v$despues_fin)),
          resumen(d - a, a, d))
  }))
  sens <- do.call(rbind, lapply(seq_len(nrow(VENTANAS)), function(i) {
    f0 <- VENTANAS$fecha[i]
    do.call(rbind, lapply(c(7, 10, 14), function(lag) {
      a <- prom(R, max(f0 - 14, fechas[1]), f0 - 1)
      d <- prom(R, f0 + lag, f0 + lag + 13)
      cbind(modelo = modelo, intervencion = VENTANAS$intervencion[i], rezago = lag,
            resumen(d - a, a, d))
    }))
  }))
  list(tesis = tesis, sens = sens)
}

cat("=== Modelo viejo (nucleo 10.5 d) ===\n")
fit <- readRDS(file.path(RES, "corrida1_ar1_ad099.rds"))
viejo <- analizar(fit$draws("Rt", format = "matrix"), "viejo_IS_10.5d")
rm(fit); gc(verbose = FALSE)

cat("=== Modelo nuevo (China-Gamma) ===\n")
fit <- readRDS(file.path(RES, "fase2_chinagamma", "f2cg_M0_ar1_gamma.rds"))
nuevo <- analizar(fit$draws("Rt", format = "matrix"), "nuevo_ChinaGamma")
rm(fit); gc(verbose = FALSE)

tab_tesis <- rbind(viejo$tesis, nuevo$tesis)
tab_sens  <- rbind(viejo$sens,  nuevo$sens)
print(tab_tesis, digits = 3, row.names = FALSE)
print(tab_sens,  digits = 3, row.names = FALSE)
write.csv(tab_tesis, file.path(OUT, "intervenciones_ventanas_tesis.csv"), row.names = FALSE)
write.csv(tab_sens,  file.path(OUT, "intervenciones_sensibilidad_rezago.csv"), row.names = FALSE)
cat("\nGuardado en resultados/respaldo_intervenciones/\n")
