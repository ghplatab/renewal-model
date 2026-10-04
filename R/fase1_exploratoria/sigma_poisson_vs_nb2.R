# =============================================================================
#  ESCALA DE LAS INNOVACIONES BAJO POISSON Y BAJO NB-2
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  Contrasta la desviacion estandar posterior de las innovaciones diarias,
#  sigma_epsilon, en las diez especificaciones preliminares: cinco intervalos
#  seriales por dos distribuciones de observacion.
#
#  La verosimilitud Poisson carece de parametro de dispersion, de modo que la
#  variabilidad de los conteos que excede la media solo puede ser absorbida
#  por el proceso latente. Si el mecanismo opera, sigma_epsilon deberia
#  resultar sistematicamente mayor bajo Poisson, en conflicto con su previa
#  N+(0, 0.2).
#
#  Entrada: resultados/fase1_exploratoria/draws_gamma_*.rds, generados por
#           analisis_bogota_ar2_gamma.R
#  Salida : resultados/fase1_exploratoria/sigma_poisson_vs_nb2.csv
#
#  No ajusta modelos: opera sobre draws ya guardados. Las especificaciones
#  preliminares son AR(2), de modo que sus escalares incluyen rho2.
# =============================================================================

library(here)

RUTA_OUTPUT <- here("resultados", "fase1_exploratoria")
if (!dir.exists(RUTA_OUTPUT)) dir.create(RUTA_OUTPUT, recursive = TRUE)

archivos <- list.files(RUTA_OUTPUT, pattern = "^draws_gamma_.*\\.rds$",
                       full.names = TRUE)
cat(sprintf("Ajustes preliminares encontrados: %d\n\n", length(archivos)))

filas <- list()
for (a in archivos) {
  nm  <- sub("^draws_gamma_", "", tools::file_path_sans_ext(basename(a)))
  dd  <- readRDS(a)
  e   <- dd$escalares
  s   <- if (!is.null(e$sigma_epsilon)) as.vector(e$sigma_epsilon) else NA_real_
  div <- if (!is.null(dd$diagnostics)) dd$diagnostics$n_div else NA
  rh  <- if (!is.null(dd$diagnostics)) dd$diagnostics$rhat  else NA
  cat(sprintf("  %-26s sigma = %.4f (sd %.4f)  rhat = %s  div = %s\n",
              nm, mean(s), sd(s), signif(rh, 5), div))
  filas[[nm]] <- data.frame(
    modelo        = nm,
    verosimilitud = ifelse(grepl("Poisson", nm), "Poisson", "NB-2"),
    sigma_media   = mean(s),
    sigma_sd      = sd(s),
    rhat          = as.numeric(rh)[1],
    n_div         = as.numeric(div)[1]
  )
  rm(dd); gc(verbose = FALSE)
}

tabla <- do.call(rbind, filas)

cat("\nPromedio por distribucion de observacion:\n")
print(aggregate(sigma_media ~ verosimilitud, tabla, mean))

write.csv(tabla, file.path(RUTA_OUTPUT, "sigma_poisson_vs_nb2.csv"),
          row.names = FALSE)
cat(sprintf("\nGuardado: %s\n",
            file.path(RUTA_OUTPUT, "sigma_poisson_vs_nb2.csv")))
