# =============================================================================
#  COBERTURA DE LOS INTERVALOS PREDICTIVOS
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  Calcula la cobertura empirica de los intervalos predictivos posteriores
#  del modelo definitivo, contrastando dos construcciones:
#
#    - a partir de y_pred, la media de la distribucion predictiva
#    - a partir de y_rep,  replicas simuladas de la distribucion predictiva
#
#  La segunda es la correcta: incorpora la variabilidad de la verosimilitud
#  NB-2 ademas de la incertidumbre sobre los parametros.
#
#  No ajusta modelos: opera sobre draws ya guardados.
#  Ejecutar primero el bloque de configuracion de rutas.
# =============================================================================


# =============================================================================
#  PASO 0 - CONFIGURACION
# =============================================================================

# --- Paquetes (todos ya instalados) ---
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# --- Raiz del proyecto -------------------------------------------------------
# Si esta linea falla por el acento en "Maestria", en RStudio usa
#   Session > Set Working Directory > Choose Directory...
# y apunta a la carpeta AJUSTES_PROFE_JCS. Luego pon: BASE <- getwd()
BASE <- file.path(Sys.getenv("USERPROFILE"),
                  "OneDrive", "Escritorio", "Tesis Maestría",
                  "AJUSTES_2026", "AJUSTES_PROFE_JCS")

STAN_ROOT <- file.path(BASE, "Documentos George", "MODELO_STAN")

# --- Las DOS copias del mismo archivo, que difieren entre si ------------------
RDS_A <- file.path(STAN_ROOT, "MODELO_REPO_GITHUB", "MODELO_STAN_DEFINITIVO",
                   "draws_AR1_Gamma_definitivo.rds")   # 37.393.962 bytes
RDS_B <- file.path(STAN_ROOT, "MODELO_STAN_DEFINITIVO",
                   "draws_AR1_Gamma_definitivo.rds")   # 35.168.812 bytes

# --- Datos observados --------------------------------------------------------
DATOS <- file.path(STAN_ROOT, "MODELO_STAN_4", "datos_agregados.xlsx")

# --- Salidas (carpeta nueva, no toca nada existente) -------------------------
OUT <- file.path(BASE, "AJUSTES_SEP_2026", "resultados")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- Verificacion ------------------------------------------------------------
cat("=== PASO 0: VERIFICACION DE RUTAS ===\n")
for (nm in c("RDS_A", "RDS_B", "DATOS")) {
  p  <- get(nm)
  ok <- file.exists(p)
  cat(sprintf("  %-6s %-3s %s\n", nm, if (ok) "OK" else "NO", basename(p)))
  if (ok) cat(sprintf("         %s bytes\n", format(file.size(p), big.mark = ".")))
}
cat(sprintf("  SALIDA %-3s %s\n", if (dir.exists(OUT)) "OK" else "NO", OUT))
cat("\nSi todo dice OK, sigue al PASO 1.\n")


# =============================================================================
#  PASO 1 - CUAL DE LAS DOS COPIAS ES LA CANONICA
#
#  Por que: hay dos archivos con el mismo nombre y contenido distinto (md5
#  diferente, 2,2 MB de diferencia). Si abrimos el equivocado, TODO lo que
#  calculemos despues esta mal y no nos enteramos.
#
#  Como decidimos: comparamos las medias posteriores de los cuatro parametros
#  escalares contra lo que reporta la Tabla 3-13 del documento. La copia que
#  coincida es la que genero los numeros publicados.
# =============================================================================

# Valores de referencia del documento (Tabla 3-13 / plan seccion 3.1)
REFERENCIA <- c(mu_rt = 0.203, rho1 = 0.962, sigma_epsilon = 0.099, phi = 3.263)

# Funcion que abre un .rds, lo resume, y lo descarga de memoria enseguida.
# Descargarlo importa: cada uno ocupa varios cientos de MB al expandirse y
# tienes 8 GB de RAM.
resumir_rds <- function(ruta, etiqueta) {
  cat(sprintf("\n--- %s ---\n%s\n", etiqueta, ruta))
  d <- readRDS(ruta)

  cat("Elementos guardados:", paste(names(d), collapse = ", "), "\n")

  # Dimensiones de las matrices de draws (filas = draws, columnas = dias)
  for (v in c("y_pred", "y_rep", "log_lik", "Rt", "f", "log_Rt")) {
    if (!is.null(d[[v]])) {
      dm <- dim(as.matrix(d[[v]]))
      cat(sprintf("  %-8s %5d draws x %3d dias\n", v, dm[1], dm[2]))
    } else {
      cat(sprintf("  %-8s AUSENTE\n", v))
    }
  }

  # Medias posteriores de los escalares
  esc <- d$escalares
  cat("\n  Parametro        media    documento   diferencia\n")
  medias <- sapply(names(REFERENCIA), function(p) {
    v <- if (!is.null(esc[[p]])) mean(as.vector(esc[[p]])) else NA_real_
    cat(sprintf("  %-14s %8.4f %10.3f %12.4f\n",
                p, v, REFERENCIA[[p]], v - REFERENCIA[[p]]))
    v
  })

  # Distancia total a los valores publicados: cuanto mas chica, mas probable
  # que sea la copia que genero el documento
  dist <- sqrt(sum((medias - REFERENCIA)^2, na.rm = TRUE))
  cat(sprintf("\n  Distancia a los valores publicados: %.5f\n", dist))

  # Diagnosticos que el script guardo junto a los draws
  if (!is.null(d$diagnostics)) {
    dg <- d$diagnostics
    cat(sprintf("  Diagnosticos guardados: rhat=%s  ess=%s  divergencias=%s\n",
                signif(dg$rhat, 5), signif(dg$ess, 5), dg$n_div))
  }
  if (!is.null(d$meta)) {
    cat(sprintf("  Meta: modelo=%s  n_draws=%s  fecha=%s\n",
                d$meta$modelo, d$meta$n_draws, as.character(d$meta$fecha)))
  }

  rm(d); gc(verbose = FALSE)
  invisible(list(medias = medias, dist = dist))
}

cat("\n=== PASO 1: INVENTARIO DE LAS DOS COPIAS ===\n")
resA <- resumir_rds(RDS_A, "COPIA A - MODELO_REPO_GITHUB")
resB <- resumir_rds(RDS_B, "COPIA B - MODELO_STAN_DEFINITIVO")

cat("\n=== VEREDICTO ===\n")
cat(sprintf("  Copia A: distancia %.5f\n", resA$dist))
cat(sprintf("  Copia B: distancia %.5f\n", resB$dist))
CANONICA <- if (resA$dist <= resB$dist) RDS_A else RDS_B
cat(sprintf("  --> Usaremos: %s\n",
            if (identical(CANONICA, RDS_A)) "COPIA A" else "COPIA B"))
cat("\nPasame esta salida antes de seguir al PASO 2.\n")


# =============================================================================
#  PASO 2 - COBERTURA POR PARTIDA DOBLE
#
#  El problema: el script PPC_CALIBRADA_TEMPORAL.R hace en su linea 36
#      y_rep <- as.matrix(draws$y_pred)
#  La variable se LLAMA y_rep pero se le asigna y_pred. Y en el modelo Stan,
#      y_pred[t] = f[t]          <- la media determinista, la curva suave
#      y_rep[t]  = neg_binomial_2_rng(f[t], phi)   <- el predictivo de verdad
#
#  Consecuencia: las bandas publicadas son intervalos de credibilidad de la
#  MEDIA, no intervalos de prediccion del CONTEO. Les falta todo el ruido de
#  observacion de la binomial negativa.
#
#  Este bloque calcula la cobertura de las dos formas:
#    (a) con y_pred -> cobertura de los intervalos de credibilidad de la media
#    (b) con y_rep  -> cobertura predictiva, que es la cifra correcta
#
#  Si (b) se aproxima a 50 / 80 / 95, la discrepancia era un artefacto del
#  calculo. Si (b) permanece por debajo, se trata de descalibracion real del
#  modelo de observacion.
# =============================================================================

cat("\n=== PASO 2: COBERTURA ===\n")

# --- Datos observados, replicando exactamente el script original -------------
df_raw <- as.data.frame(read_excel(DATOS, sheet = "NO_IMPORTADOS"))
names(df_raw) <- c("fecha", "casos")
df_raw <- df_raw[!grepl("nan", as.character(df_raw$fecha), ignore.case = TRUE), ]
df_raw <- df_raw[!is.na(df_raw$fecha), ]
df_raw$fecha <- as.Date(as.numeric(df_raw$fecha), origin = "1899-12-30")

FECHA_INICIO <- as.Date("2020-03-14")
FECHA_FIN    <- as.Date("2020-09-09")
grilla <- data.frame(fecha = seq(FECHA_INICIO, FECHA_FIN, by = "day"))
df     <- left_join(grilla, df_raw, by = "fecha")
df$casos[is.na(df$casos)] <- 0L
y_obs  <- as.integer(df$casos)

cat(sprintf("Datos: %d dias, %s casos totales, media %.1f/dia\n",
            length(y_obs), format(sum(y_obs), big.mark = "."), mean(y_obs)))

# --- Cargar la copia canonica ------------------------------------------------
d      <- readRDS(CANONICA)
M_pred <- as.matrix(d$y_pred)   # f(t): la media. Lo que uso el script original
M_rep  <- as.matrix(d$y_rep)    # el predictivo posterior de verdad
rm(d); gc(verbose = FALSE)

stopifnot(ncol(M_pred) == length(y_obs), ncol(M_rep) == length(y_obs))

# --- Funcion de cobertura ----------------------------------------------------
# Para cada dia construye el intervalo central de nivel p a partir de los
# cuantiles de la matriz, y mide en que fraccion de dias el valor observado
# cae dentro. Con 180 dias y nivel 95%, lo esperable son ~171 dias dentro.
NIVELES <- c(0.50, 0.80, 0.95)

cobertura <- function(M, y, niveles = NIVELES) {
  sapply(niveles, function(p) {
    a  <- (1 - p) / 2
    lo <- apply(M, 2, quantile, probs = a,     names = FALSE)
    hi <- apply(M, 2, quantile, probs = 1 - a, names = FALSE)
    mean(y >= lo & y <= hi)
  })
}

# Ancho medio de la banda: sirve para ver DE CUANTO es el error, no solo que lo hay
ancho_medio <- function(M, p = 0.95) {
  a  <- (1 - p) / 2
  lo <- apply(M, 2, quantile, probs = a,     names = FALSE)
  hi <- apply(M, 2, quantile, probs = 1 - a, names = FALSE)
  mean(hi - lo)
}

cob_pred <- cobertura(M_pred, y_obs)
cob_rep  <- cobertura(M_rep,  y_obs)

# --- Resultados --------------------------------------------------------------
cat("\n  Nivel   Publicado   Con y_pred   Con y_rep\n")
cat("  ------------------------------------------------\n")
publicado <- c(29.4, 46.7, 68.9)
for (i in seq_along(NIVELES)) {
  cat(sprintf("  %3.0f%%    %6.1f%%    %8.1f%%   %8.1f%%\n",
              NIVELES[i] * 100, publicado[i],
              cob_pred[i] * 100, cob_rep[i] * 100))
}

cat(sprintf("\n  Ancho medio de la banda 95%%:\n"))
cat(sprintf("    con y_pred : %8.1f casos\n", ancho_medio(M_pred)))
cat(sprintf("    con y_rep  : %8.1f casos\n", ancho_medio(M_rep)))
cat(sprintf("    factor     : %8.1f veces mas ancha la correcta\n",
            ancho_medio(M_rep) / ancho_medio(M_pred)))

cat(sprintf("\n  Dias fuera de la banda 95%% (de %d):\n", length(y_obs)))
cat(sprintf("    con y_pred : %3d dias\n", round((1 - cob_pred[3]) * length(y_obs))))
cat(sprintf("    con y_rep  : %3d dias\n", round((1 - cob_rep[3])  * length(y_obs))))

# --- Guardar -----------------------------------------------------------------
tabla <- data.frame(
  nivel_nominal   = NIVELES,
  publicado_tesis = publicado / 100,
  con_y_pred      = as.numeric(cob_pred),
  con_y_rep       = as.numeric(cob_rep)
)
write.csv(tabla, file.path(OUT, "paso2_cobertura.csv"), row.names = FALSE)
cat(sprintf("\nGuardado: %s\n", file.path(OUT, "paso2_cobertura.csv")))

cat("\n=== COMO LEER ESTO ===\n")
cat("  - Si 'Con y_pred' reproduce 29,4 / 46,7 / 68,9 -> confirmado el bug.\n")
cat("  - Si 'Con y_rep' se acerca a 50 / 80 / 95      -> obs. 8 resuelta.\n")
cat("  - Si 'Con y_rep' sigue muy por debajo          -> descalibracion real.\n")
cat("\nPasame la salida y seguimos.\n")
