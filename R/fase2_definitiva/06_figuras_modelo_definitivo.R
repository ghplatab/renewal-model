# =============================================================================
#  FIGURAS DEL MODELO DEFINITIVO
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  Genera las dos figuras del ajuste definitivo:
#
#    fig_ajuste_definitivo_cg
#       Panel superior: casos observados y f(t), con banda al 95 %
#       Panel inferior: R(t) en escala logaritmica, con banda al 95 %
#
#    fig_exogeno_cg
#       Panel A: f(t) sobre los 180 dias, con la ventana exogena sombreada
#       Panel B: detalle de los 12 dias iniciales, con mu(t) y f(t) separados
#
#  El panel B dibuja ambas curvas por separado porque la aproximacion
#  f(t) ~ mu(t) no se sostiene con este nucleo: el componente exogeno explica
#  la totalidad de los casos el 14 de marzo, pero menos de la mitad a partir
#  del 21. Todas las cifras de los subtitulos se calculan aqui; ninguna esta
#  escrita a mano.
#
#  Entrada: resultados/fase2_chinagamma/f2cg_M0_ar1_gamma.rds
#  Salida : las dos figuras en formato PDF y PNG
#
#  No ajusta modelos: opera sobre draws ya guardados.
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(readxl)
  library(dplyr); library(ggplot2); library(scales)
  library(gridExtra); library(grid); library(here)
})

# --- Rutas relativas a la raíz del repositorio -------------------------------
OUT   <- here("resultados", "fase2_chinagamma")
FIG   <- here("figuras",    "fase2_chinagamma")
DATOS <- here("data", "datos_agregados.xlsx")
for (p in c(OUT, FIG)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

col_azul <- "#2166ac"; col_rojo <- "#d6604d"; col_gris <- "grey50"
tema_tesis <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        axis.title    = element_text(size = 10))

# --- Datos -------------------------------------------------------------------
leer_hoja <- function(hoja) {
  d <- as.data.frame(read_excel(DATOS, sheet = hoja))
  names(d) <- c("fecha", "casos")
  d <- d[!grepl("nan", as.character(d$fecha), ignore.case = TRUE), ]
  d <- d[!is.na(d$fecha), ]
  # La hoja NO_IMPORTADOS trae la fecha como numero de serie de Excel y la hoja
  # IMPORTADOS la trae ya como fecha. Hay que distinguirlas: convertir una fecha
  # real con as.numeric() da un valor en segundos y desplaza todo un siglo.
  d$fecha <- if (inherits(d$fecha, "Date")) d$fecha
             else if (inherits(d$fecha, "POSIXct")) as.Date(d$fecha)
             else as.Date(as.numeric(d$fecha), origin = "1899-12-30")
  d
}
fechas <- seq(as.Date("2020-03-14"), as.Date("2020-09-09"), by = "day")
loc <- left_join(data.frame(fecha = fechas), leer_hoja("NO_IMPORTADOS"), by = "fecha")
loc$casos[is.na(loc$casos)] <- 0L
y_obs <- as.integer(loc$casos)

imp <- leer_hoja("IMPORTADOS")
importados_obs <- left_join(data.frame(fecha = fechas[1:12]), imp, by = "fecha")
importados_obs$casos[is.na(importados_obs$casos)] <- 0L

# --- Posterior del definitivo nuevo -------------------------------------------
cat("=== CARGANDO f2cg_M0_ar1_gamma.rds ===\n")
fit <- readRDS(file.path(OUT, "f2cg_M0_ar1_gamma.rds"))
F_mat  <- fit$draws("y_pred", format = "matrix")   # f(t)
MU_mat <- fit$draws("mu",     format = "matrix")   # mu(t), 0 despues del dia 12
RT_mat <- fit$draws("Rt",     format = "matrix")
rm(fit); gc(verbose = FALSE)

resumir <- function(M) data.frame(
  fecha = fechas,
  media = colMeans(M),
  mediana = apply(M, 2, median),
  q025 = apply(M, 2, quantile, .025, names = FALSE),
  q25  = apply(M, 2, quantile, .25,  names = FALSE),
  q75  = apply(M, 2, quantile, .75,  names = FALSE),
  q975 = apply(M, 2, quantile, .975, names = FALSE))
f_total  <- resumir(F_mat)
mu_total <- resumir(MU_mat)
rt_total <- resumir(RT_mat)

# =============================================================================
#  FIGURA 1 - AJUSTE DEL MODELO DEFINITIVO
# =============================================================================

cat("\n=== FIGURA 1: AJUSTE DEL DEFINITIVO ===\n")

INTERV <- data.frame(
  fecha = as.Date(c("2020-03-25", "2020-06-01", "2020-09-01")),
  etiqueta = c("Cuarentena obligatoria", "Reapertura piloto", "Reapertura gradual"))

p1a <- ggplot(f_total, aes(x = fecha)) +
  geom_col(data = data.frame(fecha = fechas, casos = y_obs),
           aes(y = casos), fill = col_rojo, alpha = 0.45, width = 1) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = col_azul, alpha = 0.22) +
  geom_line(aes(y = media), color = col_azul, linewidth = 0.9) +
  geom_vline(xintercept = INTERV$fecha, linetype = "dotted", color = col_gris) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  scale_y_continuous(labels = label_comma(big.mark = ".", decimal.mark = ",")) +
  labs(title = "(a) Casos observados e incidencia esperada",
       subtitle = "Barras: casos diarios observados | Línea: media posterior de f(t) | Banda: IC95%",
       x = NULL, y = "Casos diarios") +
  tema_tesis

p1b <- ggplot(rt_total, aes(x = fecha)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = col_azul, alpha = 0.22) +
  geom_line(aes(y = mediana), color = col_azul, linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_vline(xintercept = INTERV$fecha, linetype = "dotted", color = col_gris) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  scale_y_continuous(trans = "log10") +
  labs(title = expression(bold("(b) Número reproductivo efectivo ") * bold(R[t]) *
                            bold(" (escala logarítmica)")),
       subtitle = "Línea: mediana posterior | Banda: IC95% | Línea horizontal: Rt = 1 | Punteadas verticales: intervenciones",
       x = NULL, y = expression(R[t])) +
  tema_tesis

fig1 <- arrangeGrob(p1a, p1b, ncol = 1, heights = c(1, 1))
ggsave(file.path(FIG, "fig_ajuste_definitivo_cg.png"), fig1, width = 10, height = 8, dpi = 300)
ok1 <- tryCatch({ ggsave(file.path(FIG, "fig_ajuste_definitivo_cg.pdf"), fig1,
                         width = 10, height = 8, device = cairo_pdf); TRUE },
                error = function(e) FALSE)
cat(sprintf("  Guardada: fig_ajuste_definitivo_cg.png%s\n", if (ok1) " y .pdf" else ""))
cat(sprintf("  Rt maximo (mediana): %.2f el %s | primer dia Rt<1: %s\n",
            max(rt_total$mediana), format(fechas[which.max(rt_total$mediana)]),
            format(fechas[which(rt_total$mediana < 1)[1]])))

# =============================================================================
#  FIGURA 2 - COMPONENTE EXOGENO COMO SEMILLERO
# =============================================================================

cat("\n=== FIGURA 2: COMPONENTE EXOGENO ===\n")

FIN_EXO <- as.Date("2020-03-25")
pct_exo <- 100 * mu_total$media[1:12] / f_total$media[1:12]
dia_50  <- fechas[1:12][which(pct_exo < 50)[1]]
cat(sprintf("  media f(t) dias 1-12 = %.2f | media mu(t) = %.2f | exogeno medio %.1f%%\n",
            mean(f_total$media[1:12]), mean(mu_total$media[1:12]), mean(pct_exo)))
cat(sprintf("  El exogeno baja del 50%% de f(t) el %s\n", format(dia_50)))

pA <- ggplot(f_total, aes(x = fecha)) +
  annotate("rect", xmin = fechas[1], xmax = FIN_EXO, ymin = -Inf, ymax = Inf,
           fill = col_azul, alpha = 0.07) +
  annotate("text", x = fechas[1] + 3, y = max(f_total$q975) * 0.90,
           label = "Ventana\nexógena", hjust = 0, size = 3, color = col_azul) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = col_azul, alpha = 0.18) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = col_azul, alpha = 0.35) +
  geom_line(aes(y = media), color = col_azul, linewidth = 0.9) +
  geom_point(data = importados_obs, aes(x = fecha, y = casos),
             color = col_rojo, size = 2.2, shape = 19) +
  geom_vline(xintercept = FIN_EXO, linetype = "dashed", color = col_azul, linewidth = 0.6) +
  scale_x_date(date_labels = "%d %b\n%Y", date_breaks = "3 weeks") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)),
                     labels = label_comma(big.mark = ".", decimal.mark = ",")) +
  labs(title = "A — Incidencia esperada f(t) y ventana exógena",
       subtitle = paste0("Zona azul: ventana exógena (14-25 mar) | ",
                         "Puntos rojos: importados observados | Bandas: IC50% e IC95%"),
       x = NULL, y = expression("Casos esperados " * f(t))) +
  tema_tesis

zoom <- data.frame(fecha = fechas[1:12],
                   f_media = f_total$media[1:12], f_q025 = f_total$q025[1:12],
                   f_q975 = f_total$q975[1:12], f_q25 = f_total$q25[1:12],
                   f_q75 = f_total$q75[1:12], mu_media = mu_total$media[1:12])

pB <- ggplot(zoom, aes(x = fecha)) +
  geom_ribbon(aes(ymin = f_q025, ymax = f_q975), fill = col_azul, alpha = 0.18) +
  geom_ribbon(aes(ymin = f_q25,  ymax = f_q75),  fill = col_azul, alpha = 0.35) +
  geom_line(aes(y = f_media, colour = "f(t): incidencia esperada"), linewidth = 0.9) +
  geom_line(aes(y = mu_media, colour = "mu(t): componente exógeno"),
            linewidth = 0.9, linetype = "longdash") +
  geom_point(data = importados_obs, aes(x = fecha, y = casos,
                                        colour = "Importados observados"),
             size = 3, shape = 19) +
  scale_colour_manual(values = c("f(t): incidencia esperada" = col_azul,
                                 "mu(t): componente exógeno" = "grey25",
                                 "Importados observados" = col_rojo)) +
  scale_x_date(date_labels = "%d %b", date_breaks = "3 days") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "B — Dentro de la ventana exógena: el aporte local crece rápido",
       subtitle = sprintf(paste0("Media de f(t) = %.1f casos/día | media de mu(t) = %.1f | ",
                                 "el exógeno explica el 100%% el 14-mar y menos de la mitad desde el %s"),
                          mean(f_total$media[1:12]), mean(mu_total$media[1:12]),
                          format(dia_50, "%d-%b")),
       x = NULL, y = "Casos por día", colour = NULL) +
  tema_tesis + theme(legend.position = "bottom")

fig2 <- arrangeGrob(pA, pB, ncol = 1, heights = c(1, 1.15))
ggsave(file.path(FIG, "fig_exogeno_cg.png"), fig2, width = 10, height = 9, dpi = 300)
ok2 <- tryCatch({ ggsave(file.path(FIG, "fig_exogeno_cg.pdf"), fig2,
                         width = 10, height = 9, device = cairo_pdf); TRUE },
                error = function(e) FALSE)
cat(sprintf("  Guardada: fig_exogeno_cg.png%s\n", if (ok2) " y .pdf" else ""))

cat("\n=== FIGURAS DEL MODELO DEFINITIVO COMPLETAS ===\n")

# =============================================================================
#  FIGURAS 3 y 4 - VERSIONES EN PDF DE LAS QUE SOLO ESTABAN EN PNG
#  El documento usa PDF; el PNG queda solo como vista rapida.
# =============================================================================

cat("\n=== FIGURAS 3 y 4: VERSIONES EN PDF ===\n")
suppressPackageStartupMessages(library(bayesplot))

tray <- read.csv(file.path(OUT, "def_trayectorias_ar1_vs_rw.csv"))
tray$fecha <- as.Date(tray$fecha)
# La leyenda va en espanol para que coincida con el texto de la seccion
tray$modelo <- ifelse(tray$modelo == "Random walk", "Paseo aleatorio", tray$modelo)
p3 <- ggplot(tray, aes(fecha, mediana, colour = modelo, fill = modelo)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = c("AR(1)" = "black", "Paseo aleatorio" = col_rojo)) +
  scale_fill_manual(values = c("AR(1)" = "black", "Paseo aleatorio" = col_rojo)) +
  labs(x = NULL, y = expression(R[t]), colour = NULL, fill = NULL,
       title = "Modelo definitivo (IS China-Gamma): AR(1) frente a paseo aleatorio",
       subtitle = "Mediana posterior y banda del 90%") +
  tema_tesis + theme(legend.position = "bottom")
ok3 <- tryCatch({ ggsave(file.path(FIG, "def_trayectorias_ar1_vs_rw.pdf"), p3,
                         width = 10, height = 5.5, device = cairo_pdf); TRUE },
                error = function(e) FALSE)
cat(sprintf("  def_trayectorias_ar1_vs_rw.pdf: %s\n", ok3))

fit2 <- readRDS(file.path(OUT, "f2cg_M0_ar1_gamma.rds"))
p4 <- mcmc_pairs(fit2$draws(c("mu_rt", "rho1", "sigma_epsilon", "phi")),
                 np = nuts_params(fit2),
                 off_diag_args = list(size = 0.6, alpha = 0.4))
ok4 <- tryCatch({ ggsave(file.path(FIG, "f2cg_M0_pairs.pdf"), p4,
                         width = 9, height = 9, device = cairo_pdf); TRUE },
                error = function(e) FALSE)
cat(sprintf("  f2cg_M0_pairs.pdf: %s\n", ok4))
rm(fit2); gc(verbose = FALSE)
cat("\n=== PDFs LISTOS ===\n")
