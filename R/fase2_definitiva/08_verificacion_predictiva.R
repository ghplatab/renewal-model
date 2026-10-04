# =============================================================================
#  VERIFICACION PREDICTIVA POSTERIOR - ANALISIS COMPLEMENTARIOS
#  Modelo bayesiano de renovacion - COVID-19 Bogota
#
#  Analisis adicionales sobre la calibracion del modelo definitivo:
#
#    - Residuos de Pearson medios por dia de la semana, para evaluar la
#      presencia de un ciclo semanal de reporte
#    - Cobertura de los intervalos predictivos excluyendo los dias de
#      notificacion anomala
#    - Represas de notificacion: casos observados frente a esperados,
#      sumando cada dia de reporte bajo con los dias de recuperacion
#      inmediatamente posteriores
#    - Figura de la verificacion predictiva marginal, sobre densidades
#
#  Entrada: resultados/fase2_chinagamma/f2cg_M0_ar1_gamma.rds
#           resultados/fase2_chinagamma/def_residuos.csv
#           resultados/fase2_chinagamma/def_dias_fuera_banda95.csv
#  Salida : resultados/fase2_chinagamma/ppc_*.csv y la figura asociada
#
#  No ajusta modelos: opera sobre draws ya guardados.
# =============================================================================

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior); library(ggplot2); library(dplyr); library(here)
})
set_cmdstan_path(Sys.getenv("CMDSTAN", unset = "C:/cmdstan/cmdstan-2.39.0"))

# --- Rutas relativas a la ra\u00edz del repositorio -------------------------------
OUT <- here("resultados", "fase2_chinagamma")
FIG <- here("figuras",    "fase2_chinagamma")
for (p in c(OUT, FIG)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

res   <- read.csv(file.path(OUT, "def_residuos.csv"));        res$fecha <- as.Date(res$fecha)
fuera <- read.csv(file.path(OUT, "def_dias_fuera_banda95.csv")); fuera$fecha <- as.Date(fuera$fecha)
y_obs <- res$obs; fechas <- res$fecha; T_P <- length(y_obs)

fit   <- readRDS(file.path(OUT, "f2cg_M0_ar1_gamma.rds"))
M_rep <- fit$draws("y_rep", format = "matrix")
rm(fit); gc(verbose = FALSE)

# --- PASO 1: residuos por dia de la semana ------------------------------------
cat("=== PASO 1: RESIDUOS POR DIA DE LA SEMANA ===\n")
res$dia <- factor(format(res$fecha, "%u"), levels = as.character(1:7),
                  labels = c("lunes", "martes", "miercoles", "jueves",
                             "viernes", "sabado", "domingo"))
sem <- res %>% group_by(dia) %>%
  summarise(n = n(), resid_medio = mean(resid), resid_sd = sd(resid), .groups = "drop")
print(as.data.frame(sem), digits = 3)
cat(sprintf("  Rango de los residuos medios: [%.2f, %.2f]\n",
            min(sem$resid_medio), max(sem$resid_medio)))
write.csv(sem, file.path(OUT, "ppc_residuos_por_dia_semana.csv"), row.names = FALSE)

# --- PASO 2: cobertura sin los 9 dias de represa ------------------------------
cat("\n=== PASO 2: COBERTURA SIN LOS 9 DIAS DE REPRESA ===\n")
NIV <- c(0.50, 0.80, 0.95)
cobertura <- function(keep) sapply(NIV, function(p) {
  a <- (1 - p) / 2
  lo <- apply(M_rep[, keep, drop = FALSE], 2, quantile, a,     names = FALSE)
  hi <- apply(M_rep[, keep, drop = FALSE], 2, quantile, 1 - a, names = FALSE)
  mean(y_obs[keep] >= lo & y_obs[keep] <= hi)
})
todos  <- cobertura(1:T_P)
sin9   <- cobertura(which(!fechas %in% fuera$fecha))
cob <- data.frame(nivel = paste0("IC", NIV * 100), nominal = NIV * 100,
                  todos_180 = round(100 * todos, 1),
                  sin_9_represas = round(100 * sin9, 1))
print(cob, row.names = FALSE)
write.csv(cob, file.path(OUT, "ppc_cobertura_sin_represas.csv"), row.names = FALSE)

# --- PASO 3: represas, observado vs esperado ----------------------------------
cat("\n=== PASO 3: REPRESAS DE REPORTE ===\n")
bajos <- fuera$fecha[fuera$lado == "ABAJO"]
rep <- do.call(rbind, lapply(bajos, function(d0) {
  i <- which(fechas == d0)
  # ventana: el dia bajo y los dos dias siguientes (recuperacion)
  ix <- i:min(i + 2, T_P)
  data.frame(dia_bajo = d0, ventana = paste(format(fechas[min(ix)]), format(fechas[max(ix)])),
             obs_dia_bajo = y_obs[i], esp_dia_bajo = round(res$f[i]),
             obs_ventana = sum(y_obs[ix]), esp_ventana = round(sum(res$f[ix])),
             razon_obs_esp = round(sum(y_obs[ix]) / sum(res$f[ix]), 2))
}))
print(rep, row.names = FALSE)
cat(sprintf("  En el dia bajo, observado/esperado medio: %.2f | en la ventana de 3 dias: %.2f\n",
            mean(rep$obs_dia_bajo / rep$esp_dia_bajo), mean(rep$razon_obs_esp)))
write.csv(rep, file.path(OUT, "ppc_represas_reporte.csv"), row.names = FALSE)

# --- PASO 4: figura del PPC marginal ------------------------------------------
cat("\n=== PASO 4: FIGURA DEL PPC MARGINAL ===\n")
set.seed(42)
idx <- sample(nrow(M_rep), 200)
dens_rep <- do.call(rbind, lapply(seq_along(idx), function(k) {
  d <- density(M_rep[idx[k], ], from = 0, to = 13000, n = 512)
  data.frame(x = d$x, y = d$y, replica = k)
}))
d_obs <- density(y_obs, from = 0, to = 13000, n = 512)
p <- ggplot() +
  geom_line(data = dens_rep, aes(x, y, group = replica),
            colour = "grey70", alpha = 0.35, linewidth = 0.3) +
  geom_line(data = data.frame(x = d_obs$x, y = d_obs$y), aes(x, y),
            colour = "#1F4E79", linewidth = 1.1) +
  scale_x_continuous(labels = scales::label_comma(big.mark = ".")) +
  labs(x = "Casos diarios", y = "Densidad",
       title = "Verificación predictiva posterior marginal",
       subtitle = "Línea azul: casos observados | Líneas grises: 200 réplicas del predictivo posterior") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40"))
ggsave(file.path(FIG, "fig_ppc_marginal_cg.pdf"), p, width = 9, height = 5, device = cairo_pdf)
ggsave(file.path(FIG, "fig_ppc_marginal_cg.png"), p, width = 9, height = 5, dpi = 200)
cat("  Guardada: fig_ppc_marginal_cg.pdf (y .png de vista rapida)\n")

cat("\n=== VERIFICACION PREDICTIVA COMPLETA ===\n")
