# ==============================================================================
#  FIGURA EXPLORATORIA: JUSTIFICACIÓN RESOLUCIÓN DIARIA
#  Autor: George Harrison Plata Bello — Tesis Maestría Estadística UNAL
# ==============================================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(lubridate)
library(gridExtra)
library(grid)
library(scales)
library(here)

# ==============================================================================
# 1. CARGAR Y PREPARAR DATOS
# ==============================================================================

ruta <- here("MODELO_STAN_4", "datos_agregados.xlsx")

# Hoja NO_IMPORTADOS: fechas vienen como número serial Excel en character ("43904")
df_raw <- as.data.frame(read_excel(ruta, sheet = "NO_IMPORTADOS"))
names(df_raw) <- c("fecha", "casos")
df_raw$fecha <- as.Date(as.numeric(df_raw$fecha), origin = "1899-12-30")
df_raw$casos <- as.integer(df_raw$casos)

# Verificación
cat("Primeras filas cargadas:\n")
print(head(df_raw, 5))
cat(sprintf("\nRango fechas: %s a %s\n", min(df_raw$fecha), max(df_raw$fecha)))
cat(sprintf("Total registros: %d | Total casos: %s\n\n",
            nrow(df_raw), format(sum(df_raw$casos), big.mark=",")))

# Período de análisis: 14-marzo a 09-septiembre 2020 (180 días)
FECHA_INICIO <- as.Date("2020-03-14")
FECHA_FIN    <- as.Date("2020-09-09")

grilla <- data.frame(fecha = seq(FECHA_INICIO, FECHA_FIN, by = "day"))
df <- left_join(grilla, df_raw, by = "fecha")
df$casos[is.na(df$casos)] <- 0L
df$dia_sem <- wday(df$fecha, label = TRUE, abbr = TRUE, week_start = 1)

cat(sprintf("Período análisis: %s -> %s (%d días)\n", FECHA_INICIO, FECHA_FIN, nrow(df)))
cat(sprintf("Total casos en período: %s | Media: %.1f/día\n\n",
            format(sum(df$casos), big.mark=","), mean(df$casos)))

# ==============================================================================
# 2. FECHAS DE INTERVENCIÓN
# ==============================================================================

intervenciones <- data.frame(
  fecha    = as.Date(c("2020-03-25","2020-05-04","2020-06-19",
                       "2020-07-03","2020-07-15","2020-08-01")),
  etiqueta = c("Cuarentena\nnacional\n(Dto. 457)",
               "Apertura\nsectorial\n(CIIU)",
               "Dia sin\nIVA","Dia sin\nIVA",
               "Cuarentenas\nlocalidades\n(Dto. 169)",
               "Cuarentenas\nlocalidades\n(Dto. 173)"),
  color    = c("#D62728","#2CA02C","#FF7F0E","#FF7F0E","#9467BD","#9467BD"),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 3. PANEL 1: SERIE DE TIEMPO
# ==============================================================================

df$mm7 <- as.numeric(stats::filter(df$casos, rep(1/7,7), sides=2))

p1 <- ggplot(df, aes(x=fecha, y=casos)) +
  geom_col(fill="#AEC6CF", alpha=0.7, width=1) +
  geom_line(aes(y=mm7), color="#1F4E79", linewidth=0.8, na.rm=TRUE) +
  geom_vline(data=intervenciones, aes(xintercept=fecha, color=color),
             linetype="dashed", linewidth=0.7, show.legend=FALSE) +
  geom_text(data=intervenciones,
            aes(x=fecha, y=max(df$casos, na.rm=TRUE)*0.92,
                label=etiqueta, color=color),
            angle=90, hjust=1, vjust=-0.3, size=2.2,
            lineheight=0.85, show.legend=FALSE) +
  scale_color_identity() +
  scale_x_date(date_breaks="1 month", date_labels="%b\n%Y") +
  scale_y_continuous(labels=comma) +
  labs(title="(a) Serie de tiempo diaria con eventos de intervencion",
       subtitle="Barras: casos diarios  |  Linea: media movil 7 dias",
       x=NULL, y="Casos diarios") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(size=9, face="bold"),
        plot.subtitle=element_text(size=8, color="gray40"))

# ==============================================================================
# 4. PANEL 2: BOXPLOT DÍA DE SEMANA
# ==============================================================================

kw   <- kruskal.test(casos ~ dia_sem, data=df)
kw_p <- ifelse(kw$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", kw$p.value))

p2 <- ggplot(df, aes(x=dia_sem, y=casos, fill=dia_sem)) +
  geom_boxplot(outlier.size=0.8, outlier.alpha=0.5, alpha=0.75) +
  scale_fill_brewer(palette="Set2") +
  labs(title="(b) Distribucion de casos por dia de la semana",
       subtitle=paste0("Kruskal-Wallis: ", kw_p),
       x="Dia de la semana", y="Casos diarios") +
  theme_bw(base_size=10) +
  theme(legend.position="none",
        plot.title=element_text(size=9, face="bold"),
        plot.subtitle=element_text(size=8, color="gray40"))

# ==============================================================================
# 5. PANEL 3: ACF
# ==============================================================================

acf_data <- acf(df$casos, lag.max=21, plot=FALSE)
ic_val   <- qnorm(0.975) / sqrt(nrow(df))

acf_df <- data.frame(
  lag = as.numeric(acf_data$lag[-1]),
  acf = as.numeric(acf_data$acf[-1])
)

p3 <- ggplot(acf_df, aes(x=lag, y=acf)) +
  geom_hline(yintercept=0, color="black", linewidth=0.4) +
  geom_hline(yintercept= ic_val, linetype="dashed", color="#D62728", linewidth=0.5) +
  geom_hline(yintercept=-ic_val, linetype="dashed", color="#D62728", linewidth=0.5) +
  geom_segment(aes(x=lag, xend=lag, y=0, yend=acf), color="#1F4E79", linewidth=0.8) +
  geom_point(aes(color=abs(acf) > ic_val), size=2) +
  scale_color_manual(values=c("FALSE"="gray60","TRUE"="#D62728"), guide="none") +
  geom_vline(xintercept=c(7,14), linetype="dotted", color="orange", linewidth=0.6) +
  scale_x_continuous(breaks=seq(0,21,by=3)) +
  labs(title="(c) Funcion de autocorrelacion (ACF) - casos diarios",
       subtitle="Lineas rojas: IC 95%  |  Lineas naranjas: lag 7 y 14 dias",
       x="Lag (dias)", y="Autocorrelacion") +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(size=9, face="bold"),
        plot.subtitle=element_text(size=8, color="gray40"))

# ==============================================================================
# 6. PANEL 4: Rt DIARIO VS SEMANAL
# ==============================================================================

g_si <- numeric(30)
g_si[1] <- pgamma(1.5, shape=6.5, rate=0.62)
for (s in 2:30) {
  g_si[s] <- pgamma(s+0.5, shape=6.5, rate=0.62) - pgamma(s-0.5, shape=6.5, rate=0.62)
}
g_si <- g_si / sum(g_si)

y_vec  <- as.numeric(df$casos)
N      <- length(y_vec)
Rt_emp <- rep(NA_real_, N)

for (t in 8:N) {
  denom <- sum(sapply(1:min(t-1,30), function(s) y_vec[t-s] * g_si[s]))
  if (is.finite(denom) && denom > 5) Rt_emp[t] <- y_vec[t] / denom
}

Rt_sem <- rep(NA_real_, N)
for (i in seq(1, N-6, by=7)) {
  idx <- i:min(i+6, N)
  val <- mean(Rt_emp[idx], na.rm=TRUE)
  if (is.finite(val)) Rt_sem[idx] <- val
}

df_p4 <- data.frame(
  fecha      = df$fecha,
  Rt_diario  = as.numeric(Rt_emp),
  Rt_semanal = as.numeric(Rt_sem)
)

p4 <- ggplot(df_p4, aes(x=fecha)) +
  geom_line(aes(y=Rt_diario,  color="Diario"),            linewidth=0.5, alpha=0.7, na.rm=TRUE) +
  geom_step(aes(y=Rt_semanal, color="Semanal (promedio)"), linewidth=1.0,            na.rm=TRUE) +
  geom_hline(yintercept=1, linetype="dashed", color="black", linewidth=0.5) +
  geom_vline(data=intervenciones, aes(xintercept=fecha),
             linetype="dotted", color="gray60", linewidth=0.5) +
  scale_color_manual(values=c("Diario"="#1F4E79","Semanal (promedio)"="#D62728"),
                     name="Resolucion") +
  scale_x_date(date_breaks="1 month", date_labels="%b") +
  scale_y_continuous(limits=c(0, NA)) +
  labs(title=expression(paste("(d) ", R[t], " empirico: resolucion diaria vs semanal")),
       subtitle="IS: Bi et al. (2020)  |  Linea punteada: Rt = 1",
       x=NULL, y=expression(R[t])) +
  theme_bw(base_size=10) +
  theme(legend.position="bottom",
        legend.key.size=unit(0.4,"cm"),
        plot.title=element_text(size=9, face="bold"),
        plot.subtitle=element_text(size=8, color="gray40"))

# ==============================================================================
# 7. COMBINAR Y GUARDAR
# ==============================================================================

titulo_general <- textGrob(
  "Justificacion empirica de resolucion temporal diaria\nCOVID-19 Bogota, 14-marzo a 09-septiembre 2020",
  gp = gpar(fontsize=12, fontface="bold")
)

figura <- arrangeGrob(p1, p2, p3, p4, ncol=2, nrow=2, top=titulo_general)

ggsave(
  filename = here("MODELO_STAN_4", "fig_exploratorio_resolucion_diaria.pdf"),
  plot=figura, width=14, height=10, units="in", dpi=300
)
ggsave(
  filename = here("MODELO_STAN_4", "fig_exploratorio_resolucion_diaria.png"),
  plot=figura, width=14, height=10, units="in", dpi=300
)

cat("\n✓ Figura guardada:\n")
cat("  fig_exploratorio_resolucion_diaria.pdf\n")
cat("  fig_exploratorio_resolucion_diaria.png\n\n")

# ==============================================================================
# 8. ESTADÍSTICAS PARA REDACCIÓN
# ==============================================================================

cat("=== ESTADÍSTICAS PARA REDACCIÓN ===\n\n")
cat(sprintf("ACF lag 1: %.3f\n", acf_data$acf[2]))
cat(sprintf("ACF lag 2: %.3f\n", acf_data$acf[3]))
cat(sprintf("ACF lag 7: %.3f\n", acf_data$acf[8]))
cat(sprintf("IC 95%% ACF: +/- %.3f\n\n", ic_val))

media_dia <- df %>% group_by(dia_sem) %>%
  summarise(media=round(mean(casos),1), .groups="drop") %>%
  arrange(desc(media))
cat("Media de casos por dia de semana:\n")
print(as.data.frame(media_dia))

cat(sprintf("\nKruskal-Wallis: chi2=%.2f, df=%d, %s\n\n",
            kw$statistic, kw$parameter, kw_p))

cat("Rt en fechas de intervencion (diario vs semanal):\n")
for (f in as.character(intervenciones$fecha)) {
  idx <- which(df_p4$fecha == as.Date(f))
  if (length(idx) > 0) {
    rd <- df_p4$Rt_diario[idx];  rs <- df_p4$Rt_semanal[idx]
    cat(sprintf("  %s: Rt_diario=%s | Rt_semanal=%s\n", f,
                ifelse(is.na(rd),"NA",sprintf("%.2f",rd)),
                ifelse(is.na(rs),"NA",sprintf("%.2f",rs))))
  }
}
