library(dplyr)
library(ggplot2)
library(here)

# ── Rutas ────────────────────────────────────────────────────────────────────
RUTA <- here("MODELO_STAN_DEFINITIVO")

# ── Cargar ───────────────────────────────────────────────────────────────────
draws <- readRDS(file.path(RUTA, "draws_AR1_Gamma_definitivo.rds"))

# ── Datos observados ─────────────────────────────────────────────────────────
importados_obs <- data.frame(
  fecha = seq(as.Date("2020-03-14"), by = "day", length.out = 12),
  casos = c(7, 8, 2, 5, 4, 6, 18, 4, 1, 17, 23, 5)
)

# ── Diagnóstico: cuánto pesa el exógeno en f(t) ──────────────────────────────
f_mat <- draws$f  # 6000 × 180

f_summary_exo <- data.frame(
  fecha = seq(as.Date("2020-03-14"), by = "day", length.out = 12),
  media = apply(f_mat[, 1:12], 2, mean),
  q025  = apply(f_mat[, 1:12], 2, quantile, 0.025),
  q975  = apply(f_mat[, 1:12], 2, quantile, 0.975)
)

cat("=== f(t) durante período exógeno (días 1-12) ===\n")
print(f_summary_exo)

cat("\nMedia f(t) días  1-12:", round(mean(apply(f_mat[,  1:12], 2, mean)), 2), "\n")
cat("Media f(t) días 13-30:", round(mean(apply(f_mat[, 13:30], 2, mean)), 2), "\n")
cat("Media f(t) días 31-60:", round(mean(apply(f_mat[, 31:60], 2, mean)), 2), "\n")





############ Gráfica Exogeno

library(dplyr)
library(ggplot2)
library(patchwork)

col_azul <- "#2166ac"
col_rojo <- "#d6604d"
col_gris <- "grey50"

tema_tesis <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    axis.title       = element_text(size = 10)
  )

RUTA <- here::here("MODELO_STAN_DEFINITIVO")

draws <- readRDS(file.path(RUTA, "draws_AR1_Gamma_definitivo.rds"))
f_mat <- draws$f  # 6000 × 180

importados_obs <- data.frame(
  fecha = seq(as.Date("2020-03-14"), by = "day", length.out = 12),
  casos = c(7, 8, 2, 5, 4, 6, 18, 4, 1, 17, 23, 5)
)

# ── f(t) completo — 180 días ──────────────────────────────────────────────────
fechas_total <- seq(as.Date("2020-03-14"), by = "day", length.out = 180)

f_total <- data.frame(
  fecha = fechas_total,
  media = apply(f_mat, 2, mean),
  q025  = apply(f_mat, 2, quantile, 0.025),
  q975  = apply(f_mat, 2, quantile, 0.975),
  q25   = apply(f_mat, 2, quantile, 0.25),
  q75   = apply(f_mat, 2, quantile, 0.75)
)

# ── Fecha fin período exógeno ─────────────────────────────────────────────────
fecha_fin_exo <- as.Date("2020-03-25")

# ── PANEL A — f(t) completo con período exógeno sombreado ────────────────────
pA <- ggplot(f_total, aes(x = fecha)) +
  # Sombrear período exógeno
  annotate("rect",
           xmin = as.Date("2020-03-14"),
           xmax = fecha_fin_exo,
           ymin = -Inf, ymax = Inf,
           fill = col_azul, alpha = 0.07) +
  annotate("text",
           x = as.Date("2020-03-14") + 3, y = max(f_total$q975) * 0.90,
           label = "Período\nexógeno\n(μ(t) ≈ f(t))",
           hjust = 0, size = 3, color = col_azul) +
  # Bandas IC
  geom_ribbon(aes(ymin = q025, ymax = q975),
              fill = col_azul, alpha = 0.18) +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              fill = col_azul, alpha = 0.35) +
  geom_line(aes(y = media), color = col_azul, linewidth = 0.9) +
  # Casos importados observados
  geom_point(data = importados_obs,
             aes(x = fecha, y = casos),
             color = col_rojo, size = 2.5, shape = 19) +
  # Línea vertical fin exógeno
  geom_vline(xintercept = fecha_fin_exo,
             linetype = "dashed", color = col_azul,
             linewidth = 0.6) +
  scale_x_date(date_labels = "%d %b\n%Y", date_breaks = "3 weeks") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = expression("A — Ajuste " * f(t) *
                            ": componente exógeno como semillero epidémico"),
    subtitle = paste0(
      "Zona azul: período exógeno (14-25 mar)  |  ",
      "Puntos rojos: importados observados  |  ",
      "Banda: IC50% y IC95%"
    ),
    x = NULL,
    y = expression("Casos esperados " * f(t))
  ) +
  tema_tesis

# ── PANEL B — Zoom período exógeno: f(t) vs importados ───────────────────────
f_exo <- filter(f_total, fecha <= fecha_fin_exo)

pB <- ggplot(f_exo, aes(x = fecha)) +
  geom_ribbon(aes(ymin = q025, ymax = q975),
              fill = col_azul, alpha = 0.18) +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              fill = col_azul, alpha = 0.35) +
  geom_line(aes(y = media), color = col_azul, linewidth = 0.9) +
  geom_point(data = importados_obs,
             aes(x = fecha, y = casos),
             color = col_rojo, size = 3, shape = 19) +
  geom_segment(
    data = data.frame(importados_obs, media_f = f_exo$media),
    aes(x = fecha, xend = fecha, y = casos, yend = media_f),
    color = col_rojo, linewidth = 0.4, linetype = "dotted"
  ) +
  scale_x_date(date_labels = "%d %b", date_breaks = "3 days") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = expression("B — Zoom: " * f(t) %~~% mu(t) *
                            " durante el período exógeno"),
    subtitle = paste0(
      "Media f(t) = 9.74  |  ",
      "Puntos rojos: importados observados  |  ",
      "Línea punteada: residuo"
    ),
    x = NULL,
    y = expression(f(t) %~~% mu(t))
  ) +
  tema_tesis

# ── Combinar ─────────────────────────────────────────────────────────────────
fig5 <- pA / pB +
  plot_annotation(
    title = expression(
      "Componente exógeno " * mu(t) *
        " como condición inicial del proceso de renovación"
    ),
    subtitle = paste0(
      "Durante 14-25 mar 2020, f(t) ≈ μ(t) pues los infectados previos son mínimos. ",
      "A partir del 26 mar los casos generados se vuelven endógenos y alimentan R_t."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  )

# ── Guardar ───────────────────────────────────────────────────────────────────
ggsave(file.path(RUTA, "fig5_exogeno.pdf"), fig5,
       width = 10, height = 9, device = cairo_pdf)
ggsave(file.path(RUTA, "fig5_exogeno.png"), fig5,
       width = 10, height = 9, dpi = 300)

message("✓ fig5_exogeno.pdf y fig5_exogeno.png generados en ", RUTA)


