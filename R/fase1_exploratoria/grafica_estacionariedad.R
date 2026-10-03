library(ggplot2)
library(here)

# Función denominador varianza
calcular_D <- function(rho1, rho2) {
  ifelse(rho2 >= 1, NA, 
         1 - rho1^2 - rho2^2 - 2*rho1^2*rho2/(1-rho2))
}

# Grilla para región de estacionariedad
rho1_seq <- seq(0, 1, length.out = 200)
rho2_seq <- seq(0, 1, length.out = 200)
grid <- expand.grid(rho1 = rho1_seq, rho2 = rho2_seq)

# Evaluar condiciones
grid$estacionario <- with(grid, 
                          (rho1 + rho2 < 1) & (rho2 - rho1 < 1) & (abs(rho2) < 1)
)

grid$D <- with(grid, calcular_D(rho1, rho2))

# Cargar draws posteriores (modelo AR2+Gamma referencia del análisis preliminar)
fit_ref <- readRDS(here("MODELO_STAN_4", "EXOGENO_GAMMA_ANALISIS_PRELIMINAR",
                        "fit_gamma_referencia_NB2.rds"))
draws <- fit_ref$draws(c("rho1", "rho2"), format = "df")
draws_sample <- draws[sample(nrow(draws), 1000), ]

# Gráfico
p <- ggplot() +
  # Región de estacionariedad
  geom_raster(data = grid, 
              aes(x = rho1, y = rho2, fill = estacionario),
              alpha = 0.3) +
  scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "lightblue"),
                    guide = "none") +
  # Contorno D = 0
  geom_contour(data = grid,
               aes(x = rho1, y = rho2, z = D),
               breaks = 0, color = "red", linetype = "dashed", 
               linewidth = 1) +
  # Contorno D = 0.1
  geom_contour(data = grid,
               aes(x = rho1, y = rho2, z = D),
               breaks = 0.1, color = "orange", linetype = "dotted",
               linewidth = 0.8) +
  # Draws posteriores
  geom_point(data = draws_sample,
             aes(x = rho1, y = rho2),
             alpha = 0.3, size = 1, color = "darkblue") +
  # Frontera estacionariedad
  geom_abline(slope = -1, intercept = 1, 
              linetype = "solid", color = "black", linewidth = 0.5) +
  annotate("text", x = 0.3, y = 0.65, 
           label = expression(rho[1] + rho[2] == 1),
           size = 4, angle = -45) +
  annotate("text", x = 0.7, y = 0.05,
           label = "D = 0 (varianza indefinida)",
           size = 3.5, color = "red") +
  labs(
    x = expression(rho[1]),
    y = expression(rho[2]),
    title = "Región de estacionariedad AR(2) y masa posterior"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("Figuras/region_estacionariedad_ar2.pdf", p, 
       width = 8, height = 7)



getwd()
