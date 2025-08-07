
## Cargamos librerias 

library(readxl)
library(dplyr)
library(utils)
library(rstan)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(writexl)


# Lectura de datos
DATOS_COVID <- read_excel("C:/Users/geo04/OneDrive/Escritorio/Tesis Maestría/DATOS_COVID_BOGOTA.xlsx")

### Limpieza de Datos

names(DATOS_COVID)[names(DATOS_COVID) == "Fecha de diagnÃ³stico"] <- "Fecha de diagnóstico"

 # Filtrar los datos para excluir los importados
datos_no_importados <- filter(DATOS_COVID, `Tipo de contagio` != "Importado")

# Filtrar los datos para incluir solo los importados
datos_importados <- filter(DATOS_COVID, `Tipo de contagio` == "Importado")



# Agrupar por fecha y calcular el número total de casos
datos_no_importados$Casos <- 1

datos_agrupados_no_importados <- datos_no_importados %>%
  group_by(datos_no_importados$`Fecha de diagnóstico`) %>%
  summarise(Total_Casos = sum(Casos))


# Filtrar los datos para incluir solo los importados
datos_importados$Casos <- 1

# Agrupar por fecha y calcular el número total de casos
datos_agrupados_importados <- datos_importados %>%
  group_by(datos_importados$`Fecha de diagnóstico`) %>%
  summarise(Total_Casos = sum(Casos))


#these libraries need to be loaded

library(utils)
library(rstan)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(writexl)


#write_xlsx(serial_interval,"C:/Users/geo04/Downloads/serial_interval.xlsx")
#write_xlsx(serial_interval,"C:/Users/geo04/Downloads/serial_interval_gamma.xlsx")


d <- datos_agrupados_no_importados[1:180, ]
colnames(d)[colnames(d) == "Total_Casos"] <- "cases"
colnames(d)[colnames(d) == "datos_no_importados$`Fecha de diagnóstico`"] <- "fecha"

serial_interval <- readRDS("C:/Users/geo04/OneDrive/Escritorio/Tesis Maestría/serial-interval.rds")
pad_serial.interval <- data.frame(
  "X" = (length(serial_interval$fit)+1):200,
  "fit" = rep(1e-17, 200)
)

serial_interval = rbind(serial_interval, pad_serial.interval)
plot(serial_interval$X, serial_interval$fit, type="l")


stan_data<-list()
stan_data$N  <- nrow(d)
stan_data$cases <- d$cases
stan_data$N2 <- 60

stan_data$SI <- serial_interval[1:stan_data$N,2]

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# m <- stan_model(file.path('stan-models','bellman.stan'))
m <- stan_model(file.path("C:/Users/geo04/OneDrive/Escritorio/Tesis Maestría/bellman.stan"))

fit <- sampling(m,data=stan_data,iter=2000,warmup=1000,chains=5,thin=5,control = list(adapt_delta = 0.99, max_treedepth = 30))


########################################################################################
######################################################################################
##### Gráficos de Stan ###############################################################33
#####################################################################################

print(fit)  # Muestra los estadísticos de resumen
traceplot(fit)  # Grafica los trazos de las cadenas
stan_diag(fit)  # Diagnóstico de convergencia
traceplot(fit, pars = c("phi", "Rt[12]"))
summary(fit)
summary(fit)$summary[, c("n_eff", "Rhat")]

pairs(fit, pars = c("phi", "weekly_rho", "weekly_sd"))

summary_fit <- summary(fit, pars = c("mu[14]",
                                     "weekly_effect[14]", "weekly_effect[63]", 
                                     "weekly_effect[92]", "weekly_effect[148]", 
                                     "weekly_effect[157]", 
                                     "phi", "weekly_rho", "weekly_rho1"))$summary

summary_fit[, c("n_eff", "Rhat")]


### con ggplot
trace_data <- data.frame(
  iteration = draws_df$.iteration,
  chain = draws_df$.chain,
  value = draws_df[["weekly_effect[1]"]]
)

library(ggplot2)

ggplot(trace_data, aes(x = iteration, y = value, color = factor(chain))) +
  geom_line(alpha = 0.8, linewidth = 0.4) +
  labs(
    title = "Traceplot de weekly_effect[1]",
    x = "Iteración",
    y = "Valor",
    color = "Cadena"
  ) +
  theme_minimal(base_size = 13) +
  scale_color_brewer(palette = "Dark2")


#### Dias seleccionados para la convergencia

library(posterior)
library(dplyr)
library(ggplot2)

# Extraer draws
draws_df <- as_draws_df(fit)

# Días representativos
dias <- c(14, 63, 92, 148, 157)

# Parámetros deseados
params <- c(paste0("weekly_effect[", dias, "]"),
            paste0("mu[", dias, "]"),
            "phi", "weekly_rho", "weekly_rho1")

# Filtra solo los que existen en draws_df
params_validos <- intersect(params, colnames(draws_df))

# Construye el data frame largo manualmente
trace_data <- do.call(rbind, lapply(params_validos, function(param) {
  data.frame(
    iteration = draws_df$.iteration,
    chain = draws_df$.chain,
    parameter = param,
    value = draws_df[[param]]
  )
}))

# Graficar
ggplot(trace_data, aes(x = iteration, y = value, color = factor(chain))) +
  geom_line(alpha = 0.6, linewidth = 0.3) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  labs(
    title = "Traceplots de parámetros seleccionados",
    x = "Iteración",
    y = "Valor",
    color = "Cadena"
  ) +
  theme_minimal(base_size = 12) +
  scale_color_brewer(palette = "Dark2")



## Autocorrelacion
#stan_ac(fit, pars = c("phi"))  # Autocorrelación

library(posterior)
library(bayesplot)
library(dplyr)
library(ggplot2)

# 1. Extraer draws como draws_matrix (formato requerido)
draws_matrix <- as_draws_matrix(fit)

# 2. Días representativos
dias <- c(14, 63, 92, 148, 157)

# 3. Parámetros deseados
params_deseados <- c(paste0("weekly_effect[", dias, "]"),
                     paste0("mu[", dias, "]"),
                     "phi", "weekly_rho", "weekly_rho1")

# 4. Solo los que existen
params_validos <- intersect(params_deseados, colnames(draws_matrix))

# 5. Intentar calcular autocorrelaciones de los válidos
acf_list <- list()

for (param in params_validos) {
  # extrae como matriz con múltiples cadenas
  param_matrix <- draws_matrix[, param, drop = FALSE]
  
  # autocorrelación
  acf_result <- tryCatch({
    bayesplot::mcmc_acf(param_matrix, lags = 30)[[1]]
  }, error = function(e) NULL)
  
  # si acf_result tiene contenido lo almacenamos
  if (!is.null(acf_result) && is.matrix(acf_result)) {
    acf_list[[param]] <- acf_result
  }
}

# 6. Unir a un data.frame largo para ggplot
acf_df <- do.call(rbind, lapply(names(acf_list), function(param) {
  acf_param <- acf_list[[param]]
  do.call(rbind, lapply(1:ncol(acf_param), function(chain) {
    data.frame(
      lag = 0:(nrow(acf_param) - 1),
      acf = acf_param[, chain],
      chain = factor(chain),
      parameter = param
    )
  }))
}))


ggplot(acf_df, aes(x = lag, y = acf, color = chain)) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  labs(
    title = "Autocorrelación de parámetros seleccionados",
    x = "Retardo (lag)",
    y = "Autocorrelación",
    color = "Cadena"
  ) +
  theme_minimal(base_size = 12) +
  scale_color_brewer(palette = "Set1")



## Distribucion posterior
stan_dens(fit, pars = c("phi"))

### Con intervalos de credibilidad
stan_plot(fit, pars = c("phi", "Rt[12]"), ci_level = 0.95)

### Numero efectivo neff y R_t estadistico
stan_rhat(fit)  # Muestra solo los Rhat
stan_ess(fit)   # Número efectivo de muestras

library(posterior)
# Convertir objeto rstan a formato tidy
draws_df <- as_draws_df(fit)

library(dplyr)
library(tidyr)
install.packages("rlang")
library(tidyr)

library(ggplot2)


# Diagnóstico de convergencia con R-hat (debería estar cerca de 1)
print(summary(fit)$summary[, "Rhat"])



### Importamos los resultados de las estimaciones de los hiperparámmetros a un Excel 
## para generar gráficos en Python

library(rstan)
library(writexl)  # También se puede usar openxlsx

# Extraer resumen de los parámetros estimados
fit_summary <- summary(fit)$summary  

# Convertir en data.frame para manipulación más sencilla
df_results <- as.data.frame(fit_summary)

# Agregar nombres de parámetros como columna
df_results$Parameter <- rownames(fit_summary)
rownames(df_results) <- NULL  # Eliminar nombres de filas


# Guardar en Excel usando writexl
write_xlsx(df_results, "resultados_stan.xlsx")
getwd()



