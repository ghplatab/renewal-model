library(readr)
library(dplyr)
library(readr)
library(ggplot2)

Casos_COVID_19 <- read_csv("C:/Users/geo04/Downloads/Casos_positivos_de_COVID-19_en_Colombia._20250710.csv")


Casos_COVID_19$`Fecha de diagnóstico` <- as.Date(Casos_COVID_19$`Fecha de diagnóstico`)


### Tipo de contagio

library(ggplot2)
library(dplyr)
library(scales)

# Contar casos por tipo de contagio
conteo_tipos <- Casos_COVID_19 %>%
  count(`Tipo de contagio`) %>%
  arrange(desc(n))

# Gráfico limpio y elegante sin etiquetas
ggplot(conteo_tipos, aes(x = reorder(`Tipo de contagio`, -n), y = n)) +
  geom_bar(stat = "identity", fill = "#76AADB", width = 0.6) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Distribución casos COVID-19 por tipo de contagio",
    x = "Tipo de contagio",
    y = "Número de casos"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

## Rngo de tiempo
range(Casos_COVID_19$`Fecha de diagnóstico`, na.rm = TRUE)


## Filtramos por tipo de contagio

library(ggplot2)
library(dplyr)
library(scales)

# Filtrar y acumular casos
casos_filtrados <- Casos_COVID_19[Casos_COVID_19$`Tipo de contagio` != "Importado", ]

acumulado_por_fecha <- casos_filtrados %>%
  group_by(`Fecha de diagnóstico`) %>%
  summarise(casos_diarios = n()) %>%
  arrange(`Fecha de diagnóstico`) %>%
  mutate(casos_acumulados = cumsum(casos_diarios))

# Gráfico estilizado
ggplot(acumulado_por_fecha, aes(x = `Fecha de diagnóstico`, y = casos_acumulados)) +
  geom_line(color = "#2C3E50", size = 1.2) +
  labs(
    title = "Casos acumulados de COVID-19",
    x = "Fecha de diagnóstico",
    y = "Casos acumulados"
  ) +
  scale_y_continuous(labels = comma) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )



### Lectura de los shape files

library(sf)

# Leer shapefile desde zip
departamentos_sf <- st_read("C:/Users/geo04/Downloads/departamentos")
municipos_sf <- st_read("C:/Users/geo04/Downloads/municipios")

# Verificar el resultado
plot(departamentos_sf)
class(departamentos_sf)


### Conteo de casos por departamento

library(dplyr)

casos_por_departamento <- Casos_COVID_19 %>%
  filter(`Tipo de contagio` != "Importado") %>%
  group_by(`Nombre departamento`) %>%  # Asegúrarse que esta columna tenga el nombre correcto
  summarise(casos = n())


names(departamentos_sf)



# Agrupar casos no importados por departamento
casos_por_departamento <- Casos_COVID_19 %>%
  filter(`Tipo de contagio` != "Importado") %>%
  group_by(`Nombre departamento`) %>%
  summarise(casos = n())

# Unir con shapefile (ajusta NOMBRE_DPT según corresponda)
departamentos_mapa <- left_join(departamentos_sf, 
                                casos_por_departamento, 
                                by = c("NOMBRE_DPT" = "Nombre departamento"))


### Gráfico del mapa


library(ggplot2)

ggplot(departamentos_mapa) +
  geom_sf(aes(fill = casos), color = "gray60", size = 0.2) +
  scale_fill_gradientn(
    colours = c("#FEE5D9", "#FCAe91", "#FB6A4A", "#DE2D26", "#A50F15"),
    na.value = "gray90",
    name = "Casos"
  ) +
  labs(
    title = "Distribución de casos COVID-19 (sin importados)",
    subtitle = "Acumulado por departamento",
    caption = "Fuente: Datos COVID y shapefile oficial"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.background = element_rect(fill = "white")
  )


#### con la librería leafleft

#install.packages("leaflet")

library(leaflet)

# Crear paleta de colores
paleta <- colorQuantile("Reds", departamentos_mapa$casos, n = 5, na.color = "gray90")

# Crear mapa
leaflet(departamentos_mapa) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(
    fillColor = ~paleta(casos),
    color = "white",
    weight = 1,
    opacity = 1,
    fillOpacity = 0.7,
    smoothFactor = 0.5,
    highlightOptions = highlightOptions(color = "black", weight = 2, bringToFront = TRUE),
    label = ~paste0(NOMBRE_DPT, ": ", casos, " casos"),
    popup = ~paste("<strong>", NOMBRE_DPT, "</strong><br/>Casos: ", casos)
  ) %>%
  addLegend(
    pal = paleta,
    values = ~casos,
    title = "Casos COVID-19",
    position = "bottomright"
  )



##### Con la libreria ggspatial


library(ggplot2)
library(ggspatial)
library(sf)

# Convertimos a categorías para colores discretos
departamentos_mapa$casos_cat <- cut(
  departamentos_mapa$casos,
  breaks = c(0, 1, 10, 50, 1000),
  labels = c("1", "2-10", "11-50", "51+"),
  include.lowest = TRUE
)

ggplot(departamentos_mapa) +
  geom_sf(aes(fill = casos_cat), color = "gray70", size = 0.2) +
  scale_fill_manual(
    values = c("1" = "#a1d99b", "2-10" = "#74c476", "11-50" = "#31a354", "51+" = "#006d2c"),
    na.value = "gray90",
    name = "Casos"
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf() +
  labs(
    title = "Distribución de casos COVID-19 por departamento",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal()



###### Codigo completo con valores reales

# Cargar librerías necesarias
library(ggplot2)
library(ggspatial)
library(sf)
library(dplyr)
library(stringr)

# --- 1. Normalizar nombres de departamentos para unir datos correctamente ---
# En el shapefile (departamentos_sf) y la tabla de casos (casos_por_departamento)
casos_por_departamento <- casos_por_departamento %>%
  mutate(nombre_dpto = str_to_upper(str_trim(`Nombre departamento`)))

departamentos_sf <- departamentos_sf %>%
  mutate(nombre_dpto = str_to_upper(str_trim(NOMBRE_DPT)))

# --- 2. Unir los datos de casos al shapefile ---
departamentos_mapa <- left_join(departamentos_sf, casos_por_departamento, by = "nombre_dpto")

# --- 3. Rellenar NA con 0 (opcional, para visualización completa) ---
departamentos_mapa$casos[is.na(departamentos_mapa$casos)] <- 0

# --- 4. Crear categorías de casos ajustadas a la distribución real ---
departamentos_mapa$casos_cat <- cut(
  departamentos_mapa$casos,
  breaks = c(0, 5000, 20000, 50000, 100000, 400000),
  labels = c("≤5k", "5k–20k", "20k–50k", "50k–100k", "100k+"),
  include.lowest = TRUE
)

# --- 5. Definir colores estilo epidemiólogico (rojizos) ---
colores <- c("≤5k" = "#fee5d9", 
             "5k–20k" = "#fcae91", 
             "20k–50k" = "#fb6a4a", 
             "50k–100k" = "#de2d26", 
             "100k+" = "#a50f15")

# --- 6. Generar el mapa final ---
ggplot(departamentos_mapa) +
  geom_sf(aes(fill = casos_cat), color = "gray70", size = 0.2) +
  scale_fill_manual(
    values = colores,
    na.value = "gray90",
    name = "Casos COVID-19"
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf() +
  labs(
    title = "Distribución de casos COVID-19 por departamento",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white", color = NA)
  )


##############################################################################3
################################################################################3
#### MISMO PIPELINE PERO PARA MUNICIPIOS


### Lectura de los shape files

library(sf)

# Leer shapefile desde zip
#departamentos_sf <- st_read("C:/Users/geo04/Downloads/departamentos")
municipos_sf <- st_read("C:/Users/geo04/Downloads/municipios")

# Verifica el resultado
plot(municipos_sf)
class(municipos_sf)


library(dplyr)

casos_por_municipio <- Casos_COVID_19 %>%
  filter(`Tipo de contagio` != "Importado") %>%
  group_by(`Nombre municipio`) %>%
  summarise(casos = n())


## Normializar y unir nombres

library(stringr)

casos_por_municipio <- casos_por_municipio %>%
  mutate(nombre_mpio = str_to_upper(str_trim(`Nombre municipio`)))

municipios_sf <- municipos_sf %>%
  mutate(nombre_mpio = str_to_upper(str_trim(MPIO_CNMBR)))

municipios_mapa <- left_join(municipios_sf, casos_por_municipio, by = "nombre_mpio")
municipios_mapa$casos[is.na(municipios_mapa$casos)] <- 0



#### breaks

municipios_mapa$casos_cat <- cut(
  municipios_mapa$casos,
  breaks = c(0, 50, 200, 1000, 5000, 20000),
  labels = c("≤50", "51–200", "201–1k", "1k–5k", "5k–20k"),
  include.lowest = TRUE
)


colores_mpios <- c(
  "≤50" = "#fee5d9",
  "51–200" = "#fcae91",
  "201–1k" = "#fb6a4a",
  "1k–5k" = "#de2d26",
  "5k–20k" = "#a50f15"
)


#### Genreacion del Mapa

library(ggplot2)
library(ggspatial)

ggplot(municipios_mapa) +
  geom_sf(aes(fill = casos_cat), color = "gray70", size = 0.1) +
  scale_fill_manual(
    values = colores_mpios,
    na.value = "gray90",
    name = "Casos COVID-19"
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf() +
  labs(
    title = "Distribución de casos COVID-19 por municipio",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray95"),
    panel.background = element_rect(fill = "white", color = NA)
  )


