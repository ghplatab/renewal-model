# Modelo Bayesiano de Renovación COVID-19 — Bogotá

**Tesis de Maestría en Estadística · Universidad Nacional de Colombia**

Implementación de un modelo de ecuación de renovación derivado de un proceso de
ramificación dependiente de la edad (Bellman-Harris), para estimar el número
reproductivo efectivo R(t) de COVID-19 en Bogotá durante los primeros 180 días
del brote (14 de marzo – 9 de septiembre de 2020).

**Modelo definitivo adoptado:** AR(1) sobre log R(t) + componente exógeno Gamma
+ verosimilitud Binomial Negativa tipo 2 (NB-2) + intervalo serial China-Gamma
de Bi et al. (2020).

---

## Estructura del repositorio

```
renewal-model/
├── data/
│   └── datos_agregados.xlsx        # hojas: NO_IMPORTADOS, IMPORTADOS
│
├── stan/                           # modelos Stan
│   ├── renewal_bogota_nb2_ar1_gamma.stan        # MODELO DEFINITIVO
│   ├── renewal_bogota_nb2_ar2_gamma.stan        # AR(2), sensibilidad
│   ├── renewal_bogota_nb2_rw_gamma.stan         # paseo aleatorio, sensibilidad
│   ├── renewal_bogota_nb2_ar1_exponencial.stan  # exógeno exponencial, sensibilidad
│   └── renewal_bogota_poisson_ar2_gamma.stan    # Poisson, fase exploratoria
│
├── R/
│   ├── fase1_exploratoria/         # filtro de convergencia: 10 especificaciones
│   │   ├── analisis_bogota_ar2_gamma.R
│   │   ├── metricas_loo_sensibilidad_AR.R
│   │   └── grafica_estacionariedad.R
│   │
│   └── fase2_definitiva/           # modelo definitivo y diagnósticos
│       ├── 01_cobertura_predictiva.R
│       ├── 02_diagnosticos_posteriores.R
│       ├── 03_corridas_ar1_y_rw.R
│       ├── 04_trayectorias_ar1_vs_rw.R
│       ├── 05_consistencia_corrida_final.R
│       ├── 06_fase2_china_gamma.R
│       ├── 07_cifras_modelo_definitivo.R
│       ├── 08_figuras_modelo_definitivo.R
│       ├── 09_respaldo_intervenciones.R
│       └── 10_verificacion_predictiva.R
│
├── resultados/                     # salidas numéricas (CSV)
└── README.md
```

> **Sobre los archivos grandes:** los ajustes MCMC (`.rds`, ~300 MB en total) y
> las figuras generadas no se versionan; están en `.gitignore` porque se
> regeneran ejecutando los scripts. Los resultados numéricos sí se incluyen,
> en `resultados/`, para que cada cifra del documento sea rastreable sin
> necesidad de volver a correr el modelo.

---

## Datos

**Archivo:** `data/datos_agregados.xlsx`

| Hoja | Contenido |
|------|-----------|
| `NO_IMPORTADOS` | Casos locales diarios (entran en la verosimilitud) |
| `IMPORTADOS` | Casos importados diarios (componente exógeno, días 1–12) |

**Período:** 14 de marzo – 9 de septiembre de 2020 (180 días, 235 030 casos).
**Fuente:** portal de datos abiertos del Estado colombiano; Sistema de Vigilancia
en Salud Pública del Distrito Capital de Bogotá.

---

## Modelo estadístico

```
y(t) ~ NB-2(f(t), φ)

f(t) = μ(t) + R_t · Σ_{s=1}^{t-1} f(s) · g(t-s)

log R_t = μ_R + ε_t,     ε_t = ρ₁ ε_{t-1} + η_t,    η_t ~ N(0, σ_η²)

μ(t) ~ Gamma(1.401, 0.168)    [exógeno, activo solo en los días 1–12]
```

`g(·)` es el intervalo serial discretizado (longitud máxima 30 días) y `φ` el
parámetro de sobredispersión de la NB-2. Nótese que la convolución se construye
sobre `f(s)` —la media del proceso— y no sobre los conteos observados `y(s)`:
es la forma recursiva de la ecuación de renovación.

El componente exógeno se calibró por **método de momentos** sobre los 12 días
previos a la cuarentena (14–25 de marzo), con media empírica 8,33 casos/día y
varianza 49,52.

### Distribuciones previas

| Parámetro | Previa | Interpretación |
|-----------|--------|----------------|
| `μ_R` | N(0, 0.3) | Nivel de log R(t); centrado en R(t) ≈ 1 |
| `ρ₁` | N(0.7, 0.15) truncada en [0,1] | Persistencia diaria |
| `σ_η` | N⁺(0, 0.2) | Escala de las innovaciones diarias |
| `φ` | N⁺(0, 5) | Sobredispersión NB-2 |

En todos los casos el segundo argumento es una **desviación estándar**, no una
varianza, siguiendo la convención de Stan.

---

## Resultados del modelo definitivo

**Configuración del muestreo:** 4 cadenas × (1500 calentamiento + 1500 muestreo)
= 6000 muestras posteriores; `adapt_delta = 0.99`, `max_treedepth = 12`,
`seed = 42`.

### Posteriores

| Parámetro | Media | IC 95 % |
|-----------|-------|---------|
| μ_R | 0,152 | [−0,088, 0,324] |
| ρ₁ | 0,893 | [0,689, 0,995] |
| σ_η | 0,092 | [0,030, 0,192] |
| φ | 3,249 | [2,586, 3,992] |

> **Advertencia sobre μ_R.** Se conserva como parámetro de la especificación,
> pero **no se interpreta como nivel de equilibrio de la transmisión**: la
> persistencia estimada se sitúa próxima a la raíz unitaria, y un modelo de
> paseo aleatorio —en el que μ_R ni siquiera está definido— resulta
> predictivamente indistinguible.

### Diagnósticos MCMC

| Métrica | Valor |
|---------|-------|
| R̂ máximo | 1,0021 |
| ESS mínimo (bulk / tail) | 1923 / 2236 |
| Divergencias | 48 / 6000 (0,8 %) |
| Treedepth saturadas | 0 |

### Trayectoria de R(t)

| | |
|---|---|
| Máximo | **1,70** el 1 de abril de 2020 |
| Primer día por debajo de 1 | **17 de agosto de 2020** |
| Cobertura empírica de los intervalos predictivos | 70,0 % (IC 50) · 87,2 % (IC 80) · 95,0 % (IC 95) |

Los intervalos centrales resultan conservadores: las bandas predictivas son más
anchas de lo que los datos requieren en la región central de la distribución.

---

## Comparación de especificaciones del proceso latente

| Modelo | elpd_loo | SE | Δ vs AR(1) | SE(Δ) |
|--------|----------|-----|-----------|-------|
| AR(1) + Gamma **(adoptado)** | −1233,01 | 29,2 | 0,00 | — |
| AR(2) + Gamma | −1232,75 | 29,1 | +0,26 | 0,31 |
| Paseo aleatorio + Gamma | −1231,34 | 29,0 | +1,67 | 0,67 |
| AR(1) + exógeno exponencial | −1232,83 | 29,3 | +0,19 | 0,39 |

**Las cuatro especificaciones son predictivamente equivalentes.** La diferencia
máxima es de 1,67 unidades de elpd, muy inferior al error estándar de cada elpd
(≈ 29). El AR(1) **no se adopta por superioridad predictiva** —que el LOO-CV no
permite establecer— sino por parsimonia e interpretabilidad frente al AR(2), y
por admitir un nivel explícito frente al paseo aleatorio.

---

## Intervalos seriales evaluados

| Nombre | Distribución | Parámetros | Media (días) | Fuente |
|--------|-------------|------------|--------------|--------|
| **China-Gamma** | Gamma | shape 2,29 · rate 0,36 | **6,4** | **Bi et al. (2020) — seleccionado** |
| Referencia | Gamma | shape 6,5 · rate 0,62 | 10,5 | Mishra et al. (2020) ᵃ |
| China-Normal | Normal truncada en 0 | μ 5,0 · σ 5,2 | 6,6 ᵇ | Xu et al. (2020) |
| Colombia | Gamma | shape 1,96 · rate 0,51 | 3,8 | Estrada-Álvarez et al. (2020) |
| Burkina Faso | Gamma | shape 1,04 · rate 0,18 | 5,8 | Somda et al. (2022) |

ᵃ Mishra et al. escriben `Gamma(6.5, 0.62)` sin declarar la parametrización.
Interpretada como (forma, tasa) implica una media de 10,5 días; el mismo grupo
reporta para esa expresión una media de 6,5 días, coherente con Bi et al. Por eso
el modelo definitivo adopta directamente la distribución de Bi et al.

ᵇ La normal de Xu et al. admite intervalos seriales negativos; al truncarla en
cero se descarta el 16,8 % de su masa y la media efectiva pasa de 5,0 a 6,6 días.

Todos los intervalos se discretizan por integración numérica antes de pasarse a
Stan como vector fijo `g[max_si]`.

---

## Requisitos de software

```r
# R >= 4.2
install.packages(c(
  "cmdstanr", "posterior", "bayesplot", "loo",
  "readxl", "dplyr", "tidyr", "ggplot2",
  "patchwork", "lubridate", "gridExtra", "scales"
))

cmdstanr::install_cmdstan()
```

**Versiones usadas en el análisis original:**

| Fase | CmdStan | cmdstanr |
|------|---------|----------|
| Exploratoria (10 especificaciones) | 2.38.0 | 0.8.0 |
| Definitiva | 2.39.0 | 0.9.0 |

Sistema: Windows 11, R 4.2.2 con Rtools correspondiente.

---

## Orden de ejecución

> Los scripts contienen rutas absolutas del equipo en que se ejecutaron. Antes de
> correrlos hay que ajustar la variable de ruta al inicio de cada archivo.

### Fase 1 — Filtro exploratorio de convergencia

Ajusta 10 modelos: 5 intervalos seriales × 2 distribuciones de observación
(Poisson y NB-2), todos con AR(2).

```r
source("R/fase1_exploratoria/analisis_bogota_ar2_gamma.R")      # ≈ 4-5 h
source("R/fase1_exploratoria/metricas_loo_sensibilidad_AR.R")   # < 1 min
source("R/fase1_exploratoria/grafica_estacionariedad.R")
```

**Resultado:** las cinco especificaciones NB-2 superan los criterios de R̂ y ESS;
ninguna de las cinco Poisson lo hace. Sin un parámetro de dispersión, el proceso
latente debe absorber toda la variabilidad de los conteos: la escala de las
innovaciones se multiplica por 8,6 y la persistencia cae a menos de una cuarta
parte.

### Fase 2 — Modelo definitivo

```r
source("R/fase2_definitiva/06_fase2_china_gamma.R")        # ≈ 30-45 min por modelo
source("R/fase2_definitiva/07_cifras_modelo_definitivo.R")
source("R/fase2_definitiva/08_figuras_modelo_definitivo.R")
source("R/fase2_definitiva/10_verificacion_predictiva.R")
```

Los scripts `01`–`05` contienen las verificaciones previas sobre la corrida
con `adapt_delta = 0.99`: cobertura de los intervalos predictivos,
diagnósticos de la posterior, ajuste de las corridas AR(1) y paseo aleatorio,
comparación de trayectorias y consistencia de las cifras derivadas. El `09`
recoge el análisis exploratorio de intervenciones, que no se reporta entre los
hallazgos y se conserva como respaldo de esa decisión.

---

## Resultados numéricos incluidos

La carpeta `resultados/` contiene las salidas en CSV que respaldan cada tabla y
cifra del documento, entre ellas:

| Archivo | Contenido |
|---------|-----------|
| `f2cg_posteriores.csv` | Posteriores de los cuatro modelos de fase 2 |
| `f2cg_diagnosticos.csv` | R̂, ESS, divergencias, EBFMI |
| `f2cg_loo.csv` | Comparación LOO-CV |
| `f2cg_trayectorias_resumen.csv` | R(t) máximo, fecha de cruce del umbral |
| `ppc_cobertura_sin_represas.csv` | Cobertura de los intervalos predictivos |
| `ppc_represas_reporte.csv` | Días de notificación anómala |
| `paso5_residuos.csv` | Residuos de Pearson y casos diarios |

---

## Referencias

- Bi, Q., *et al.* (2020). Epidemiology and transmission of COVID-19 in 391 cases
  and 1286 of their close contacts in Shenzhen, China. *The Lancet Infectious
  Diseases*, 20(8), 911–919.
- Estrada-Álvarez, J. M., *et al.* (2020). Estimación del intervalo serial y
  número reproductivo básico para los casos importados de COVID-19.
  *Revista de Salud Pública*, 22(2). https://doi.org/10.15446/rsap.v22n2.87492
- Mishra, S., *et al.* (2020). A COVID-19 Model for Local Authorities of the
  United Kingdom. *Imperial College London*. arXiv:2006.16487.
- Somda, S. M. A., Ouedraogo, B., Paré, C. B., & Kouanda, S. (2022). Estimation
  of the Serial Interval and the Effective Reproductive Number of COVID-19
  Outbreak Using Contact Data in Burkina Faso, a Sub-Saharan African Country.
  *Computational and Mathematical Methods in Medicine*.
  https://doi.org/10.1155/2022/8239915
- Xu, X., *et al.* (2020). Reconstruction of transmission pairs for novel
  coronavirus disease 2019 (COVID-19) in mainland China.
- Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
  evaluation using leave-one-out cross-validation and WAIC.
  *Statistics and Computing*, 27(5), 1413–1432.
- Stan Development Team (2024). *Stan Modeling Language Users Guide and
  Reference Manual*.

---

## Licencia

Código bajo licencia [MIT](LICENSE). Los datos epidemiológicos provienen de
fuente oficial del Distrito Capital de Bogotá y su uso debe respetar las
condiciones de la fuente original. El autor no se responsabiliza por
modificaciones realizadas por terceros.
