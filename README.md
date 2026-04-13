# Modelo Bayesiano de Renovación COVID-19 — Bogotá

**Tesis de Maestría en Estadística · Universidad Nacional de Colombia**

Implementación de un modelo de ecuaciones de renovación dependiente de la edad
para estimar el número reproductivo efectivo R(t) de COVID-19 en Bogotá
(14 mar – 09 sep 2020).

**Modelo definitivo seleccionado:** AR(1) + componente exógeno Gamma + verosimilitud
Binomial Negativa tipo 2 (NB-2) + intervalo serial de Bi et al. (2020) vía
Mishra et al. (2020).

---

## Estructura del repositorio

```
renewal-model/
├── data/
│   └── datos_agregados.xlsx          # datos de entrada (hojas: NO_IMPORTADOS, IMPORTADOS)
│
├── MODELO_STAN_4/
│   ├── figura_exploratoria_resolucion_diaria_1.R   # figura exploratoria inicial
│   ├── grafica_estacionariedad.R                   # diagnóstico de estacionariedad AR
│   └── EXOGENO_GAMMA_ANALISIS_PRELIMINAR/          # Fase 1: análisis preliminar
│       ├── analisis_bogota_ar2_gamma.R              # ajuste 10 modelos (5 IS × NB2/Poisson)
│       ├── metricas_loo_sensibilidad_AR.R           # comparación LOO y sensibilidad
│       ├── renewal_bogota_nb2_ar2_gamma.stan        # modelo NB-2 + AR(2) + Gamma
│       └── renewal_bogota_poisson_ar2_gamma.stan    # modelo Poisson + AR(2) + Gamma
│
└── MODELO_STAN_DEFINITIVO/                          # Fase 2: modelo definitivo
    ├── analisis_fase2_definitivo.R                  # sensibilidad AR1 vs AR2, Gamma vs Exp
    ├── figura_exogeno.R                             # figura del componente exógeno
    ├── renewal_bogota_nb2_ar1_gamma.stan            # MODELO GANADOR: AR(1) + NB-2 + Gamma
    ├── renewal_bogota_nb2_ar2_gamma.stan            # competidor: AR(2) + NB-2 + Gamma
    └── renewal_bogota_nb2_ar2_exponencial.stan      # competidor: AR(2) + NB-2 + Exponencial
```

> **Nota sobre outputs:** los archivos `.rds` (fits MCMC), `.pdf`/`.png` (figuras),
> `.csv` (tablas intermedias) y `.exe` (binarios Stan compilados) están en `.gitignore`
> porque se regeneran al ejecutar los scripts. Ver sección *Reproducibilidad* abajo.

---

## Datos

**Archivo:** `data/datos_agregados.xlsx`

| Hoja | Contenido |
|------|-----------|
| `NO_IMPORTADOS` | Casos locales diarios de COVID-19 en Bogotá (usados en la verosimilitud) |
| `IMPORTADOS` | Casos importados diarios (componente exógeno del modelo) |

**Período:** 14 de marzo de 2020 – 09 de septiembre de 2020 (180 días).

---

## Requisitos de software

```r
# R >= 4.2.0
install.packages(c(
  "here",       # rutas relativas a la raíz del repositorio
  "cmdstanr",   # interfaz R a CmdStan (MCMC)
  "posterior",  # manejo de draws
  "bayesplot",  # diagnósticos y figuras
  "loo",        # comparación de modelos (LOO-CV)
  "readxl",     # lectura de datos Excel
  "dplyr",
  "ggplot2",
  "tidyr",
  "patchwork",  # composición de figuras (figura_exogeno.R)
  "lubridate",  # manejo de fechas (figura exploratoria)
  "gridExtra",  # composición de figuras (figura exploratoria)
  "scales"      # formateo de ejes
))

# CmdStan (compilador Stan) — instalar vía cmdstanr:
cmdstanr::install_cmdstan()
```

**Versiones usadas en el análisis original:**
- R 4.4.x, CmdStan 2.35.x, cmdstanr 0.8.x
- Sistema: Windows 11 con Rtools 4.4

---

## Orden de ejecución

### Paso 0 — Datos exploratorios (opcional, para entender los datos)

```r
# Figura exploratoria de casos diarios
source("MODELO_STAN_4/figura_exploratoria_resolucion_diaria_1.R")

# Diagnóstico de estacionariedad del proceso AR
source("MODELO_STAN_4/grafica_estacionariedad.R")
```

### Paso 1 — Fase 1: Análisis preliminar de convergencia

> **Tiempo estimado: 4–5 horas.** Se ajustan 10 modelos en secuencia,
> cada uno con 4 cadenas paralelas × (1000 warmup + 1000 sampling).
> **Se recomienda correr overnight.**

Ajusta **10 modelos** AR(2)+Gamma con 5 especificaciones del intervalo serial (IS):

| IS | Fuente |
|----|--------|
| Referencia | Bi et al. (2020) vía Mishra et al. (2020) — **IS seleccionado** |
| Burkina Faso | Somda et al. (2022) |
| China gamma | Li et al. (2020) — parametrización Gamma |
| China normal | Li et al. (2020) — parametrización Normal |
| Colombia | Estrada-Álvarez et al. (2020) |

```r
# 1a. Ajustar los 10 modelos — correr overnight (≈ 4-5 h)
source("MODELO_STAN_4/EXOGENO_GAMMA_ANALISIS_PRELIMINAR/analisis_bogota_ar2_gamma.R")

# 1b. Comparación LOO y análisis de sensibilidad al IS (< 1 min)
source("MODELO_STAN_4/EXOGENO_GAMMA_ANALISIS_PRELIMINAR/metricas_loo_sensibilidad_AR.R")
```

**Outputs generados en** `MODELO_STAN_4/EXOGENO_GAMMA_ANALISIS_PRELIMINAR/`:
- `fit_gamma_*.rds` / `draws_gamma_*.rds` — fits MCMC y draws por modelo
- `diagnosticos_convergencia_gamma.csv` — R-hat, ESS, divergencias
- `loo_comparison_gamma.csv` — tabla de comparación LOO
- `fig1_traceplots_gamma.pdf` … `fig6_posterior_params_gamma.pdf` — figuras diagnósticas

### Paso 2 — Fase 2: Modelo definitivo y análisis de sensibilidad

> **Tiempo estimado: 30–45 minutos.** Se ajusta 1 modelo definitivo AR(1)
> con 4 cadenas × (1500 warmup + 1500 sampling), más 2 modelos de
> sensibilidad con la misma configuración del Paso 1.

Compara AR(1) vs AR(2) y prior Gamma vs Exponencial para el exógeno.
Usa el IS de referencia (Bi et al. 2020) en todos los modelos.

```r
# 2a. Ajustar modelos de sensibilidad y modelo definitivo (≈ 30-45 min)
source("MODELO_STAN_DEFINITIVO/analisis_fase2_definitivo.R")

# 2b. Figura del componente exógeno (requiere fit del paso 2a)
source("MODELO_STAN_DEFINITIVO/figura_exogeno.R")
```

**Outputs generados en** `MODELO_STAN_DEFINITIVO/`:
- `fit_fase2_M1_AR1_Gamma.rds` — modelo AR(1)+Gamma
- `fit_fase2_M4_AR2_Exp.rds` — modelo AR(2)+Exponencial
- `fit_AR1_Gamma_definitivo.rds` / `draws_AR1_Gamma_definitivo.rds` — modelo definitivo final
- `loo_fase2_comparacion.csv` — comparación LOO Fase 2
- `posterior_modelo_definitivo.csv` — resumen posterior de parámetros
- `fig1_traceplots_definitivo.pdf` … `fig5_exogeno.pdf` — figuras del modelo final

---

## Intervalos seriales usados

| Nombre | Distribución | Parámetros | Media (días) | Fuente |
|--------|-------------|------------|--------------|--------|
| **Referencia** | Gamma | shape=6.5, rate=0.62 | 10.5 | Bi et al. (2020) vía Mishra et al. (2020) |
| China gamma | Gamma | shape=2.29, rate=0.36 | 6.4 | Li et al. (2020) |
| China normal | Normal truncada en 0 | μ=5.0, σ=5.2 | ~5.0 | Li et al. (2020) |
| Burkina Faso | Gamma | shape=1.04, rate=0.18 | 5.8 | Somda et al. (2022) |
| Colombia | Gamma | shape=1.96, rate=0.51 | 3.8 | Estrada-Álvarez et al. (2020) |

> Todos los intervalos seriales se discretizan mediante integración numérica
> antes de pasarlos a Stan como vector fijo `SI[N]` (longitud máxima = 30 días).

---

## Modelo estadístico (resumen)

El número de casos en el día t sigue:

```
y(t) ~ NB-2(f(t), φ)

f(t) = μ(t) + R_t · Σ_{s=1}^{t-1} y(s) · g(t-s)

log R_t = μ_R + ρ₁ · (log R_{t-1} − μ_R) + ε_t,   ε_t ~ N(0, σ_ε²)   [AR(1)]

μ(t) ~ Gamma(1.401, 0.168)   [exógeno: casos importados, activo días 1-12]
```

donde `g(s)` es el intervalo serial discretizado y `φ` es el parámetro de
sobredispersión de la NB-2. El componente exógeno `μ(t)` modela los casos
importados durante el período pre-cuarentena (14–25 marzo 2020) y está
calibrado por método de momentos sobre los datos observados
(media = 8.33 casos/día, varianza = 49.52).

**Priors del modelo definitivo:**

| Parámetro | Prior | Interpretación |
|-----------|-------|----------------|
| `μ_R` | N(0, 0.3) | Media de largo plazo de log(R_t); centrado en R_t ≈ 1 |
| `ρ₁` | N(0.7, 0.15) truncado en [0,1] | Autocorrelación diaria de log(R_t) |
| `σ_ε` | N⁺(0, 0.2) | Desviación estándar de innovaciones diarias |
| `φ` | N⁺(0, 5) | Sobredispersión NB-2 |

**Resultados posteriores del modelo definitivo:**

| Parámetro | Media | IC 95% |
|-----------|-------|--------|
| μ_R | 0.203 | [−0.258, 0.552] |
| ρ₁ | 0.962 | [0.870, 0.997] |
| σ_ε | 0.099 | [0.047, 0.180] |
| φ | 3.263 | [2.609, 4.015] |

**Diagnósticos MCMC:**

| Métrica | Valor |
|---------|-------|
| R̂ máximo | 1.0036 |
| ESS mínimo | 1,718 |
| Divergencias | 120 / 6,000 (2.0%) |
| Treedepth saturadas | 0 (0.0%) |

---

## Ajuste de rutas antes de correr

Los scripts usan rutas absolutas. Antes de ejecutar, editar en cada script la
variable `RUTA_BASE` al inicio del archivo:

```r
# Reemplazar con la ruta local al repositorio clonado:
RUTA_BASE <- "C:/tu_ruta/renewal-model"
```

---

## Descripción de archivos Stan

| Archivo | Modelo |
|---------|--------|
| `renewal_bogota_nb2_ar2_gamma.stan` | AR(2) + exógeno Gamma + NB-2 (análisis preliminar) |
| `renewal_bogota_poisson_ar2_gamma.stan` | AR(2) + exógeno Gamma + Poisson (análisis preliminar) |
| `renewal_bogota_nb2_ar1_gamma.stan` | **AR(1) + exógeno Gamma + NB-2 (MODELO DEFINITIVO)** |
| `renewal_bogota_nb2_ar2_gamma.stan` | AR(2) + exógeno Gamma + NB-2 (sensibilidad Fase 2) |
| `renewal_bogota_nb2_ar2_exponencial.stan` | AR(2) + exógeno Exponencial + NB-2 (sensibilidad Fase 2) |

---

## Referencias

- Bi Q, et al. (2020). Epidemiology and transmission of COVID-19 in 391 cases and 1286 of their close contacts in Shenzhen, China. *Lancet Infect Dis*, 20(8), 911–919.
- Estrada-Álvarez JM, et al. (2020). Estimación del intervalo serial y número reproductivo básico para los casos importados de COVID-19. *Rev Salud Pública*, 22(2). https://doi.org/10.15446/rsap.v22n2.87492
- Li Q, et al. (2020). Early Transmission Dynamics in Wuhan, China, of Novel Coronavirus–Infected Pneumonia. *NEJM*, 382(13), 1199–1207.
- Mishra S, et al. (2020). A COVID-19 Model for Local Authorities of the United Kingdom. *Imperial College London*, arXiv:2006.16487.
- Somda ZC, et al. (2022). Serial interval of COVID-19 in Burkina Faso. *[fuente por completar]*.
- Stan Development Team (2024). *Stan Modeling Language Users Guide and Reference Manual*.

---

## Datos

Los datos (`data/datos_agregados.xlsx`) están incluidos en el repositorio.
Fuente: Sistema de Vigilancia en Salud Pública del Distrito Capital de Bogotá.

---

## Licencia

Código bajo licencia [MIT](LICENSE). El autor no se responsabiliza por
modificaciones realizadas por terceros. Los datos epidemiológicos son de
fuente oficial del Distrito Capital de Bogotá y su uso debe respetar
las condiciones de la fuente original.
