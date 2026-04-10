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

**Archivo:** `MODELO_STAN_4/datos_agregados.xlsx`

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

Ajusta **10 modelos** AR(2)+Gamma con 5 especificaciones del intervalo serial (IS):

| IS | Fuente |
|----|--------|
| Referencia | Bi et al. (2020) vía Mishra et al. (2020) — **IS final** |
| Burkina | Datos locales |
| China gamma | Li et al. (2020) — parámetros gamma |
| China normal | Li et al. (2020) — parámetros normales |
| Colombia | Datos Colombia |

```r
# 1a. Ajustar los 10 modelos (≈ 30-60 min con 4 cadenas paralelas)
source("MODELO_STAN_4/EXOGENO_GAMMA_ANALISIS_PRELIMINAR/analisis_bogota_ar2_gamma.R")

# 1b. Comparación LOO y análisis de sensibilidad al IS
source("MODELO_STAN_4/EXOGENO_GAMMA_ANALISIS_PRELIMINAR/metricas_loo_sensibilidad_AR.R")
```

**Outputs generados en** `MODELO_STAN_4/EXOGENO_GAMMA_ANALISIS_PRELIMINAR/`:
- `fit_gamma_*.rds` / `draws_gamma_*.rds` — fits MCMC y draws por modelo
- `diagnosticos_convergencia_gamma.csv` — R-hat, ESS, divergencias
- `loo_comparison_gamma.csv` — tabla de comparación LOO
- `fig1_traceplots_gamma.pdf` … `fig6_posterior_params_gamma.pdf` — figuras diagnósticas

### Paso 2 — Fase 2: Modelo definitivo y análisis de sensibilidad

Compara AR(1) vs AR(2) y prior Gamma vs Exponencial para el exógeno.
Usa el IS de referencia (Bi et al. 2020) en todos los modelos.

```r
# 2a. Ajustar modelos de sensibilidad y modelo definitivo (≈ 20-40 min)
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

| Nombre | Distribución | Media (días) | Fuente |
|--------|-------------|--------------|--------|
| **Referencia** | Gamma(μ=4.7, σ=2.9) | 4.7 | Bi et al. (2020) vía Mishra et al. (2020) |
| Burkina | Empírico | — | Datos locales |
| China gamma | Gamma | — | Li et al. (2020) |
| China normal | Normal | — | Li et al. (2020) |
| Colombia | Empírico | — | Datos Colombia |

---

## Modelo estadístico (resumen)

El número de casos en el día t sigue:

```
y_t ~ NB-2(μ_t, φ)

μ_t = R_t · Σ_{s=1}^{t-1} y_{t-s} · w_s  +  λ_t

log R_t = α · log R_{t-1} + ε_t,    ε_t ~ Normal(0, σ²)   [AR(1)]

λ_t ~ Gamma(α_λ, β_λ)   [exógeno: casos importados]
```

donde `w_s` es el intervalo serial discretizado y `φ` es el parámetro de
sobredispersión de la NB-2.

**Priors:**
- `α ~ Uniform(-1, 1)` (estacionariedad AR)
- `σ ~ HalfNormal(0, 0.5)`
- `R_1 ~ LogNormal(log(2.5), 0.5)`
- `φ ~ HalfNormal(0, 10)`

---

## Ajuste de rutas antes de correr

Los scripts usan rutas absolutas. Antes de ejecutar, editar en cada script la
sección `SECCIÓN 1: CONFIGURACIÓN`:

```r
# En 01_analisis_bogota_ar2_gamma.R y 01_analisis_fase2_definitivo.R:
RUTA_DATOS  <- "<tu_ruta>/data/datos_agregados.xlsx"
RUTA_STAN   <- "<tu_ruta>/MODELO_STAN_4/EXOGENO_GAMMA_ANALISIS_PRELIMINAR"
# etc.
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
- Mishra S, et al. (2020). A COVID-19 Model for Local Authorities of the United Kingdom. *Imperial College London*, Report 33.
- Li Q, et al. (2020). Early Transmission Dynamics in Wuhan, China, of Novel Coronavirus–Infected Pneumonia. *NEJM*, 382(13), 1199–1207.
- Stan Development Team (2024). *Stan Modeling Language Users Guide and Reference Manual*.
- Gabry J, et al. (2024). *cmdstanr: R Interface to CmdStan*.

---

## Licencia

Código bajo licencia MIT. Datos epidemiológicos de fuente oficial del Distrito Capital de Bogotá.
