library(loo)
library(cmdstanr)
library(here)

RUTA_OUTPUT <- here("MODELO_STAN_4", "EXOGENO_GAMMA_ANALISIS_PRELIMINAR")

N_CHAINS   <- 4
N_SAMPLING <- 1000

si_nombres <- c("referencia", "china_norm", "china_gamma", "colombia", "burkina")

loo_lista <- list()
for (nm in si_nombres) {
  fit   <- readRDS(file.path(RUTA_OUTPUT, paste0("fit_gamma_", nm, "_NB2.rds")))
  ll    <- fit$draws("log_lik", format="matrix")
  r_eff <- relative_eff(exp(ll),
                        chain_id = rep(1:N_CHAINS, each = N_SAMPLING))
  loo_lista[[nm]] <- loo(ll, r_eff=r_eff)
  cat(sprintf("%-15s elpd=%.1f (SE=%.1f)\n", nm,
              loo_lista[[nm]]$estimates["elpd_loo","Estimate"],
              loo_lista[[nm]]$estimates["elpd_loo","SE"]))
}

comp <- loo_compare(loo_lista)
print(comp)
