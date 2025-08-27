# ----------------------------------------------
# Apply KS and Kuiper test to all fitted model
# ----------------------------------------------
# Note: This script requires the parameters of all fitted model. 
# Ensure all modelling scripts have been run first:
#library(here)
#here::i_am("analysis/KS&Kuiper_test.R")
#source(here("analysis","model_S&P_500_Meixner.R"))
#source(here("analysis","model_S&P_500_GH.R"))
#source(here("analysis","model_S&P_500_Kou_DEJE.R"))
#source(here("analysis","model_S&P_500_CGMY(COS_method).R"))
#source(here("analysis","model_S&P_500_alpha-stable(COS_method).R"))

# --- Import required package ---
library(KSgeneral)
library(twosamples)

# --- Simulate samples of Meixner ---
set.seed(7914) # For reproducibility
meixner_samples <- AcceptReject::accept_reject(
  n = 10000L,
  f = meixner_pdf_m,
  continuous = TRUE,
  args_f = list(m=m_Meixner, a=alpha_Meixner, b=beta_Meixner, d=delta_Meixner),
  xlim = c(-100, 100)
)

# Computes p-value of two-sided KS test: log-returns vs fitted Meixner samples
KS2sample(log_returns, meixner_samples, conservative = T)

# Apply Two-sample Kuiper test
kuiper_mexiner = kuiper_test(log_returns, meixner_samples)
summary(kuiper_mexiner)

# --- Simulate samples of GH ---
set.seed(7914) # For reproducibility
gh_samples <- rghyp(10000, 
                    object = ghyp(lambda = gh_params[1], 
                                  alpha = gh_params[2], 
                                  mu = gh_params[3], 
                                  sigma = gh_params[4], 
                                  gamma = gh_params[5])
                    )

# Computes p-value of two-sided KS test: log-returns vs fitted GH samples
KS2sample(log_returns, gh_samples, conservative = F)

# Apply Two-sample Kuiper test
kuiper_gh = kuiper_test(log_returns, gh_samples)
summary(kuiper_gh)

# --- Simulate samples of CGMY ---
set.seed(7914) # For reproducibility
cgmy_samples <- cgmy_quantile(runif(10000))

# Apply KS test
KS2sample(log_returns, cgmy_samples, conservative = TRUE)

# Apply Kuiper test
kuiper_cgmy <- kuiper_test(log_returns, cgmy_samples)
summary(kuiper_cgmy)

# --- COS inversion for α-stable ---
cf_stable <- function(u) stable_cf(u, alpha_stable, c_stable, beta_stable, mu_stable)

# Compute density and CDF via COS
x_vals <- seq(-1, 1, length.out = 10000)
dens_vals <- cos_density(x_vals, cf_stable, a = -1, b = 1, N = 1024)
dx <- diff(x_vals)[1]
cdf_vals <- cumsum(dens_vals) * dx
cdf_vals <- cdf_vals / max(cdf_vals)
stable_quantile <- approxfun(cdf_vals, x_vals, rule = 2)

# --- Simulate samples of α-stable ---
set.seed(7914) # For reproducibility
stable_samples <- stable_quantile(runif(10000))

# Apply KS test
KS2sample(log_returns, stable_samples, conservative = TRUE)

# Apply Kuiper test
kuiper_stable <- kuiper_test(log_returns, stable_samples)
summary(kuiper_stable)

# --- Simulate samples of Kou-DEJE ---
set.seed(7914) # For reproducibility
dej_samples <- de_quantile(runif(10000), p, eta1, eta2)

# Apply KS test
KS2sample(log_returns, dej_samples, conservative = TRUE)

# Apply Kuiper test
kuiper_dejd <- kuiper_test(log_returns, dej_samples)
summary(kuiper_dejd)
