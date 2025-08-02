library(KSgeneral)
library(twosamples)
# --- Accept-Reject Sampling with AcceptReject package ---
set.seed(7914)
meixner_samples <- AcceptReject::accept_reject(
  n = 10000L,
  f = meixner_pdf_m,
  continuous = TRUE,
  args_f = list(m=m, a=alpha, b=beta, d=delta),
  xlim = c(-100, 100)
)

# Computes p-value of two-sided KS test: log-returns vs fitted Meixner samples
KS2sample(log_returns, meixner_samples, conservative = T)

# Two-sample Kuiper test
kuiper_mexiner = kuiper_test(log_returns, meixner_samples)
summary(kuiper_mexiner)
#Kuiper2sample(log_returns, meixner_samples, conservative = T)

gh_samples <- rghyp(10000, 
                    object = ghyp(lambda = gh_params[1], 
                                  alpha = gh_params[2], 
                                  mu = gh_params[3], 
                                  sigma = gh_params[4], 
                                  gamma = gh_params[5])
                    )
# Computes p-value of two-sided KS test: log-returns vs fitted GH samples
KS2sample(log_returns, gh_samples, conservative = F)

# Two-sample Kuiper test
kuiper_gh = kuiper_test(log_returns, gh_samples)
summary(kuiper_gh)

# --- CGMY inverse CDF already built from your COS fit ---
cgmy_quantile <- approxfun(cdf_vals, x_vals, rule = 2)

# Sample from CGMY via inverse CDF
set.seed(7914)
cgmy_samples <- cgmy_quantile(runif(10000))

# KS test: log-returns vs CGMY
KS2sample(log_returns, cgmy_samples, conservative = TRUE)

# Kuiper test
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

# Sample from α-stable via inverse CDF
set.seed(7914)
stable_samples <- stable_quantile(runif(10000))

# KS test
KS2sample(log_returns, stable_samples, conservative = TRUE)

# Kuiper test
kuiper_stable <- kuiper_test(log_returns, stable_samples)
summary(kuiper_stable)

# KS test: log-returns vs Kou-DEJD
dej_samples <- de_quantile(runif(10000), p, eta1, eta2)
ks_dejd <- KS2sample(log_returns, dej_samples, conservative = TRUE)
print(ks_dejd)

# Kuiper test
kuiper_dejd <- kuiper_test(log_returns, dej_samples)
summary(kuiper_dejd)
