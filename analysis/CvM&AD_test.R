# ----------------------------------------------
# Apply Cramér–von Mises and AD test to all fitted model
# ----------------------------------------------
# Note: This script requires the parameters of all fitted model. 
# Ensure all modelling scripts have been run first:
#library(here)
#here::i_am("analysis/CvM&AD_test.R")
#source(here("analysis","model_S&P_500_Meixner.R"))
#source(here("analysis","model_S&P_500_GH.R"))
#source(here("analysis","model_S&P_500_Kou_DEJE.R"))
#source(here("analysis","model_S&P_500_CGMY(COS_method).R"))
#source(here("analysis","model_S&P_500_alpha-stable(COS_method).R"))

# --- Import required package and function ---
library(goftest)
source(here("src","AD_test.R"))

# Meixner CDF with fitted parameters
meixner_cdf_fitted <- function(x) meixner_cdf_m(x, 
                                                m = m_Meixner, 
                                                a = alpha_Meixner, 
                                                b = beta_Meixner, 
                                                d = delta_Meixner
                                                )

# Meixner RNG via inverse CDF (quantile function)
meixner_rgen <- function(n) {
  u <- runif(n)
  meixner_qf_m(u)
}

# Apply Cramér–von Mises test to Meixner
set.seed(7914) # For reproducibility
cvm.test(log_returns, "meixner_cdf_m", 
         m=m_Meixner, a=alpha_Meixner, b=beta_Meixner, d=delta_Meixner, 
         estimated=TRUE
         )

# Apply Anderson–Darling test to Meixner
result_meixner <- ad_test_monte_carlo(
  sample   = log_returns,
  cdf_fun  = meixner_cdf_fitted,
  rgen_fun = meixner_rgen,
  B        = 1000,
  seed     = 7914
)

# Print results
cat("Anderson–Darling statistic (Meixner):", result_meixner$ad_statistic, "\n")
cat("p-value:", result_meixner$p_value, "\n")
cat("Critical values (90%, 95%, 99%):\n")
print(result_meixner$critical_values)


# Define CDF and RNG for fitted GH model
gh_cdf <- function(x) pghyp(x, object = gh_fit)
gh_rgen <- function(n) rghyp(n, object = gh_fit)

# Apply Cramér–von Mises test to GH
set.seed(7914) # For reproducibility
cvm.test(log_returns, null = gh_cdf, estimated = TRUE)

# Apply Anderson–Darling test to GH
result_gh <- ad_test_monte_carlo(
  sample = log_returns,
  cdf_fun = gh_cdf,
  rgen_fun = gh_rgen,
  B = 3000,
  seed = 7914
)

# Print results
cat("Anderson–Darling statistic (GH):", result_gh$ad_statistic, "\n")
cat("p-value:", result_gh$p_value, "\n")
cat("Critical values (90%, 95%, 99%):\n")
print(result_gh$critical_values)

# Grid and density
x_vals <- seq(-0.5, 0.5, length.out = 10000)
dens_vals <- cos_density(x_vals, cf_fit, a, b, N = 1024)

# Compute CDF from density
dx <- diff(x_vals)[1]
cdf_vals <- cumsum(dens_vals) * dx
cdf_vals <- cdf_vals / max(cdf_vals)  # normalize to [0,1]

# Define CGMY CDF
cgmy_cdf_func <- function(x) {
  sapply(x, function(xi) {
    if (xi <= min(x_vals)) return(0)
    if (xi >= max(x_vals)) return(1)
    approx(x_vals, cdf_vals, xout = xi, rule = 2)$y
  })
}

# RNG via inverse transform sampling
cgmy_rgen <- function(n) {
  u <- runif(n)
  cgmy_quantile(u)
}


# Apply CvM test
set.seed(7914) # For reproducibility
cvm.test(log_returns, null = cgmy_cdf_func, estimated = TRUE)

# Apply AD test
result_cgmy <- ad_test_monte_carlo(
  sample   = log_returns,
  cdf_fun  = cgmy_cdf_func,
  rgen_fun = cgmy_rgen,
  B        = 3000,
  seed     = 7914
)

# Print results
cat("Anderson–Darling statistic (CGMY):", result_cgmy$ad_statistic, "\n")
cat("p-value:", result_cgmy$p_value, "\n")
cat("Critical values (90%, 95%, 99%):\n")
print(result_cgmy$critical_values)


# Define CDF and RNG for α-stable
stable_cdf <- function(x) pstable(x, alpha_stable, beta_stable, c_stable, mu_stable)
stable_rgen <- function(n) rstable(n, alpha_stable, beta_stable, c_stable, mu_stable)

# Apply CvM test
set.seed(7914) # For reproducibility
cvm.test(log_returns, null = stable_cdf, estimated = TRUE)

# Apply AD test
result_stable <- ad_test_monte_carlo(
  sample   = log_returns,
  cdf_fun  = stable_cdf,
  rgen_fun = stable_rgen,
  B        = 3000,
  seed     = 7914
)

# Print results
cat("Anderson–Darling statistic (α-Stable):", result_stable$ad_statistic, "\n")
cat("p-value:", result_stable$p_value, "\n")
cat("Critical values (90%, 95%, 99%):\n")
print(result_stable$critical_values)

# Define DEJD CDF
dejd_cdf <- function(x) {
  sapply(x, function(xi) {
    if (xi >= 0) {
      1 - p * exp(-eta1 * xi)
    } else {
      (1 - p) * exp(eta2 * xi)
    }
  })
}

# Define RNG function using inverse CDF
dejd_rgen <- function(n) {
  u <- runif(n)
  de_quantile(u, p = p, eta1 = eta1, eta2 = eta2)
}

# Apply CvM test
set.seed(7914) # For reproducibility
cvm.test(log_returns, null = dejd_cdf, estimated = TRUE)

# Apply AD test
result_dejd <- ad_test_monte_carlo(
  sample   = log_returns,
  cdf_fun  = dejd_cdf,
  rgen_fun = dejd_rgen,
  B        = 3000,
  seed     = 7914
)

# Print results
cat("Anderson–Darling statistic (DEJD):", result_dejd$ad_statistic, "\n")
cat("p-value:", result_dejd$p_value, "\n")
cat("Critical values (90%, 95%, 99%):\n")
print(result_dejd$critical_values)


