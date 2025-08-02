library(goftest)

# Fixed-parameter Meixner CDF
meixner_cdf_fixed <- function(x) meixner_cdf_m(x, m = m, a = alpha, b = beta, d = delta)

# Meixner RNG via inverse CDF (quantile function)
meixner_rgen <- function(n) {
  u <- runif(n)
  meixner_qf_m(u)
}

# Apply Cramér–von Mises test to Meixner
cvm.test(log_returns, "meixner_cdf_m", m=m, a=alpha, b=beta, d=delta, estimated=TRUE)

# Apply Anderson–Darling test to Meixner
# Run AD test for Meixner
result_meixner <- ad_test_monte_carlo(
  sample   = log_returns,
  cdf_fun  = meixner_cdf_fixed,
  rgen_fun = meixner_rgen,
  B        = 1000,
  seed     = 7914
)

# Print results
cat("Anderson–Darling statistic (Meixner):", result_meixner$ad_statistic, "\n")
cat("p-value:", result_meixner$p_value, "\n")
cat("Critical values (90%, 95%, 99%):\n")
print(result_meixner$critical_values)


# Define CDF and RNG from fitted GH model
gh_cdf <- function(x) pghyp(x, object = gh_fit)
gh_rgen <- function(n) rghyp(n, object = gh_fit)

# Apply Cramér–von Mises test to GH
cvm.test(log_returns, null = gh_cdf, estimated = TRUE)

# Apply Anderson–Darling test to GH
result_gh <- ad_test_monte_carlo(
  sample = log_returns,
  cdf_fun = gh_cdf,
  rgen_fun = gh_rgen,
  B = 3000,
  seed = 7914
)

# View output
result_gh$p_value
result_gh$critical_values


# Define CGMY CDF wrapper from COS table
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


# CvM test
cvm.test(log_returns, null = cgmy_cdf_func, estimated = TRUE)

# AD test
result_cgmy <- ad_test_monte_carlo(
  sample   = log_returns,
  cdf_fun  = cgmy_cdf_func,
  rgen_fun = cgmy_rgen,
  B        = 3000,
  seed     = 7914
)

# Show results
cat("Anderson–Darling statistic (CGMY):", result_cgmy$ad_statistic, "\n")
cat("p-value:", result_cgmy$p_value, "\n")
cat("Critical values (90%, 95%, 99%):\n")
print(result_cgmy$critical_values)


# Define CDF and RNG for α-stable
stable_cdf <- function(x) pstable(x, alpha_stable, beta_stable, c_stable, mu_stable)
stable_rgen <- function(n) rstable(n, alpha_stable, beta_stable, c_stable, mu_stable)

# CvM test
cvm.test(log_returns, null = stable_cdf, estimated = TRUE)

# AD test
result_stable <- ad_test_monte_carlo(
  sample   = log_returns,
  cdf_fun  = stable_cdf,
  rgen_fun = stable_rgen,
  B        = 3000,
  seed     = 7914
)

# Output results
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


# CvM test
cvm.test(log_returns, null = dejd_cdf, estimated = TRUE)

# AD test
result_dejd <- ad_test_monte_carlo(
  sample   = log_returns,
  cdf_fun  = dejd_cdf,
  rgen_fun = dejd_rgen,
  B        = 3000,
  seed     = 7914
)

# Output the results
cat("Anderson–Darling statistic (DEJD):", result_dejd$ad_statistic, "\n")
cat("p-value:", result_dejd$p_value, "\n")
cat("Critical values (90%, 95%, 99%):\n")
print(result_dejd$critical_values)


