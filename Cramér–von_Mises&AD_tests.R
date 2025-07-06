library(goftest)

# Apply Cramér–von Mises test to Meixner
cvm.test(log_returns, "meixner_cdf_m", m=m, a=alpha, b=beta, d=delta, estimated=TRUE)

# Apply Anderson–Darling test to Meixner
ad.test(log_returns, "meixner_cdf_m", m=m, a=alpha, b=beta, d=delta, estimated=TRUE)

# Define the GH CDF function wrapper
gh_cdf <- function(x) {
  pghyp(x, object = gh_fit)
}

# Apply Cramér–von Mises test to GH
cvm.test(log_returns, null = gh_cdf, estimated = TRUE)

# Apply Anderson–Darling test to GH
ad.test(log_returns, null = gh_cdf, estimated = TRUE)
