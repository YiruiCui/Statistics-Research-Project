# ----------------------------------------------
# Apply Chi Squared test to all fitted model
# ----------------------------------------------
# Note: This script requires the parameters of all fitted model. 
# Ensure all modelling scripts have been run first:
#library(here)
#here::i_am("analysis/chi2_test.R")
#source(here("analysis","model_S&P_500_Meixner.R"))
#source(here("analysis","model_S&P_500_GH.R"))
#source(here("analysis","model_S&P_500_Kou_DEJE.R"))
#source(here("analysis","model_S&P_500_CGMY(COS_method).R"))
#source(here("analysis","model_S&P_500_alpha-stable(COS_method).R"))

# Divide log-returns into bins 
breaks <- quantile(log_returns, probs = seq(0, 1, length.out = 101))  # 100 bins
observed_counts <- hist(log_returns, breaks = breaks, plot = FALSE)$counts

set.seed(7914) # For reproducibility
# Compute expected counts from fitted Meixner
midpoints <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
expected_probs <- sapply(1:100, function(i) {
  integrate(meixner_pdf_m, lower = breaks[i], upper = breaks[i+1],
            m = m_Meixner, a = alpha_Meixner, b = beta_Meixner, d = delta_Meixner)$value
})

expected_counts <- expected_probs * length(log_returns)

# Apply Chi-squared test
chisq_mexiner <- chisq.test(x = observed_counts, p = expected_probs, rescale.p = TRUE)

# Extract the chi-squared statistic from the test result
chi_statistic <- chisq_mexiner$statistic

# Calculate the correct degrees of freedom
df <- length(observed_counts) - 1 - 4 # k - 1 - p

# Calculate the correct p-value
p_value <- pchisq(chi_statistic, df = df, lower.tail = FALSE)

# Print the corrected results
cat("Chi-squared Statistic:", chi_statistic, "\n")
cat("Correct Degrees of Freedom:", df, "\n")
cat("Correct p-value:", p_value, "\n")


set.seed(7914) # For reproducibility
# Compute expected probabilities using fitted GH model
# Get cumulative probabilities at bin edges
gh_cdf_vals <- sapply(breaks, function(x) pghyp(x, object = gh_fit))

# Compute bin probabilities as CDF differences
expected_probs <- diff(gh_cdf_vals)
expected_counts <- expected_probs * length(log_returns)

# Apply Chi-squared test ---
chisq_gh <- chisq.test(x = observed_counts, p = expected_probs, rescale.p = TRUE)

# Extract the chi-squared statistic from the test result
chi_statistic <- chisq_gh$statistic

# Calculate the correct degrees of freedom
df <- length(observed_counts) - 1 - 5 # k - 1 - p

# Calculate the correct p-value
p_value <- pchisq(chi_statistic, df = df, lower.tail = FALSE)

# Print the corrected results
cat("Chi-squared Statistic:", chi_statistic, "\n")
cat("Correct Degrees of Freedom:", df, "\n")
cat("Correct p-value:", p_value, "\n")

set.seed(7914) # For reproducibility
# Get cumulative probabilities using COS-derived CDF
cgmy_cdf_vals <- sapply(breaks, function(x) {
  if (x <= min(x_vals)) return(0)
  if (x >= max(x_vals)) return(1)
  approx(x_vals, cdf_vals, xout = x, rule = 2)$y
})
expected_probs <- diff(cgmy_cdf_vals)
expected_counts <- expected_probs * length(log_returns)

# Apply Chi-squared test
chisq_cgmy <- chisq.test(observed_counts, p = expected_probs, rescale.p = TRUE)

# Extract the chi-squared statistic from the test result
chi_statistic <- chisq_cgmy$statistic

# Calculate the correct degrees of freedom
df <- length(observed_counts) - 1 - 5 # k - 1 - p

# Calculate the correct p-value
p_value <- pchisq(chi_statistic, df = df, lower.tail = FALSE)

# Print the corrected results
cat("Chi-squared Statistic:", chi_statistic, "\n")
cat("Correct Degrees of Freedom:", df, "\n")
cat("Correct p-value:", p_value, "\n")

set.seed(7914) # For reproducibility
# Compute stable CDF at breakpoints
stable_cdf_vals <- sapply(breaks, function(x) {
  if (x <= min(x_vals)) return(0)
  if (x >= max(x_vals)) return(1)
  approx(x_vals, cdf_vals, xout = x, rule = 2)$y
})
expected_probs <- diff(stable_cdf_vals)
expected_counts <- expected_probs * length(log_returns)

# Apply Chi-squared test
chisq_stable <- chisq.test(observed_counts, p = expected_probs, rescale.p = TRUE)

# Extract the chi-squared statistic from the test result
chi_statistic <- chisq_stable$statistic

# Calculate the correct degrees of freedom
df <- length(observed_counts) - 1 - 4 # k - 1 - p

# Calculate the correct p-value
p_value <- pchisq(chi_statistic, df = df, lower.tail = FALSE)

# Print the corrected results
cat("Chi-squared Statistic:", chi_statistic, "\n")
cat("Correct Degrees of Freedom:", df, "\n")
cat("Correct p-value:", p_value, "\n")

# CDF of DEJD
de_cdf <- function(x, p, eta1, eta2) {
  ifelse(x >= 0,
         1 - p * exp(-eta1 * x),
         (1 - p) * exp(eta2 * x))
}

set.seed(7914) # For reproducibility
# Compute cumulative probabilities at bin edges
dejd_cdf_vals <- sapply(breaks, function(x) de_cdf(x, p, eta1, eta2))
expected_probs <- diff(dejd_cdf_vals)
expected_counts <- expected_probs * length(log_returns)

# Apply Chi-squared test
chisq_dejd <- chisq.test(observed_counts, p = expected_probs, rescale.p = TRUE)

# Extract the chi-squared statistic from the test result
chi_statistic <- chisq_dejd$statistic

# Calculate the correct degrees of freedom
df <- length(observed_counts) - 1 - 3 # k - 1 - p

# Calculate the correct p-value
p_value <- pchisq(chi_statistic, df = df, lower.tail = FALSE)

# Print the corrected results
cat("Chi-squared Statistic:", chi_statistic, "\n")
cat("Correct Degrees of Freedom:", df, "\n")
cat("Correct p-value:", p_value, "\n")

