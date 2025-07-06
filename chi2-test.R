# Binning
breaks <- quantile(log_returns, probs = seq(0, 1, length.out = 101))  # 100 bins
observed_counts <- hist(log_returns, breaks = breaks, plot = FALSE)$counts

# Compute expected counts from fitted Meixner
midpoints <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
expected_probs <- sapply(1:100, function(i) {
  integrate(meixner_pdf_m, lower = breaks[i], upper = breaks[i+1],
            m = m, a = alpha, b = beta, d = delta)$value
})

expected_counts <- expected_probs * length(log_returns)

# Chi-squared test
chisq_mexiner <- chisq.test(x = observed_counts, p = expected_probs, rescale.p = TRUE)
print(chisq_mexiner)

# Divide log-returns into bins 
breaks <- quantile(log_returns, probs = seq(0, 1, length.out = 101))  # 10 bins
observed_counts <- hist(log_returns, breaks = breaks, plot = FALSE)$counts

# Compute expected probabilities using fitted GH model
# Get cumulative probabilities at bin edges
gh_cdf_vals <- sapply(breaks, function(x) pghyp(x, object = gh_fit))

# Compute bin probabilities as CDF differences
expected_probs <- diff(gh_cdf_vals)
expected_counts <- expected_probs * length(log_returns)

# Chi-squared test ---
chisq_gh <- chisq.test(x = observed_counts, p = expected_probs, rescale.p = TRUE)
print(chisq_gh)

