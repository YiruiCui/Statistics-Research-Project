# Meixner CDF via numerical integration of PDF
meixner_cdf <- function(q, m, a, b, d) {
  sapply(q, function(x) {
    integrate(meixner_pdf_m, lower = -10, upper = x, m = m, a = a, b = b, d = d)$value
  })
}
log_returns_jittered <- jitter(log_returns, amount = 1e-8)
# KS Test: log-returns vs fitted Meixner CDF
ks_result <- ks.test(log_returns_jittered, function(x) meixner_cdf(x, m = m, a = alpha, b = beta, d = delta))

# Print result
print(ks_result)
