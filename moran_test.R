# Moran Goodness-of-Fit Test (Cheng & Stephens 1989)
moran_test <- function(u_sorted, k) {
  n <- length(u_sorted)
  m <- n + 1
  gamma <- 0.5772156649  # Euler–Mascheroni constant
  
  # Compute spacings
  u_ext <- c(0, u_sorted, 1)
  spacings <- diff(u_ext)
  spacings[spacings == 0] <- 1e-12  # apply tie-adjusted method
  
  # Moran statistic
  M <- -sum(log(spacings))
  
  # Compute mu and sigma2
  mu <- m * (log(m) + gamma) - 1/2 - 1/(12*m)
  sigma2 <- (pi^2 / 6 - 1) * m - 1/2 - 1/(6*m)
  
  # Compute C1 and C2 (Cheng & Stephens 1989)
  C1 <- mu - sqrt(n/2) * sigma2
  C2 <- sigma2 / sqrt(2*n)
  
  # Standardized test statistic T
  T <- (M + k/2 - C1) / C2
  p_value <- 1 - pchisq(T, df = n)
  
  return(list(
    Moran_statistic = M,
    T_statistic = T,
    df = k,
    p_value = p_value
  ))
}

# Apply Moran test to Meixner fit
u_Meixner <- sort(meixner_cdf_m(log_returns,m,alpha,beta,delta))
result_Meixner <- moran_test(u_sorted = u_Meixner, k = 4)
print(result_Meixner)

# Apply Moran test to GH fit
u_gh <- sort(pghyp(log_returns, object = gh_fit))
result_gh <- moran_test(u_sorted = u_gh, k = 5)
print(result_gh)
