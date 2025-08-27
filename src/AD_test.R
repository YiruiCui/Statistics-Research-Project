# This file contain the function that compute critical values for AD test

ad_test_monte_carlo <- function(sample,
                                cdf_fun,    # function(x)  CDF
                                rgen_fun,   # function(n)  samples
                                B = 1000,   # Monte Carlo sample count
                                n = 100,
                                plot = TRUE,
                                probs = c(0.90, 0.95, 0.99),
                                seed = 7914) {
  
  # Internal: compute Anderson–Darling statistic for U(0,1) transformed data
  ad_statistic_uniform <- function(u) {
    u <- sort(u)
    n <- length(u)
    i <- seq_len(n)
    A2 <- -n - mean((2 * i - 1) * (log(u) + log(1 - rev(u))))
    return(A2)
  }
  
  # Step 1: transform real sample to uniform
  u_vals <- cdf_fun(sample)
  obs_ad <- ad_statistic_uniform(u_vals)
  
  # Step 2: simulate B AD statistics under null model
  set.seed(seed)
  ad_stats_sim <- replicate(B, {
    sim_data <- rgen_fun(n)
    u_sim <- cdf_fun(sim_data)
    ad_statistic_uniform(u_sim)
  })
  
  # Step 3: compute p-value and critical values
  p_value <- mean(ad_stats_sim >= obs_ad)
  crit_vals <- quantile(ad_stats_sim, probs = probs)
  
  # Step 4: optional plot
  if (plot) {
    hist(ad_stats_sim, breaks = 50, col = "lightgray",
         main = "Monte Carlo Null of AD Test",
         xlab = "AD Statistic")
    abline(v = obs_ad, col = "red", lwd = 2)
    legend("topright", legend = sprintf("Observed AD = %.4f\np-value = %.4f", obs_ad, p_value),
           col = "red", lwd = 2, bty = "n")
  }
  
  # Return result as list
  return(list(
    ad_statistic = obs_ad,
    p_value = p_value,
    critical_values = crit_vals,
    simulated_stats = ad_stats_sim
  ))
}
