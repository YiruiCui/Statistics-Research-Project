library(here)
library(tidyverse)
library(ghyp)
library(gridExtra)
library(SuppDists) # For rinvGauss function

set.seed(7914)

# Identify project location
here::i_am("analysis/model_SP500_BNS-IG(MLE-APF).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- 2. Auxiliary Particle Filter for BNS with IG-SV ---
apf_log_likelihood_ig <- function(params, data, n_particles) {
  # Unpack parameters
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  a_ig      <- params[4]
  b_ig      <- params[5]
  sigma0_sq <- params[6]
  
  n_obs <- length(data)
  log_likelihood <- 0
  dt <- 1
  
  # Initialization
  particles_sigma_sq <- rep(sigma0_sq, n_particles)
  
  # Main Filtering Loop
  for (t in 1:n_obs) {
    # --- Stage 1: Auxiliary Step (Look-Ahead) ---
    
    # 1a. Make a simplified prediction of the state at time t, based on particles at t-1.
    # This prediction ignores the random jumps to create a tractable "first-guess".
    predicted_sigma_sq <- pmax(1e-9, particles_sigma_sq * (1 - lambda * dt))
    # Note: The mean of the return also ignores the jump term for this first stage.
    predicted_means <- (mu - 0.5 * predicted_sigma_sq) * dt
    predicted_sds <- sqrt(predicted_sigma_sq * dt)
    
    # 1b. Calculate first-stage weights using the CURRENT observation data[t].
    log_first_stage_weights <- dnorm(data[t], mean = predicted_means, sd = predicted_sds, log = TRUE)
    
    # Normalize
    max_log_weight <- max(log_first_stage_weights, na.rm = TRUE)
    if (!is.finite(max_log_weight)) return(1e9)
    first_stage_weights <- exp(log_first_stage_weights - max_log_weight)
    
    # 1c. Resample ancestor indices from t-1 based on these "look-ahead" weights.
    if (sum(first_stage_weights) == 0 || !is.finite(sum(first_stage_weights))) return(1e9)
    ancestor_indices <- sample(1:n_particles, size = n_particles, replace = TRUE, prob = first_stage_weights)
    
    # --- Stage 2: Propagation and Correction ---
    
    # 2a. Propagate ONLY the chosen ancestors using the FULL model dynamics.
    propagated_particles <- particles_sigma_sq[ancestor_indices]
    
    # Simulate the BDLP increment for the chosen particles.
    bdlp_increments <- simulate_ig_ou_bdlp_increment(n_particles, lambda * dt, a_ig, b_ig)
    
    new_particles_sigma_sq <- pmax(1e-9, propagated_particles * (1 - lambda * dt) + bdlp_increments)
    
    # 2b. Calculate second-stage (correction) weights.
    # CRITICAL FIX: The mean now includes the jump term rho * dz.
    new_means <- (mu - 0.5 * new_particles_sigma_sq) * dt + rho * bdlp_increments
    new_sds <- sqrt(new_particles_sigma_sq * dt)
    log_numerator_weights <- dnorm(data[t], mean = new_means, sd = new_sds, log = TRUE)
    
    # Denominator is the likelihood from the simplified prediction for the CHOSEN ancestor.
    log_denominator_weights <- dnorm(data[t], 
                                     mean = predicted_means[ancestor_indices], 
                                     sd = predicted_sds[ancestor_indices], 
                                     log = TRUE)
    
    log_second_stage_weights <- log_numerator_weights - log_denominator_weights
    
    # Normalize
    max_log_weight2 <- max(log_second_stage_weights, na.rm = TRUE)
    if (!is.finite(max_log_weight2)) return(1e9)
    second_stage_weights <- exp(log_second_stage_weights - max_log_weight2)
    
    # 2c. Update total log-likelihood.
    if (mean(second_stage_weights) == 0 || !is.finite(mean(second_stage_weights))) return(1e9)
    log_likelihood <- log_likelihood + 
      log(mean(first_stage_weights)) + 
      max_log_weight + max_log_weight2 + 
      log(mean(second_stage_weights))
    
    # 2d. Final resampling for the next iteration.
    final_normalized_weights <- second_stage_weights / sum(second_stage_weights)
    final_indices <- sample(1:n_particles, size = n_particles, replace = TRUE, prob = final_normalized_weights)
    
    particles_sigma_sq <- new_particles_sigma_sq[final_indices]
  }
  
  return(-log_likelihood)
}

# --- 3. Define and Run the Optimization ---
mle_objective_function_apf_ig <- function(scaled_params, data, n_particles) {
  params <- numeric(6)
  params[1] <- scaled_params[1]
  params[2] <- -exp(scaled_params[2])
  params[3] <- exp(scaled_params[3])
  params[4] <- exp(scaled_params[4]) # a_ig
  params[5] <- exp(scaled_params[5]) # b_ig
  params[6] <- exp(scaled_params[6])
  
  neg_log_lik <- apf_log_likelihood_ig(params, data, n_particles)
  
  cat("Testing Params (IG):", round(params, 6), " | -LogLik (APF):", round(neg_log_lik, 6), "\n")
  
  if (!is.finite(neg_log_lik)) return(1e9)
  return(neg_log_lik)
}

# Initial parameters (plausible guesses for IG-OU)
initial_scaled_params_ig <- c(
  mu = mean(log_returns),
  rho_trans = log(0.9), 
  log_lambda = log(0.01),
  log_a_ig = log(1.0),
  log_b_ig = log(0.05),
  log_sigma0_sq = log(var(log_returns))
)

n_particles <- 1000

cat("\nStarting MLE optimization for BNS-IG via APF (this will be slow)...\n")
mle_results_apf_ig <- optim(
  par = initial_scaled_params_ig,
  fn = mle_objective_function_apf_ig,
  data = log_returns,
  n_particles = n_particles,
  method = "Nelder-Mead",
  control = list(maxit = 200, trace = 1) 
)

# --- 4. Display Results ---
cat("\n--- APF MLE Optimization Finished ---\n")
print("Optimal Scaled Parameters Found:")
print(mle_results_apf_ig$par)

# Transform parameters back to their original scale
estimated_params_mle_apf_ig <- numeric(6)
estimated_params_mle_apf_ig[1] <- mle_results_apf_ig$par[1]
estimated_params_mle_apf_ig[2] <- -exp(mle_results_apf_ig$par[2])
estimated_params_mle_apf_ig[3] <- exp(mle_results_apf_ig$par[3])
estimated_params_mle_apf_ig[4] <- exp(mle_results_apf_ig$par[4])
estimated_params_mle_apf_ig[5] <- exp(mle_results_apf_ig$par[5])
estimated_params_mle_apf_ig[6] <- exp(mle_results_apf_ig$par[6])
names(estimated_params_mle_apf_ig) <- c("mu", "rho", "lambda", "a", "b", "sigma0_sq")

print("Optimal Interpretable Parameters (MLE):")
print(estimated_params_mle_apf_ig)

cat("\nFinal Minimized Negative Log-Likelihood:", mle_results_apf_ig$value, "\n")

cat("\nSimulating final BNS model path with estimated parameters...\n")

bns_returns <- simulate_bns_ig_sv(
  params = estimated_params_mle_apf_ig,
  n_steps = length(log_returns),
  dt = 1
)
cat("BNS simulation complete.\n")


# --- Comparison Plots ---
cat("Generating comparison plots...\n")

# A. Density Overlay Plot
x_grid <- seq(min(log_returns), max(log_returns), length.out = 500)
gh_density_data <- data.frame(x = x_grid, y = dghyp(x_grid, object = gh_fit))

density_plot_hist_style <- ggplot(data.frame(value = log_returns), aes(x = value)) +
  # Empirical data as a histogram.
  geom_histogram(aes(y = after_stat(density), fill = "S&P 500 Empirical"), 
                 binwidth = 0.001, 
                 alpha = 0.7, 
                 color = "white") +
  # BNS simulated density as a line
  geom_density(data = data.frame(value = bns_returns), 
               aes(color = "BNS (Simulated)"), linewidth = 1.2) +
  # GH fitted density as a line
  geom_line(data = gh_density_data, 
            aes(x = x, y = y, color = "GH (Fitted)"), 
            linewidth = 1.2) +
  # Labels and Titles
  labs(
    title = "Static Comparison: Unconditional Return Distributions", 
    x = "Log Return", 
    y = "Density"
  ) +
  # Colors and Theme
  scale_fill_manual(values = c("S&P 500 Empirical" = "lightblue")) +
  scale_color_manual(values = c("BNS (Simulated)" = "blue", "GH (Fitted)" = "red")) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = quantile(log_returns, c(0.005, 0.995))) +
  guides(fill = guide_legend(title = NULL), color = guide_legend(title = "Model Overlays"))

print(density_plot_hist_style)

# B. Q-Q Plot Comparison
# Generate random sample from the fitted GH model
gh_samples <- rghyp(length(log_returns), object = gh_fit)

# Prepare data for Q-Q plots
qq_data <- data.frame(
  empirical = sort(log_returns), 
  gh = sort(gh_samples), 
  bns = sort(bns_returns))

qq_plot_gh <- ggplot(qq_data, aes(x = gh, y = empirical)) +
  geom_point(alpha = 0.5, color = "red") +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(
    title = "Q-Q Plot: Empirical vs. Fitted GH", 
    x = "Theoretical Quantiles (GH)", 
    y = "Empirical Quantiles (S&P 500)"
  ) +
  theme_minimal(base_size = 14)

qq_plot_bns <- ggplot(qq_data, aes(x = bns, y = empirical)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(
    title = "Q-Q Plot: Empirical vs. Simulated BNS", 
    x = "Theoretical Quantiles (BNS)", 
    y = "Empirical Quantiles (S&P 500)"
  ) +
  theme_minimal(base_size = 14)

# Arrange plots side-by-side
QQplots <- grid.arrange(qq_plot_gh, qq_plot_bns, ncol = 2)
print(QQplots)

# C. Dynamic Comparison (Volatility Clustering)
# Calculate ACF for Absolute returns for all three series
acf_empirical <- acf(abs(log_returns), plot = FALSE, lag.max = 50)
acf_bns <- acf(abs(bns_returns), plot = FALSE, lag.max = 50)
acf_gh <- acf(abs(gh_samples), plot = FALSE, lag.max = 50)

# Create a data frame for the ACF plot
acf_data <- data.frame(
  Lag = acf_empirical$lag,
  ACF = c(acf_empirical$acf, acf_bns$acf, acf_gh$acf),
  Model = factor(rep(c("S&P 500 Empirical", "BNS (Simulated)", "GH (Static)"), 
                     each = length(acf_empirical$lag)))
)

# Confidence interval for ACF plot (for white noise)
ci <- qnorm(0.975) / sqrt(length(log_returns))

acf_plot <- ggplot(acf_data, aes(x = Lag, y = ACF, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = c(ci, -ci), linetype = "dashed", color = "black") +
  labs(
    title = "Dynamic Comparison: Autocorrelation of Absolute Returns", 
    subtitle = "Demonstrating Volatility Clustering", 
    x = "Lag (Days)", 
    y = "Autocorrelation"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("S&P 500 Empirical" = "black", 
                               "BNS (Simulated)" = "blue", 
                               "GH (Static)" = "red")) +
  theme(legend.position = "bottom")

print(acf_plot)

## Save plots
ggsave(
  filename = here("outputs", "BNS-IG(MLE-APF)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-IG(MLE-APF)&GH_QQplot.png"),
  plot = QQplots,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-IG(MLE-APF)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

