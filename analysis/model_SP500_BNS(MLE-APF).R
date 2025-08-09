library(here)
library(tidyverse)
library(ghyp)

set.seed(7914)

# Identify project location
here::i_am("analysis/model_SP500_BNS(MLE-APF).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))


# --- 2. Auxiliary Particle Filter Log-Likelihood Function ---

auxiliary_particle_filter_log_likelihood <- function(params, data, n_particles) {
  # Unpack parameters
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  a         <- params[4]
  b         <- params[5]
  sigma0_sq <- params[6]
  
  n_obs <- length(data)
  log_likelihood <- 0
  
  # --- Initialization ---
  particles_sigma_sq <- rep(sigma0_sq, n_particles)
  
  # --- Main Loop (Filtering) ---
  for (t in 1:n_obs) {
    # --- Stage 1: Auxiliary Step ---
    
    # 1a. Predict the auxiliary variable: E[sigma_t^2 | sigma_{t-1}^2]
    # For the Gamma-OU process, this is the mean-reverting part.
    # We ignore the jump part here as it's an approximation.
    predicted_sigma_sq <- pmax(1e-9, particles_sigma_sq * (1 - lambda))
    
    # 1b. Calculate first-stage weights based on the observation at time t.
    # These weights tell us which of the current particles are most likely
    # to produce a good prediction for the next step.
    predicted_means <- mu - 0.5 * predicted_sigma_sq
    predicted_sds <- sqrt(predicted_sigma_sq)
    
    log_first_stage_weights <- dnorm(data[t], mean = predicted_means, sd = predicted_sds, log = TRUE)
    
    # Normalize weights to avoid numerical issues
    max_log_weight <- max(log_first_stage_weights)
    first_stage_weights <- exp(log_first_stage_weights - max_log_weight)
    
    # 1c. First-stage resampling.
    # We sample the *indices* of the particles that will be propagated.
    # This focuses our computational effort on promising particles.
    ancestor_indices <- sample(1:n_particles, 
                               size = n_particles, 
                               replace = TRUE, 
                               prob = first_stage_weights)
    
    # --- Stage 2: Propagation and Final Weighting ---
    
    # 2a. Propagate the chosen ancestor particles forward one step.
    propagated_particles <- particles_sigma_sq[ancestor_indices]
    
    jump_intensity <- a * lambda * 1 # dt = 1
    jumps <- rpois(n_particles, jump_intensity)
    jump_sizes <- sapply(jumps, function(nj) {
      if (nj == 0) return(0)
      sum(rgamma(nj, shape = 1, rate = b))
    })
    
    # This is our new set of particles for time t
    new_particles_sigma_sq <- pmax(1e-9, propagated_particles * (1 - lambda) + jump_sizes)
    
    # 2b. Calculate second-stage weights (importance correction).
    # These weights correct for the approximation made in Stage 1.
    new_means <- mu - 0.5 * new_particles_sigma_sq
    new_sds <- sqrt(new_particles_sigma_sq)
    
    log_numerator_weights <- dnorm(data[t], mean = new_means, sd = new_sds, log = TRUE)
    
    # The denominator is the weight from the particle we chose in the first stage
    log_denominator_weights <- dnorm(data[t], 
                                     mean = predicted_means[ancestor_indices], 
                                     sd = predicted_sds[ancestor_indices], 
                                     log = TRUE)
    
    log_second_stage_weights <- log_numerator_weights - log_denominator_weights
    
    # Normalize
    max_log_weight2 <- max(log_second_stage_weights)
    second_stage_weights <- exp(log_second_stage_weights - max_log_weight2)
    
    # 2c. Update the total log-likelihood.
    # This combines the likelihood from both stages.
    log_likelihood <- log_likelihood + 
      log(mean(first_stage_weights)) + 
      max_log_weight + max_log_weight2 + 
      log(mean(second_stage_weights))
    
    # 2d. Final resampling.
    final_normalized_weights <- second_stage_weights / sum(second_stage_weights)
    final_indices <- sample(1:n_particles, 
                            size = n_particles, 
                            replace = TRUE, 
                            prob = final_normalized_weights)
    
    # Update the particle cloud for the next iteration
    particles_sigma_sq <- new_particles_sigma_sq[final_indices]
  }
  
  return(-log_likelihood) # Return negative log-likelihood for minimization
}

# --- 3. Define and Run the Optimization ---
# (The objective function and optimization call would be the same as in the
# previous script, but would call `auxiliary_particle_filter_log_likelihood`
# instead of the standard one).

# We can define the objective function
mle_objective_function_apf <- function(scaled_params, data, n_particles) {
  params <- numeric(6)
  params[1] <- scaled_params[1]
  params[2] <- scaled_params[2]
  params[3] <- exp(scaled_params[3])
  params[4] <- exp(scaled_params[4])
  params[5] <- exp(scaled_params[5])
  params[6] <- exp(scaled_params[6])
  
  neg_log_lik <- auxiliary_particle_filter_log_likelihood(params, data, n_particles)
  
  cat("Testing Params:", round(params, 4), " | -LogLik (APF):", round(neg_log_lik, 4), "\n")
  
  if (!is.finite(neg_log_lik)) return(1e9)
  return(neg_log_lik)
}

# Initial parameters (same as before)
initial_scaled_params <- c(
  mu = mean(log_returns),
  rho_trans = -0.5, 
  log_lambda = log(0.05),
  log_a = log(0.1),
  log_b = log(10),
  log_sigma0_sq = log(var(log_returns))
)

n_particles <- 1000

cat("\nStarting MLE optimization via Auxiliary Particle Filter (this will be very slow)...\n")
mle_results_apf_gamma <- optim(
  par = initial_scaled_params,
  fn = mle_objective_function_apf,
  data = log_returns,
  n_particles = n_particles,
  method = "Nelder-Mead",
  control = list(maxit = 200, trace = 1)
)

# --- 4. Display Results ---
cat("\n--- APF MLE Optimization Finished ---\n")
print("Optimal Scaled Parameters Found:")
print(mle_results_apf_gamma$par)

# Transform parameters back to their original scale
estimated_params_mle_apf_gamma <- numeric(6)
estimated_params_mle_apf_gamma[1] <- mle_results_apf_gamma$par[1]
estimated_params_mle_apf_gamma[2] <- mle_results_apf_gamma$par[2]
estimated_params_mle_apf_gamma[3] <- exp(mle_results_apf_gamma$par[3])
estimated_params_mle_apf_gamma[4] <- exp(mle_results_apf_gamma$par[4])
estimated_params_mle_apf_gamma[5] <- exp(mle_results_apf_gamma$par[5])
estimated_params_mle_apf_gamma[6] <- exp(mle_results_apf_gamma$par[6])
names(estimated_params_mle_apf_gamma) <- c("mu", "rho", "lambda", "a", "b", "sigma0_sq")

print("Optimal Interpretable Parameters (MLE):")
print(estimated_params_mle_apf_gamma)

cat("\nFinal Minimized Negative Log-Likelihood:", mle_results_apf_gamma$value, "\n")

cat("\nSimulating final BNS model path with estimated parameters...\n")

bns_returns <- simulate_bns_gamma_sv(
  params = estimated_params_mle_apf_gamma,
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
  filename = here("outputs", "BNS-gamma(MLE-APF)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-gamma(MLE-APF)&GH_QQplot.png"),
  plot = QQplots,
  width = 3000,
  height = 1500,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-gamma(MLE-APF)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

