# ----------------------------------------------
# Estimate BNS-IG parameters from S&P 500 data using Particle Filter based MLE
# ----------------------------------------------

# NOTE: Particle filtering is computationally intensive. (usually take ~8h to run)
# For a quick test, change line 147 n_particles to 100 and line 155 maxit to 20 (~1min run, poor result)

# --- Import required package ---
library(here)
library(tidyverse)
library(ghyp)
library(gridExtra)
library(SuppDists)

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_SP500_BNS-IG(MLE-APF).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Fit Generalized Hyperbolic distribution for plot comparison ---
gh_fit <- fit.ghypuv(log_returns, lambda = -0.5, symmetric = FALSE)

# --- Import BNS-IG Simulation Function ---
source(here("src","BNS-IG.R"))

# --- Auxiliary Particle Filter (APF) for BNS with IG-SV ---
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
    predicted_sigma_sq <- pmax(1e-9, particles_sigma_sq * (1 - lambda * dt))
    predicted_means <- (mu - 0.5 * predicted_sigma_sq) * dt
    predicted_sds <- sqrt(predicted_sigma_sq * dt)
    
    # 1b. Calculate first-stage weights using the current observation data[t].
    # These weights determine which particles at t-1 are good "ancestors".
    log_first_stage_weights <- dnorm(data[t], mean = predicted_means, sd = predicted_sds, log = TRUE)
    
    # Normalize weights
    max_log_weight <- max(log_first_stage_weights, na.rm = TRUE)
    if (!is.finite(max_log_weight)) return(1e9)
    first_stage_weights <- exp(log_first_stage_weights - max_log_weight)
    
    # 1c. Resample the ancestor indices from t-1 based on the first-stage weights.
    if (sum(first_stage_weights) == 0 || !is.finite(sum(first_stage_weights))) return(1e9)
    ancestor_indices <- sample(1:n_particles, size = n_particles, replace = TRUE, prob = first_stage_weights)
    
    # --- Stage 2: Propagation and Correction ---
    
    # 2a. Propagate only the chosen ancestors using the full model dynamics.
    propagated_particles <- particles_sigma_sq[ancestor_indices]
    
    # Simulate the BDLP increment for the chosen particles.
    bdlp_increments <- simulate_ig_ou_bdlp_increment(n_particles, lambda * dt, a_ig, b_ig)
    
    # Generate new cloud of particles at time t.
    new_particles_sigma_sq <- pmax(1e-9, propagated_particles * (1 - lambda * dt) + bdlp_increments)
    
    # 2b. Calculate second-stage (correction) weights.
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

# --- Define the MLE Objective Function ---
mle_objective_function_apf_ig <- function(scaled_params, data, n_particles) {
  params <- numeric(6)
  params[1] <- scaled_params[1]
  params[2] <- -exp(scaled_params[2])
  params[3] <- exp(scaled_params[3])
  params[4] <- exp(scaled_params[4]) 
  params[5] <- exp(scaled_params[5]) 
  params[6] <- exp(scaled_params[6])
  
  neg_log_lik <- apf_log_likelihood_ig(params, data, n_particles)
  
  cat("Testing Params (IG):", round(params, 6), " | -LogLik (APF):", round(neg_log_lik, 6), "\n")
  
  if (!is.finite(neg_log_lik)) return(1e9)
  return(neg_log_lik)
}

# --- Run the Optimization ---
# Starting values on the transformed scale
initial_scaled_params_ig <- c(
  mu = mean(log_returns),
  rho_trans = log(0.5), 
  log_lambda = log(0.05),
  log_a_ig = log(0.05),
  log_b_ig = log(100),
  log_sigma0_sq = log(var(log_returns))
)

n_particles <- 10000

mle_results_apf_ig <- optim(
  par = initial_scaled_params_ig,
  fn = mle_objective_function_apf_ig,
  data = log_returns,
  n_particles = n_particles,
  method = "Nelder-Mead",
  control = list(maxit = 200, trace = 1) 
)

# Transform parameters back to their original scale
estimated_params_mle_apf_ig <- numeric(6)
estimated_params_mle_apf_ig[1] <- mle_results_apf_ig$par[1]
estimated_params_mle_apf_ig[2] <- -exp(mle_results_apf_ig$par[2])
estimated_params_mle_apf_ig[3] <- exp(mle_results_apf_ig$par[3])
estimated_params_mle_apf_ig[4] <- exp(mle_results_apf_ig$par[4])
estimated_params_mle_apf_ig[5] <- exp(mle_results_apf_ig$par[5])
estimated_params_mle_apf_ig[6] <- exp(mle_results_apf_ig$par[6])
names(estimated_params_mle_apf_ig) <- c("mu", "rho", "lambda", "a", "b", "sigma0_sq")

# --- Display Results ---
print("Optimal Interpretable Parameters (MLE):")
print(estimated_params_mle_apf_ig)

cat("\nFinal Minimized Negative Log-Likelihood:", mle_results_apf_ig$value, "\n")

# Simulate BNS process for comparison
set.seed(7914)
bns_returns <- simulate_bns_ig_sv(
  params = estimated_params_mle_apf_ig,
  n_steps = length(log_returns),
  dt = 1
)


# --- Comparison Plots ---

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
               aes(color = "BNS (IG)"), linewidth = 1.2) +
  # GH fitted density as a line
  geom_line(data = gh_density_data, 
            aes(x = x, y = y, color = "GH"), 
            linewidth = 1.2) +
  # Labels and Titles
  labs(
    title = "Unconditional Return Distributions", 
    x = "Log Return", 
    y = "Density"
  ) +
  # Colors and Theme
  scale_fill_manual(values = c("S&P 500 Empirical" = "lightblue")) +
  scale_color_manual(values = c("BNS (IG)" = "blue", "GH" = "red")) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = quantile(log_returns, c(0.005, 0.995))) +
  guides(fill = guide_legend(title = NULL), color = guide_legend(title = "Model Overlays"))

print(density_plot_hist_style)

# B. Q-Q Plot Comparison
# Generate samples from the fitted GH model
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
    title = "Q-Q Plot: Empirical vs. Fitted BNS (IG)", 
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
  Model = factor(rep(c("S&P 500 Empirical", "BNS (IG)", "GH"), 
                     each = length(acf_empirical$lag)))
)

# Confidence interval for ACF plot (for white noise)
ci <- qnorm(0.975) / sqrt(length(log_returns))

acf_plot <- ggplot(acf_data, aes(x = Lag, y = ACF, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = c(ci, -ci), linetype = "dashed", color = "black") +
  labs(
    title = "Autocorrelation of Absolute Returns", 
    subtitle = "Demonstrating Volatility Clustering", 
    x = "Lag (Days)", 
    y = "Autocorrelation"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("S&P 500 Empirical" = "black", 
                               "BNS (IG)" = "blue", 
                               "GH" = "red")) +
  theme(legend.position = "bottom")

print(acf_plot)

## Save plots
ggsave(
  filename = here("outputs", "BNS-IG(MLE-APF-10000p)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-IG(MLE-APF-10000p)&GH_QQplot.png"),
  plot = QQplots,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-IG(MLE-APF-10000p)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

