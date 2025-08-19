library(here)
library(tidyverse)
library(ghyp)
library(SuppDists)
library(gridExtra)
library(moments)

set.seed(7914)

# --- Load and Prepare Data ---
here::i_am("analysis/model_SP500_BNS-IG(MLE-PF).R")
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]
log_returns <- diff(log(price_vector))

# --- Fit GH model for comparison plots ---
gh_fit <- fit.ghypuv(log_returns, lambda = -0.5, symmetric = FALSE, silent = TRUE)


# --- 2. Standard Particle Filter (SIR) Log-Likelihood Function ---
particle_filter_log_likelihood_ig <- function(params, data, n_particles) {
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
    # --- Step 1: Prediction/Propagation (Propose) ---
    # Move all particles forward one step using the state equation (IG-OU process).
    
    # Simulate the BDLP for IG-OU
    z1_increments <- rinvGauss(n_particles, nu = (a_ig / 2) * (lambda*dt), lambda = (b_ig^2) * (a_ig / 2) * (lambda*dt))
    jump_intensity_z2 <- (a_ig * b_ig / 2) * (lambda*dt)
    jumps_z2 <- rpois(n_particles, jump_intensity_z2)
    jump_sizes_z2 <- sapply(jumps_z2, function(nj) {
      if (nj == 0) return(0)
      sum((1 / b_ig^2) * (rnorm(nj)^2))
    })
    bdlp_increments <- z1_increments + jump_sizes_z2
    
    particles_sigma_sq <- pmax(1e-9, particles_sigma_sq * (1 - lambda * dt) + bdlp_increments)
    
    # --- Step 2: Weighting (Correct) ---
    # Calculate the likelihood (weight) of the observation data[t] for each new particle position.
    particle_means <- mu - 0.5 * particles_sigma_sq
    particle_sds <- sqrt(particles_sigma_sq)
    
    log_weights <- dnorm(data[t], mean = particle_means, sd = particle_sds, log = TRUE)
    
    # Normalize weights to prevent numerical underflow
    max_log_weight <- max(log_weights)
    weights <- exp(log_weights - max_log_weight)
    
    # --- Step 3: Update Log-Likelihood ---
    log_likelihood <- log_likelihood + max_log_weight + log(mean(weights))
    
    # --- Step 4: Resampling ---
    normalized_weights <- weights / sum(weights)
    if(any(is.na(normalized_weights)) || sum(normalized_weights) == 0) {
      indices <- sample(1:n_particles, size = n_particles, replace = TRUE)
    } else {
      indices <- sample(1:n_particles, size = n_particles, replace = TRUE, prob = normalized_weights)
    }
    particles_sigma_sq <- particles_sigma_sq[indices]
  }
  
  return(-log_likelihood)
}

# --- 3. Define and Run the Optimization ---
mle_objective_function_pf_ig <- function(scaled_params, data, n_particles) {
  params <- numeric(6)
  params[1] <- scaled_params[1]
  params[2] <- -exp(scaled_params[2])
  params[3] <- exp(scaled_params[3])
  params[4] <- exp(scaled_params[4])
  params[5] <- exp(scaled_params[5])
  params[6] <- exp(scaled_params[6])
  
  neg_log_lik <- particle_filter_log_likelihood_ig(params, data, n_particles)
  
  cat("Testing Params (IG):", round(params, 4), " | -LogLik (PF):", round(neg_log_lik, 4), "\n")
  
  if (!is.finite(neg_log_lik)) return(1e9)
  return(neg_log_lik)
}

# Initial parameters
initial_scaled_params_ig <- c(
  mu = mean(log_returns),
  rho_trans = log(0.7), 
  log_lambda = log(0.05),
  log_a_ig = log(1e-5),
  log_b_ig = log(1.0),
  log_sigma0_sq = log(var(log_returns))
)

n_particles <- 1000

cat("\nStarting MLE optimization for BNS-IG via Standard Particle Filter...\n")
mle_results_pf_ig <- optim(
  par = initial_scaled_params_ig,
  fn = mle_objective_function_pf_ig,
  data = log_returns,
  n_particles = n_particles,
  method = "Nelder-Mead",
  control = list(maxit = 200, trace = 1) 
)

# --- 4. Display Results ---
cat("\n--- APF MLE Optimization Finished ---\n")
print("Optimal Scaled Parameters Found:")
print(mle_results_pf_ig$par)

# Transform parameters back to their original scale
estimated_params_mle_pf_ig <- numeric(6)
estimated_params_mle_pf_ig[1] <- mle_results_pf_ig$par[1]
estimated_params_mle_pf_ig[2] <- -exp(mle_results_pf_ig$par[2])
estimated_params_mle_pf_ig[3] <- exp(mle_results_pf_ig$par[3])
estimated_params_mle_pf_ig[4] <- exp(mle_results_pf_ig$par[4])
estimated_params_mle_pf_ig[5] <- exp(mle_results_pf_ig$par[5])
estimated_params_mle_pf_ig[6] <- exp(mle_results_pf_ig$par[6])
names(estimated_params_mle_pf_ig) <- c("mu", "rho", "lambda", "a", "b", "sigma0_sq")

print("Optimal Interpretable Parameters (MLE):")
print(estimated_params_mle_pf_ig)

cat("\nFinal Minimized Negative Log-Likelihood:", mle_results_pf_ig$value, "\n")

cat("\nSimulating final BNS model path with estimated parameters...\n")

bns_returns <- simulate_bns_ig_sv(
  params = estimated_params_mle_pf_ig,
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
  filename = here("outputs", "BNS-IG(MLE-PF)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-IG(MLE-PF)&GH_QQplot.png"),
  plot = QQplots,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-IG(MLE-PF)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

