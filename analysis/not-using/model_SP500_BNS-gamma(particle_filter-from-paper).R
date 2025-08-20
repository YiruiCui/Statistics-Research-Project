# This script implements the exact particle filter for the Barndorff-Nielsen & Shephard (BNS)
# model with Gamma-OU stochastic volatility, with the weighting formula corrected to
# match the one specified in Section 5.4.4 of the 2001 paper.

# --- 1. Load Libraries ---
library(here)
library(tidyverse)
library(gridExtra)

# Set seed for reproducibility
set.seed(7914)

# Identify project location
here::i_am("analysis/model_SP500_BNS-gamma(particle_filter-from-paper).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- 2. Helper Functions (Corrected and Finalized) ---

simulate_gamma_ou_innovation <- function(nu, alpha, lambda, t, K = 1) {
  poisson_rate <- nu * lambda * t
  if (poisson_rate <= 0) {
    return(matrix(0, nrow = K, ncol = 2, dimnames = list(NULL, c("eta1", "eta2"))))
  }
  
  innovations <- lapply(seq_len(K), function(k) {
    n_to_sim <- ceiling(poisson_rate + 10 * sqrt(poisson_rate) + 5)
    inter_arrival_times <- rexp(n_to_sim, rate = poisson_rate)
    arrival_times <- cumsum(inter_arrival_times)
    c_i <- arrival_times[arrival_times <= 1]
    num_jumps <- length(c_i)
    
    if (num_jumps == 0) return(c(eta1 = 0, eta2 = 0))
    
    r_i <- runif(num_jumps)
    sum_for_eta1 <- sum(-log(c_i) * exp(lambda * t * r_i))
    sum_for_eta2 <- sum(-log(c_i))
    
    eta1 <- (1 / alpha) * exp(-lambda * t) * sum_for_eta1
    eta2 <- (1 / alpha) * sum_for_eta2
    
    return(c(eta1 = eta1, eta2 = eta2)) 
  })
  
  return(do.call(rbind, innovations))
}

update_bns_state <- function(previous_state, eta_n, lambda, delta_t) {
  new_sigma_sq <- exp(-lambda * delta_t) * previous_state[1] + eta_n[1]
  new_z <- previous_state[2] + eta_n[2]
  return(c(sigma_sq = new_sigma_sq, z = new_z))
}

calculate_integrated_vol <- function(current_state, previous_state, lambda) {
  delta_z <- current_state[2] - previous_state[2]
  delta_sigma_sq <- current_state[1] - previous_state[1]
  integrated_vol <- (1 / lambda) * (delta_z - delta_sigma_sq)
  return(integrated_vol)
}

# --- 3. The Particle Filter for a Model WITH LEVERAGE ---

run_bns_particle_filter_leverage <- function(data, params, M, K, delta_t = 1) {
  n_obs <- length(data)
  
  # Initialization
  particles_state <- data.frame(
    sigma_sq = rep(params$sigma0_sq, times=M),
    z = 0
  )
  
  log_likelihood <- 0
  
  # --- Main Loop ---
  for (t in 1:n_obs) {
    offspring_states <- data.frame(
      sigma_sq = numeric(M * K), z = numeric(M * K), 
      integrated_vol = numeric(M * K), eta2 = numeric(M*K) 
    )
    
    for (m in 1:M) {
      previous_state_m <- c(
        sigma_sq = particles_state$sigma_sq[m],
        z = particles_state$z[m]
      )
      
      innovations_mk <- simulate_gamma_ou_innovation(
        params$nu, params$alpha, params$lambda, delta_t, K = K
      )
      for (k in 1:K) {
        idx <- (m - 1) * K + k
        eta_n <- innovations_mk[k, ]
        current_state_mk <- update_bns_state(previous_state_m, eta_n, params$lambda, delta_t)
        integrated_vol_mk <- calculate_integrated_vol(current_state_mk, previous_state_m, params$lambda)
        
        offspring_states$sigma_sq[idx] <- current_state_mk[1]
        offspring_states$z[idx] <- current_state_mk[2]
        offspring_states$integrated_vol[idx] <- pmax(1e-9, integrated_vol_mk)
        offspring_states$eta2[idx] <- eta_n[2]
      }
    }
    
    expected_eta2 <- (params$nu / params$alpha) * params$lambda * delta_t
    leverage_effect <- params$rho * offspring_states$eta2
    particle_means <- params$mu * delta_t - 0.5 * offspring_states$integrated_vol + leverage_effect

    log_liks <- dnorm(data[t], 
                         mean = particle_means, 
                         sd = sqrt(offspring_states$integrated_vol), 
                         log = TRUE)
    
    log_weights <- -0.5 * log(offspring_states$integrated_vol) - 
      (data[t]^2) / (2 * offspring_states$integrated_vol)
    
    weights_unnormalized <- exp(log_weights)
    normalized_weights <- weights_unnormalized / sum(weights_unnormalized, na.rm = TRUE)
    
    log_likelihood <- log_likelihood + sum(normalized_weights * log_liks)
    
    if(any(is.na(normalized_weights)) || sum(normalized_weights) == 0) {
      indices <- sample(1:(M * K), size = M, replace = TRUE)
    } else {
      indices <- sample(1:(M * K), size = M, replace = TRUE, prob = normalized_weights)
    }
    
    particles_state <- offspring_states[indices, c("sigma_sq", "z")]
  }
  return(-log_likelihood)
}

# --- 4. Corrected MLE Objective Function ---
mle_objective_function_leverage <- function(scaled_params, data, M, K) {
  
  # Unscale parameters to their natural, constrained domain
  params <- list(
    mu     = scaled_params[1],
    nu     = exp(scaled_params[2]),      # nu > 0
    alpha  = exp(scaled_params[3]),      # alpha > 0
    lambda = exp(scaled_params[4]),      # lambda > 0
    rho    = -exp(scaled_params[5]),     # rho <= 0
    sigma0_sq = exp(scaled_params[6])
  )
  
  # Calculate the negative log-likelihood
  neg_log_lik <- run_bns_particle_filter_leverage(data, params, M, K)
  
  cat("Params (mu, nu, alpha, lambda, rho, sigma0_sq):", round(unlist(params), 6), 
      " | -LogLik:", round(neg_log_lik, 4), "\n")
  
  if (!is.finite(neg_log_lik)) {
    return(1e9) 
  }
  
  return(neg_log_lik)
}

# --- 5. Corrected Optimization Run ---
initial_scaled_params_leverage <- c(
  mu = mean(log_returns),
  log_nu = log(1.0),
  log_alpha = log(10),
  log_lambda = log(0.5),
  log_abs_rho = log(1.0),
  log_sigma0_sq = log(var(log_returns))
)

M_particles <- 10
K_offspring <- 3

cat("\nStarting MLE optimization with leverage (this may be slow)...\n")
mle_results <- optim(
  par = initial_scaled_params_leverage,
  fn = mle_objective_function_leverage,
  data = log_returns,
  M = M_particles,
  K = K_offspring,
  method = "Nelder-Mead",
  control = list(maxit = 100) # Use a small maxit for this example
)

# --- 5. Display Results ---
cat("\n--- MLE Optimization Finished ---\n")
print("Optimal Scaled Parameters Found:")
print(mle_results$par)

# Transform parameters back to their original scale
estimated_params_pf_gamma <- list(
  mu     = mle_results$par[1],
  nu     = exp(mle_results$par[2]),
  alpha  = exp(mle_results$par[3]),
  lambda = exp(mle_results$par[4]),
  rho    = exp(mle_results$par[5]),
  sigma0_sq = exp(mle_results$par[6])
)

print("Optimal Interpretable Parameters (MLE):")
print(estimated_params_pf_gamma)

cat("\nFinal Minimized Negative Log-Likelihood:", mle_results$value, "\n")

# --- Main BNS Simulation Function ---
simulate_bns_gamma_sv <- function(params, S0, n_steps, delta_t) {
  
  # --- 1. Initialization ---
  # Create a data frame to store the results
  path <- tibble(
    time = seq(0, n_steps, by = delta_t),
    log_price = 0.0,
    asset_price = 0.0,
    inst_volatility = 0.0,
    z_process = 0.0
  )
  
  # Set initial values at time 0
  path$log_price[1] <- log(S0)
  path$asset_price[1] <- S0
  # Draw initial volatility from its stationary Gamma distribution
  path$inst_volatility[1] <- params$sigma0_sq
  path$z_process[1] <- 0
  
  # Calculate the constant expected value of the jump component (eta2)
  expected_eta2 <- (params$nu / params$alpha) * params$lambda * delta_t
  
  # --- 2. Simulation Loop ---
  for (i in 2:(n_steps + 1)) {
    # Get state from the previous time step
    prev_sigma_sq <- path$inst_volatility[i - 1]
    prev_z <- path$z_process[i - 1]
    prev_log_price <- path$log_price[i - 1]
    
    # a) Simulate the random innovation for the interval
    eta_n <- simulate_gamma_ou_innovation(
      params$nu, params$alpha, params$lambda, delta_t
    )
    
    # b) Update the hidden volatility states (based on Eq. 49)
    current_sigma_sq <- exp(-params$lambda * delta_t) * prev_sigma_sq + eta_n[1]
    current_z <- prev_z + eta_n[2]
    
    # c) Update the log price (based on the discrete version of Eq. 8)
    drift_term <- (params$mu - 1/2 * prev_sigma_sq) * delta_t
    diffusion_term <- sqrt(prev_sigma_sq * delta_t) * rnorm(1)
    leverage_term <- params$rho * eta_n[2]
    
    current_log_price <- prev_log_price + drift_term + diffusion_term + leverage_term
    
    # Store the new values
    path$log_price[i] <- current_log_price
    path$asset_price[i] <- exp(current_log_price)
    path$inst_volatility[i] <- current_sigma_sq
    path$z_process[i] <- current_z
  }
  
  return(path)
}

cat("\nSimulating final BNS model path with estimated parameters...\n")

bns_process <- simulate_bns_gamma_sv(
  params = estimated_params_pf_gamma,
  S0 = data$Last_Price[1],
  n_steps = length(log_returns),
  delta_t = 1
)

bns_returns <- diff(bns_process$log_price)
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
  filename = here("outputs", "BNS-gamma(PF-paper)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-gamma(PF-paper)&GH_QQplot.png"),
  plot = QQplots,
  width = 3000,
  height = 1500,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-gamma(PF-paper)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)