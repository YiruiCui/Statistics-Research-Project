library(here)
library(tidyverse)
library(ghyp)
library(stabledist)
library(gridExtra)
library(moments)
library(SuppDists) # For rinvGauss function

set.seed(7914)

# Identify project location
here::i_am("analysis/model_SP500_BNS-TS(MLE-APF).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Function to Simulate Tempered Stable Variates ---
rts <- function(n, dt, kappa, a, b, K = 1000) {
  
  # Use sapply to generate n independent values
  sapply(1:n, function(i) {
    # For each increment, we need a fresh set of random numbers
    e <- rexp(K)
    u_tilde <- runif(K)
    b_i <- cumsum(rexp(K, rate = 1))
    
    # Calculate the jump sizes. Here, T from the formula is our time step 'dt'.
    # Since we want the value at the end of the interval, the indicator is always 1.
    jumps <- pmin(
      2 * (a * dt / (b_i * gamma(1 - kappa)))^(1 / kappa),
      (2 * e * u_tilde^(1 / kappa) / b^(1 / kappa))
    )
    
    # The value of the increment is the sum of all K potential jumps
    sum(jumps)
  })
}


# --- BNS Simulation Function with TS Stochastic Volatility ---
simulate_bns_ts_sv <- function(params, n_steps, dt) {
  # Unpack parameters for the BNS model with TS-SV
  # TS-OU is parameterized by (kappa, a, b) and lambda for mean reversion.
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3] # Mean-reversion speed for volatility
  kappa     <- params[4] # TS parameter kappa (0 < kappa < 1)
  a         <- params[5] # TS parameter 'a'
  b         <- params[6] # TS parameter 'b'
  sigma0_sq <- params[7] # Initial variance
  
  # Initialize vectors for the simulation path
  log_S <- numeric(n_steps + 1)
  sigma_sq <- numeric(n_steps + 1)
  
  # Set initial values
  log_S[1] <- 0
  sigma_sq[1] <- sigma0_sq
  
  # --- Simulate the Background Driving Lévy Process (BDLP) for TS-OU ---
  # According to Schoutens (2003), Section 5.5.3, the BDLP z is the sum
  # of two independent Lévy processes: z = z^(1) + z^(2).
  
  # Component 1: A TS(kappa, kappa*a, b) Lévy process, z^(1)
  # An increment over dt*lambda has the law of a TS(kappa, kappa*a*lambda*dt, b) variable.
  # We use our new rts function to simulate these increments.
  z1_increments <- rts(n_steps, dt=1, kappa = kappa, a = kappa * a * lambda * dt, b = b)
  
  # --- Component 2: A compound Poisson process z^(2) ---
  # Jump intensity is a*b*kappa
  jump_intensity_z2 <- a * b * kappa * (lambda * dt)
  jumps_z2 <- rpois(n_steps, jump_intensity_z2)
  
  # Jump sizes are Gamma(1 - kappa, b^(1/kappa) / 2)
  jump_sizes_z2 <- sapply(jumps_z2, function(nj) {
    if (nj == 0) return(0)
    sum(rgamma(nj, shape = 1 - kappa, rate = b^(1/kappa) / 2))
  })
  
  # The total increment of the BDLP at each step
  bdlp_increments <- z1_increments + jump_sizes_z2
  
  # Generate standard normal random variables for the main Brownian motion
  dW <- rnorm(n_steps, mean = 0, sd = sqrt(dt))
  
  # --- Euler-Maruyama Discretization for the BNS System ---
  for (t in 1:n_steps) {
    sigma_sq_prev <- max(1e-9, sigma_sq[t]) # Ensure variance is non-negative
    
    # Update variance process (TS-OU)
    # d(sigma_t^2) = -lambda * sigma_t^2 * dt + d(z_{lambda*t})
    d_sigma_sq <- -lambda * sigma_sq_prev * dt + bdlp_increments[t]
    sigma_sq[t+1] <- sigma_sq_prev + d_sigma_sq
    
    # Update log-return process (BNS)
    # d(log S_t) = (mu - 0.5*sigma_t^2)dt + sigma_t*dW_t + rho*d(z_{lambda*t})
    d_log_S <- (mu - 0.5 * sigma_sq_prev) * dt + 
      sqrt(sigma_sq_prev) * dW[t] + 
      rho * bdlp_increments[t]
    
    log_S[t+1] <- log_S[t] + d_log_S
  }
  
  # Return the simulated log-returns
  return(diff(log_S))
}

# --- 2. Auxiliary Particle Filter for BNS with TS-SV ---

apf_log_likelihood_ts <- function(params, data, n_particles) {
  # Unpack parameters for BNS with TS-SV (7 parameters)
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  kappa     <- params[4]
  a         <- params[5]
  b         <- params[6]
  sigma0_sq <- params[7]
  
  n_obs <- length(data)
  log_likelihood <- 0
  dt <- 1 # Assuming daily data
  
  # Initialization
  particles_sigma_sq <- rep(sigma0_sq, n_particles)
  
  # Main Filtering Loop
  for (t in 1:n_obs) {
    # --- Stage 1: Auxiliary Step ---
    predicted_sigma_sq <- pmax(1e-9, particles_sigma_sq * (1 - lambda * dt))
    predicted_means <- mu - 0.5 * predicted_sigma_sq
    predicted_sds <- sqrt(predicted_sigma_sq)
    
    log_first_stage_weights <- dnorm(data[t], mean = predicted_means, sd = predicted_sds, log = TRUE)
    max_log_weight <- max(log_first_stage_weights)
    first_stage_weights <- exp(log_first_stage_weights - max_log_weight)
    
    ancestor_indices <- sample(1:n_particles, 
                               size = n_particles, 
                               replace = TRUE, 
                               prob = first_stage_weights)
    
    # --- Stage 2: Propagation and Final Weighting ---
    propagated_particles <- particles_sigma_sq[ancestor_indices]
    
    # --- Simulate the BDLP for TS-OU ---
    # Component 1: A TS(kappa, kappa*a, b) Lévy process, z^(1)
    z1_increments <- rts(n_particles, dt=1, kappa = kappa, a = kappa * a * lambda * dt, b = b)
    
    # Component 2: A compound Poisson process z^(2)
    jump_intensity_z2 <- a * b * kappa * (lambda * dt)
    jumps_z2 <- rpois(n_particles, jump_intensity_z2)
    jump_sizes_z2 <- sapply(jumps_z2, function(nj) {
      if (nj == 0) return(0)
      sum(rgamma(nj, shape = 1 - kappa, rate = b^(1/kappa) / 2))
    })
    
    bdlp_increments <- z1_increments + jump_sizes_z2
    
    # Update variance process (TS-OU)
    new_particles_sigma_sq <- pmax(1e-9, propagated_particles * (1 - lambda * dt) + bdlp_increments)
    
    # Calculate second-stage weights (importance correction)
    new_means <- mu - 0.5 * new_particles_sigma_sq
    new_sds <- sqrt(new_particles_sigma_sq)
    
    log_numerator_weights <- dnorm(data[t], mean = new_means, sd = new_sds, log = TRUE)
    log_denominator_weights <- dnorm(data[t], 
                                     mean = predicted_means[ancestor_indices], 
                                     sd = predicted_sds[ancestor_indices], log = TRUE)
    
    log_second_stage_weights <- log_numerator_weights - log_denominator_weights
    max_log_weight2 <- max(log_second_stage_weights)
    second_stage_weights <- exp(log_second_stage_weights - max_log_weight2)
    
    # Update total log-likelihood
    log_likelihood <- log_likelihood + 
      log(mean(first_stage_weights)) + 
      max_log_weight + max_log_weight2 + 
      log(mean(second_stage_weights))
    
    # Final resampling
    final_normalized_weights <- second_stage_weights / sum(second_stage_weights)
    if(any(is.na(final_normalized_weights)) || sum(final_normalized_weights) == 0) {
      final_indices <- sample(1:n_particles, size = n_particles, replace = TRUE)
    } else {
      final_indices <- sample(1:n_particles, 
                              size = n_particles, 
                              replace = TRUE, 
                              prob = final_normalized_weights)
    }
    
    particles_sigma_sq <- new_particles_sigma_sq[final_indices]
  }
  
  return(-log_likelihood)
}

# --- 3. Define and Run the Optimization ---
mle_objective_function_apf_ts <- function(scaled_params, data, n_particles) {
  params <- numeric(7)
  params[1] <- scaled_params[1]
  params[2] <- scaled_params[2]
  params[3] <- exp(scaled_params[3]) # lambda
  params[4] <- scaled_params[4]      # kappa is in (0,1)
  params[5] <- exp(scaled_params[5]) # a
  params[6] <- exp(scaled_params[6]) # b
  params[7] <- exp(scaled_params[7]) # sigma0_sq
  
  neg_log_lik <- apf_log_likelihood_ts(params, data, n_particles)
  
  cat("Testing Params (TS):", round(params, 4), " | -LogLik (APF):", round(neg_log_lik, 4), "\n")
  
  if (!is.finite(neg_log_lik)) return(1e9)
  return(neg_log_lik)
}

# Initial parameters (plausible guesses for TS-OU)
initial_scaled_params_ts <- c(
  mu = mean(log_returns),
  rho_trans = -0.7, 
  log_lambda = log(0.05),
  kappa = 0.7, # Bounded between 0 and 1
  log_a = log(1.5),
  log_b = log(10.0),
  log_sigma0_sq = log(var(log_returns))
)

n_particles <- 100

cat("\nStarting MLE optimization for BNS-TS via APF (this will be very slow)...\n")
# Using optim with bounds for kappa
mle_results_apf_ts <- optim(
  par = initial_scaled_params_ts,
  fn = mle_objective_function_apf_ts,
  data = log_returns,
  n_particles = n_particles,
  method = "L-BFGS-B", # Use a method that supports bounds
  lower = c(-Inf, -Inf, -Inf, 1e-6, -Inf, -Inf, -Inf),
  upper = c(Inf, Inf, Inf, 1-1e-6, Inf, Inf, Inf),
  control = list(maxit = 200, trace = 1) 
)

# --- 4. Display Results ---
cat("\n--- APF MLE Optimization Finished ---\n")
print("Optimal Scaled Parameters Found:")
print(mle_results_apf_ts$par)

# Transform parameters back to their original scale
estimated_params_mle_apf_ts <- numeric(6)
estimated_params_mle_apf_ts[1] <- mle_results_apf_ts$par[1]
estimated_params_mle_apf_ts[2] <- mle_results_apf_ts$par[2]
estimated_params_mle_apf_ts[3] <- exp(mle_results_apf_ts$par[3])
estimated_params_mle_apf_ts[4] <- mle_results_apf_ts$par[4]
estimated_params_mle_apf_ts[5] <- exp(mle_results_apf_ts$par[5])
estimated_params_mle_apf_ts[6] <- exp(mle_results_apf_ts$par[6])
estimated_params_mle_apf_ts[7] <- exp(mle_results_apf_ts$par[7])
names(estimated_params_mle_apf_ts) <- c("mu", "rho", "lambda", "kappa", "a", "b", "sigma0_sq")

print("Optimal Interpretable Parameters (MLE):")
print(estimated_params_mle_apf_ts)

cat("\nFinal Minimized Negative Log-Likelihood:", mle_results_apf_ts$value, "\n")

cat("\nSimulating final BNS model path with estimated parameters...\n")

bns_returns <- simulate_bns_ts_sv(
  params = estimated_params_mle_apf_ts,
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
  filename = here("outputs", "BNS-TS(MLE-APF)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-TS(MLE-APF)&GH_QQplot.png"),
  plot = QQplots,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-TS(MLE-APF)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

