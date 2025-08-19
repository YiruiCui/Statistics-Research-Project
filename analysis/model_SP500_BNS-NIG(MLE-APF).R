# --- Load Libraries ---
library(here)
library(tidyverse)
library(ghyp)
library(fBasics) # For rNIG function
library(gridExtra)
library(moments)

# Set a seed for reproducibility
set.seed(7914)

# --- Load and Prepare Data ---
# This assumes your project is set up with a 'data' folder
# and your script is in an 'analysis' folder.
here::i_am("analysis/model_SP500_BNS-NIG(MLE-APF).R")
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]
log_returns <- diff(log(price_vector))

# --- Helper Functions for COS-based Simulation (z3) ---
cf_z3 <- function(omega, dt, alpha, beta, delta) {
  rho_param <- beta / alpha
  u <- 1i * omega
  term1 <- beta * sqrt((alpha - beta) / (alpha + beta))
  term2 <- (u + beta) * sqrt((alpha - u - beta) / (alpha + u + beta))
  exponent <- dt * rho_param * delta * (term1 - term2)
  return(exp(exponent))
}

r_from_cf_cos <- function(n, cf_func, a_int, b_int, N = 2^14) {
  cos_density <- function(x) {
    k <- 0:(N - 1)
    u <- k * pi / (b_int - a_int)
    F_k <- Re(cf_func(u) * exp(-1i * u * a_int))
    F_k[1] <- F_k[1] / 2
    x_scaled <- (x - a_int) / (b_int - a_int)
    cos_matrix <- cos(outer(x_scaled, k * pi))
    density <- (2 / (b_int - a_int)) * (cos_matrix %*% F_k)
    return(pmax(0, as.vector(density)))
  }
  cos_cdf <- function(x) {
    sapply(x, function(val) {
      if (val <= a_int) return(0)
      if (val >= b_int) return(1)
      integrate(cos_density, lower = a_int, upper = val)$value
    })
  }
  cos_inverse_cdf <- function(p) {
    sapply(p, function(val) {
      if (val <= 0) return(a_int)
      if (val >= 1) return(b_int)
      uniroot(function(x) cos_cdf(x) - val, interval = c(a_int, b_int))$root
    })
  }
  uniform_samples <- runif(n)
  return(cos_inverse_cdf(uniform_samples))
}

# --- Helper Function to Simulate the NIG-OU BDLP Increment ---
# Implements the decomposition z = z^(1) + z^(2) + z^(3) from Schoutens (2003).
# The z^(3) component requires advanced simulation (like the COS method) and is
# computational intensive, can omit(default)/include control by parameter full
simulate_nig_ou_bdlp_increment <- function(n, dt, alpha, beta, delta, full=FALSE) {
  
  rho_param <- beta / alpha
  
  # Component 1: An NIG-Lévy process z^(1)
  z1_delta <- (1 - rho_param) * delta
  z1_increments <- fBasics::rnig(n, alpha = alpha, beta = beta, delta = z1_delta * dt, mu = 0)
  
  # Component 2: A compound Poisson process z^(2)
  poisson_intensity <- (1 / (delta * alpha * sqrt((1 - rho_param) / (1 + rho_param)))) * dt
  num_jumps <- rpois(n, lambda = poisson_intensity)
  
  z2_increments <- sapply(num_jumps, function(nj) {
    if (nj == 0) return(0)
    v_sq <- rnorm(nj)^2
    v_tilde_sq <- rnorm(nj)^2
    return(sum(v_sq - v_tilde_sq) / (alpha * sqrt(1 - rho_param^2)))
  })
  
  if (full){
    # Component 3: z^(3)
    cf_z3_param <- function(omega) cf_z3(omega, dt, alpha, beta, delta)
    integration_range <- c(-0.5, 0.5) 
    z3_increments <- r_from_cf_cos(n, cf_func = cf_z3_param, 
                                   a_int = integration_range[1], 
                                   b_int = integration_range[2])
    
    return(z1_increments + z2_increments + z3_increments)
  }
  else{
    return(z1_increments + z2_increments)
  }
  
}


# --- BNS Simulation Function with NIG Stochastic Volatility ---
simulate_bns_nig_sv <- function(params, n_steps, dt, full=FALSE) {
  # Unpack parameters
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  alpha_nig <- params[4]
  beta_nig  <- params[5]
  delta_nig <- params[6]
  sigma0_sq <- params[7]
  
  # Initialize vectors
  log_S <- numeric(n_steps + 1); log_S[1] <- 0
  sigma_sq <- numeric(n_steps + 1); sigma_sq[1] <- sigma0_sq
  
  # Simulate increments for the driving processes
  dW <- rnorm(n_steps, mean = 0, sd = sqrt(dt))
  bdlp_increments <- simulate_nig_ou_bdlp_increment(
    n = n_steps, 
    dt = lambda * dt,
    alpha = alpha_nig, 
    beta = beta_nig, 
    delta = delta_nig, 
    full = full
  )
  
  # Euler-Maruyama Discretization
  for (t in 1:n_steps) {
    sigma_sq_prev <- max(1e-9, sigma_sq[t])
    
    # Update variance (NIG-OU)
    sigma_sq[t+1] <- sigma_sq[t] - lambda * sigma_sq_prev * dt + bdlp_increments[t]
    
    # Update log-price (BNS)
    log_S[t+1] <- log_S[t] + (mu - 0.5 * sigma_sq_prev) * dt + 
      sqrt(sigma_sq_prev) * dW[t] + 
      rho * bdlp_increments[t]
  }
  
  return(diff(log_S))
}

# --- Auxiliary Particle Filter for BNS with NIG-SV ---
apf_log_likelihood_nig <- function(params, data, n_particles, full=FALSE) {
  # Unpack parameters
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  alpha_nig <- params[4]
  beta_nig  <- params[5]
  delta_nig <- params[6]
  sigma0_sq <- params[7]
  
  n_obs <- length(data)
  log_likelihood <- 0
  dt <- 1 # Daily data
  
  # Initialization
  particles_sigma_sq <- rep(sigma0_sq, n_particles)
  
  # Main Filtering Loop
  for (t in 1:n_obs) {
    # Stage 1: Auxiliary Step
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
    
    # Stage 2: Propagation and Weighting
    propagated_particles <- particles_sigma_sq[ancestor_indices]
    
    # Simulate the BDLP for NIG-OU for each particle
    bdlp_increments <- simulate_nig_ou_bdlp_increment(
      n = n_particles, 
      dt = lambda * dt,
      alpha = alpha_nig, 
      beta = beta_nig, 
      delta = delta_nig, 
      full = full
    )
    
    # Update variance particles
    new_particles_sigma_sq <- pmax(1e-9, propagated_particles * (1 - lambda * dt) + bdlp_increments)
    
    # Calculate second-stage weights
    new_means <- mu - 0.5 * new_particles_sigma_sq
    new_sds <- sqrt(new_particles_sigma_sq)
    log_numerator_weights <- dnorm(data[t], mean = new_means, sd = new_sds, log = TRUE)
    log_denominator_weights <- dnorm(data[t], 
                                     mean = predicted_means[ancestor_indices], 
                                     sd = predicted_sds[ancestor_indices], 
                                     log = TRUE)
    
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

# --- Define and Run the Optimization ---
mle_objective_function_apf_nig <- function(scaled_params, data, n_particles) {
  params <- numeric(7)
  params[1] <- scaled_params[1]          # mu
  params[2] <- -exp(scaled_params[2])    # rho
  params[3] <- exp(scaled_params[3])     # lambda
  params[4] <- exp(scaled_params[4])     # alpha_nig
  params[5] <- scaled_params[5]          # beta_nig
  params[6] <- exp(scaled_params[6])     # delta_nig
  params[7] <- exp(scaled_params[7])     # sigma0_sq
  
  # Constraint: |beta| < alpha
  # Add a small buffer (e.g., 1e-6) to prevent floating point issues at the boundary.
  if (params[4] - abs(params[5]) < 1e-6) {
    return(1e9) # Return a large value if constraint is violated
  }
  
  neg_log_lik <- apf_log_likelihood_nig(params, data, n_particles)
  
  cat("Testing Params (NIG):", round(params, 4), " | -LogLik (APF):", round(neg_log_lik, 4), "\n")
  
  if (!is.finite(neg_log_lik)) return(1e9)
  return(neg_log_lik)
}

# Initial parameters (plausible guesses for NIG-OU)
initial_scaled_params_nig <- c(
  mu = mean(log_returns),
  rho = log(0.7), 
  log_lambda = log(0.01),
  log_alpha_nig = log(10.0),
  beta_nig = -1.0,
  log_delta_nig = log(0.005),
  log_sigma0_sq = log(var(log_returns))
)

n_particles <- 1000 # Increase for accuracy, decrease for speed

cat("\nStarting MLE optimization for BNS-NIG via APF (this will be very slow)...\n")
mle_results_apf_nig <- optim(
  par = initial_scaled_params_nig,
  fn = mle_objective_function_apf_nig,
  data = log_returns,
  n_particles = n_particles,
  method = "Nelder-Mead", # Nelder-Mead is robust for complex likelihood surfaces
  control = list(maxit = 200, trace = 1)
)

# --- Display Results ---
cat("\n--- APF MLE Optimization for BNS-NIG Finished ---\n")
print("Optimal Scaled Parameters Found:")
print(mle_results_apf_nig$par)

# Transform parameters back to their original scale
estimated_params_mle_apf_nig <- numeric(7)
estimated_params_mle_apf_nig[1] <- mle_results_apf_nig$par[1]
estimated_params_mle_apf_nig[2] <- -exp(mle_results_apf_nig$par[2])
estimated_params_mle_apf_nig[3] <- exp(mle_results_apf_nig$par[3])
estimated_params_mle_apf_nig[4] <- exp(mle_results_apf_nig$par[4])
estimated_params_mle_apf_nig[5] <- mle_results_apf_nig$par[5]
estimated_params_mle_apf_nig[6] <- exp(mle_results_apf_nig$par[6])
estimated_params_mle_apf_nig[7] <- exp(mle_results_apf_nig$par[7])
names(estimated_params_mle_apf_nig) <- c("mu", "rho", "lambda", "alpha", "beta", "delta", "sigma0_sq")

print("Optimal Interpretable Parameters (MLE for BNS-NIG):")
print(estimated_params_mle_apf_nig)

cat("\nFinal Minimized Negative Log-Likelihood:", mle_results_apf_nig$value, "\n")

cat("\nSimulating final BNS model path with estimated parameters...\n")

bns_returns <- simulate_bns_nig_sv(
  params = estimated_params_mle_apf_nig,
  n_steps = length(log_returns),
  dt = 1,
  full = FALSE
)
cat("BNS simulation complete.\n")


# --- Comparison Plots ---
cat("Generating comparison plots...\n")

# A. Density Overlay Plot
x_grid <- seq(min(log_returns), max(log_returns), length.out = 500)
gh_density_data <- data.frame(x = x_grid, y = ghyp::dghyp(x_grid, object = gh_fit))

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
  filename = here("outputs", "BNS-NIG(MLE-APF)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-NIG(MLE-APF)&GH_QQplot.png"),
  plot = QQplots,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-NIG(MLE-APF)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)