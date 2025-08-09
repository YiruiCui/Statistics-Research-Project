# --- Load Libraries and Data ---
library(here)
library(tidyverse)
library(moments)
library(ghyp)
library(gridExtra)

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_SP500_BNS(GMM).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Define the Moment Calculation Function ---
calculate_moments <- function(returns, num_acf_lags = 10) {
  m1 <- mean(returns)
  m2 <- var(returns)
  m3 <- skewness(returns)
  m4 <- kurtosis(returns)
  acf_abs_returns <- acf(abs(returns), lag.max = num_acf_lags, plot = FALSE)$acf
  return(c(mean = m1, variance = m2, skewness = m3, kurtosis = m4, acf_abs = acf_abs_returns))
}

# --- Calculate Empirical Moments from Real Data ---
cat("Calculating empirical moments from S&P 500 data...\n")
empirical_moments <- calculate_moments(log_returns, num_acf_lags = 10)
print("Empirical Moments:")
print(empirical_moments)


# --- BNS Simulation Function ---
simulate_bns_gamma_sv <- function(params, n_steps, dt) {
  # Unpack parameters from the input vector
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  a         <- params[4]
  b         <- params[5]
  sigma0_sq <- params[6]
  
  log_S <- numeric(n_steps + 1)
  sigma_sq <- numeric(n_steps + 1)
  log_S[1] <- 0
  sigma_sq[1] <- sigma0_sq
  
  jump_intensity <- a * lambda * dt
  jumps <- rpois(n_steps, jump_intensity)
  max_jumps <- max(1, jumps)
  jump_sizes <- matrix(0, nrow = n_steps, ncol = max_jumps)
  if (max_jumps > 0) {
    for(i in 1:n_steps) {
      if(jumps[i] > 0) {
        jump_sizes[i, 1:jumps[i]] <- rgamma(jumps[i], shape = 1, rate = b)
      }
    }
  }
  dW <- rnorm(n_steps, mean = 0, sd = sqrt(dt))
  
  for (t in 1:n_steps) {
    sigma_sq_prev <- max(1e-8, sigma_sq[t])
    d_sigma_sq <- -lambda * sigma_sq_prev * dt + sum(jump_sizes[t,])
    sigma_sq[t+1] <- sigma_sq_prev + d_sigma_sq
    d_log_S <- (mu - 0.5 * sigma_sq_prev) * dt + 
      sqrt(sigma_sq_prev) * dW[t] + 
      rho * sum(jump_sizes[t,])
    log_S[t+1] <- log_S[t] + d_log_S
  }
  return(diff(log_S))
}

# --- Define the GMM Objective Function ---
gmm_objective_function_scaled <- function(scaled_params, empirical_moments, n_sim_steps) {
  
  # --- Unscale Parameters (with Parameter Scaling) ---
  # The optimizer works with 'scaled_params'. We transform them back to their
  # natural scale before using them in the simulation.
  # This is a powerful technique to improve optimizer stability.
  params <- numeric(6)
  params[1] <- scaled_params[1]              # mu (no transform)
  params[2] <- tanh(scaled_params[2])         # rho is transformed to be in (-1, 1)
  params[3] <- exp(scaled_params[3])          # lambda = exp(log_lambda)
  params[4] <- exp(scaled_params[4])          # a = exp(log_a)
  params[5] <- exp(scaled_params[5])          # b = exp(log_b)
  params[6] <- exp(scaled_params[6])          # sigma0_sq = exp(log_sigma0_sq)
  
  # --- Simulation ---
  simulated_returns <- simulate_bns_gamma_sv(params, n_steps = n_sim_steps, dt = 1)
  
  # --- Moment Calculation ---
  simulated_moments <- calculate_moments(simulated_returns, 
                                         num_acf_lags = (length(empirical_moments) - 5))
  
  # --- Distance Calculation ---
  moment_differences <- simulated_moments - empirical_moments
  W <- diag(length(moment_differences))
  distance <- t(moment_differences) %*% W %*% moment_differences
  
  # Print progress
  cat("Testing Params:", round(params, 4), " | Distance:", round(distance, 4), "\n")
  
  # Handle non-finite results from simulation/moment calculation
  if (!is.finite(distance)) {
    return(1e9)
  }
  
  return(distance)
}

# --- Run the Optimization ---
# Starting values on the transformed scale
initial_scaled_params <- c(
  mu = mean(log_returns),
  rho_trans = atanh(-0.5), 
  log_lambda = log(0.05),
  log_a = log(0.1),
  log_b = log(10),
  log_sigma0_sq = log(var(log_returns))
)

n_simulation_steps <- 50000  #50000

cat("\nStarting GMM optimization...\n")
gmm_results <- optim(
  par = initial_scaled_params,
  fn = gmm_objective_function_scaled,
  # Additional arguments
  empirical_moments = empirical_moments,
  n_sim_steps = n_simulation_steps,
  method = "Nelder-Mead",
  control = list(maxit = 1000) 
)

# --- Display Results and Simulate Final Path ---
cat("\n--- GMM Optimization Finished ---\n")
print("Optimal Scaled Parameters Found:")
print(gmm_results$par)

# Transform parameters back to their original scale for interpretation
estimated_params <- numeric(6)
estimated_params[1] <- gmm_results$par[1]
estimated_params[2] <- tanh(gmm_results$par[2])
estimated_params[3] <- exp(gmm_results$par[3])
estimated_params[4] <- exp(gmm_results$par[4])
estimated_params[5] <- exp(gmm_results$par[5])
estimated_params[6] <- exp(gmm_results$par[6])
names(estimated_params) <- c("mu", "rho", "lambda", "a", "b", "sigma0_sq")

print("Optimal Interpretable Parameters:")
print(estimated_params)

cat("\nFinal Minimized Distance (Objective Value):", gmm_results$value, "\n")


cat("\nSimulating final BNS model path with estimated parameters...\n")

bns_returns <- simulate_bns_gamma_sv(
  params = estimated_params,
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
  filename = here("outputs", "BNS-gamma(GMM)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-gamma(GMM)&GH_QQplot.png"),
  plot = QQplots,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-gamma(GMM)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)
