# Required package for higher moments
if (!require(e1071)) install.packages("e1071")
library(e1071)
library(here)  # Load {here} package for file path management

# Identify project location
here::i_am("data/model_S&P_500_ghyp.R")

# --- Step 1: Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Step 2: Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Empirical moments ---
mu_emp <- mean(log_returns)
sigma2_emp <- var(log_returns)
skew_emp <- skewness(log_returns)
kurt_emp <- kurtosis(log_returns) + 3  # e1071 returns excess kurtosis

# Display
cat("Estimated Moments:\n")
cat("Mean     :", mu_emp, "\n")
cat("Variance :", sigma2_emp, "\n")
cat("Skewness :", skew_emp, "\n")
cat("Kurtosis :", kurt_emp, "\n")

# --- Meixner PDF using complex gamma from pracma ---
meixner_pdf <- function(x, alpha, beta, delta) {
  z <- complex(real = delta, imaginary = x / alpha)
  gamma_vals <- pracma::gammaz(z) 
  factor <- ((2 * cos(beta / 2))^(2 * delta)) / (2 * alpha * pi * gamma(2 * delta))
  density <- factor * exp(beta * x / alpha) * (Mod(gamma_vals)^2)
  return(density)
}

# --- Moment matching for Meixner parameters ---
estimate_meixner <- function(m, s2, skew, kurt) {
  loss_fn <- function(par) {
    alpha <- par[1]
    beta <- par[2]
    delta <- par[3]
    if (delta <= 0 || alpha <= 0 || abs(beta) >= pi) return(1e6)
    
    mu_th <- alpha * delta * tan(beta / 2)
    var_th <- 0.5 * alpha^2 * delta / cos(beta / 2)^2
    skew_th <- sin(beta / 2) * sqrt(2 / delta)
    kurt_th <- 3 + (2 - cos(beta)) / delta
    
    sum((c(mu_th, var_th, skew_th, kurt_th) - c(m, s2, skew, kurt))^2)
  }
  
  # Initial guesses
  start <- c(0.0559, -0.15065, 0.09638)
  result <- optim(start, loss_fn, method = "L-BFGS-B",
                  lower = c(1e-4, -pi + 1e-4, 1e-4),
                  upper = c(100, pi - 1e-4, 100))
  
  names(result$par) <- c("alpha", "beta", "delta")
  return(result$par)
}

# Run estimation
meixner_par <- estimate_meixner(mu_emp, sigma2_emp, skew_emp, kurt_emp)
print(meixner_par)

# --- Estimate Meixner parameters ---
delta <- meixner_par[3] # (2-0.9886723)/(kurt_emp-3)  
beta <- meixner_par[2] # 2*asin(skew_emp/sqrt(2/delta))  
alpha <- meixner_par[1] # sqrt(2*sigma2_emp / delta)*cos(beta/2)#0.06 

# --- Step 5: Plot comparison ---
hist(log_returns, breaks = 200, probability = TRUE, col = "lightblue",
     main = "S&P 500 Log-Returns vs Meixner Fit", xlab = "Log-Return")

# Compute Meixner PDF over the range of histogram
x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
pdf_vals <- meixner_pdf(x_vals, alpha, beta, delta)

# Overlay the Meixner PDF
lines(x_vals, pdf_vals, col = "red", lwd = 2)

legend("topright", legend = c("Empirical", "Meixner Fit"), col = c("lightblue", "red"), lwd = 2, cex = 0.5)

