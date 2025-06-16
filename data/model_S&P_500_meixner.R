library(here)  # Load {here} package for file path management

# Identify project location
here::i_am("data/model_S&P_500_ghyp.R")

# --- Step 1: Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[23112:nrow(data)]

# --- Step 2: Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- Empirical moments ---
mu_emp <- mean(log_returns)
sigma2_emp <- var(log_returns)
sigma_emp <- sqrt(sigma2_emp)
mu3 <- mean((log_returns - mu_emp)^3)
mu4 <- mean((log_returns - mu_emp)^4)
skew_emp <- mu3 / sigma2_emp^(3/2)
kurt_emp <- mu4 / sigma2_emp^2

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

# --- Meixner PDF with location parameter m ---
meixner_pdf_m <- function(x, m, a, b, d) {
  z <- (x - m) / a
  density <- meixner_pdf(z, alpha = 1, beta = b, delta = d)
  return(density / a)
}

# --- Estimate Meixner parameters ---
delta <- 1/ (kurt_emp - skew_emp^2 - 3)  
beta <- sign(skew_emp) * acos(2 - delta * (kurt_emp - 3))  
alpha <- sigma_emp * sqrt((cos(beta)+1)/delta)
m <- mu_emp - alpha * delta * tan(beta/2)

# --- Step 5: Plot comparison ---
hist(log_returns, breaks = 100, probability = TRUE, col = "lightblue",
     main = "S&P 500 Log-Returns vs Meixner Fit", xlab = "Log-Return")

# Compute Meixner PDF over the range of histogram
x_vals <- seq(min(log_returns), max(log_returns), length.out = 1000)
pdf_vals <- meixner_pdf_m(x_vals, m, alpha, beta, delta)

# Overlay the Meixner PDF
lines(x_vals, pdf_vals, col = "red", lwd = 2)

legend("topright", legend = c("Empirical", "Meixner Fit"), col = c("lightblue", "red"), lwd = 2, cex = 0.5)

