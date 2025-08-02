# Required package for higher moments
if (!require(e1071)) install.packages("e1071")
library(e1071)
library(here)  # Load {here} package for file path management

# Identify project location
here::i_am("analysis/model_S&P_500_CGMY.R")

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

# --- Moment matching for CGMY parameters ---
estimate_CGMY <- function(m, s2, skew, kurt) {
  loss_fn <- function(par) {
    C <- par[1]
    G <- par[2]
    M <- par[3]
    Y <- par[4]
    #if (C <= 0 || G <= 0 || M <= 0 || Y <= 0) return(1e6)
    
    mu_th <- C*(M^(Y-1)-G^(Y-1))*gamma(1-Y)
    var_th <- C*(M^(Y-2)+G^(Y-2))*gamma(2-Y)
    skew_th <- C*(M^(Y-3)-G^(Y-3))*gamma(3-Y)/(var_th^(3/2))
    kurt_th <- 3 + C*(M^(Y-4)-G^(Y-4))*gamma(4-Y)/var_th^2
    
    sum((c(mu_th, var_th, skew_th, kurt_th) - c(m, s2, skew, kurt))^2)
  }
  
  # Initial guesses
  start <- c(0.0000001, 2, 1, 1.997)
  result <- optim(start, loss_fn, method = "L-BFGS-B",
                  lower = c(1e-6, 1e-3, 1e-3, 0.01),
                  upper = c(100, 10, 10, 1.99))
  
  names(result$par) <- c("C", "G", "M", "Y")
  return(result$par)
}

# Run estimation
CGMY_par <- estimate_CGMY(mu_emp, sigma2_emp, skew_emp, kurt_emp)
print(CGMY_par)

C <- 0.0000001  #CGMY_par[1]
G <- 2.07543122 #CGMY_par[2]
M <- 1.06335469 #CGMY_par[3]
Y <- 1.9972  #CGMY_par[4]
Delta <- 1      # time increment

# --- Numerical parameters ---
N <- 2^12       # number of grid points (FFT resolution)
L <- 555        # Fourier domain truncation
D <- 30.1          # support of real space [-D/2, D/2]

h <- L / N                      # spacing in Fourier domain
delta <- D / N                 # spacing in real domain
u <- (-N/2):(N/2 - 1) * h       # Fourier grid
x <- (-N/2):(N/2 - 1) * delta   # real grid

# --- CGMY characteristic function φ(u) ---
cgmy_cf <- function(u, C, G, M, Y, t = Delta) {
  psi <- C * gamma(-Y) * ((M - 1i * u)^Y - M^Y + (G + 1i * u)^Y - G^Y)
  return(exp(t * psi))
}

phi_u <- cgmy_cf(u, C, G, M, Y)

# --- Regularization φ^R(u) = (1 - cos(uL)) / (iu) * φ(u) ---
phi_reg <- rep(0, length(u))
nonzero_u <- u != 0
phi_reg[nonzero_u] <- -(1 - cos(u[nonzero_u] * D)) * phi_u[nonzero_u] / (1i * u[nonzero_u])
phi_reg[!nonzero_u] <- 0

# Numerical evaluation of CDF using FFT
F_d = 0
for (l in 1:N) {
  sum = 0
  for (j in 1:N) {
    sum = sum + exp(-1i * (j-1) * (l-1) * h * delta) * exp(-1i * u[j] * x[1]) * phi_reg[j]
  }
  F_d[l] <- (h/(2*pi)) * exp(-1i * u[1] * (x[l] - x[1])) * sum
}

# Compute CDF approximations
F_cdf <- Re(F_d) + 0.5

# --- Step 5: Plot comparison ---
hist(log_returns, breaks = 100, probability = TRUE, col = "lightblue",
     main = "S&P 500 Log-Returns vs CGMY Fit", xlab = "Log-Return")

# Overlay the CGMY PDF
lines(x, diff(c(0, F_cdf)) / delta, col = "red", lwd = 2)

legend("topright", legend = c("Empirical", "NIG Fit"), col = c("lightblue", "red"), lwd = 2, cex = 0.5)
