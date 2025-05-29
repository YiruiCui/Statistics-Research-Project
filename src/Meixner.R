# ----------------------------------------------
# Simulate Meixner Process via AcceptReject package
# using book-based Meixner(α, β, δ) density
# ----------------------------------------------

require(AcceptReject)
require(pracma)

# --- Meixner PDF using complex gamma from pracma ---
meixner_pdf <- function(x, alpha, beta, delta) {
  z <- complex(real = delta, imaginary = x / alpha)
  gamma_vals <- pracma::gammaz(z) 
  factor <- ((2 * cos(beta / 2))^(2 * delta)) / (2 * alpha * pi * gamma(2 * delta))
  density <- factor * exp(beta * x / alpha) * (Mod(gamma_vals)^2)
  return(density)
}


# --- Parameters ---
alpha <- 1.7
beta <- 0.1
delta <- 0.9

# --- Accept-Reject Sampling with AcceptReject package ---
set.seed(7914)
meixner_samples <- AcceptReject::accept_reject(
  n = 10000L,
  f = meixner_pdf,
  continuous = TRUE,
  args_f = list(alpha=alpha, beta=beta, delta=delta),
  xlim = c(-1000, 1000)
)

# --- Plot 1: Histogram with True PDF Overlay ---
hist_data <- hist(meixner_samples, breaks = 100, probability = TRUE,
                  main = "Meixner Distribution: Samples + True PDF",
                  xlab = "x", col = "lightgray", border = "white")

# Compute true PDF over the range of histogram
x_vals <- seq(min(hist_data$breaks), max(hist_data$breaks), length.out = 1000)
pdf_vals <- meixner_pdf(x_vals, alpha, beta, delta)

# Overlay the true PDF
lines(x_vals, pdf_vals, col = "red", lwd = 2)

# Add legend
legend("topright", legend = c("True PDF"), col = "red", lwd = 2, bty = "n")


# --- Plot 2: Meixner Process Path ---
T <- 100
dt <- 0.1
n_steps <- length(seq(0, T, by = dt))
meix_increments <- meixner_samples[1:(n_steps - 1)]
meix_process <- c(0, cumsum(meix_increments))
time_grid <- seq(0, T, by = dt)

plot(time_grid, meix_process, type = "s", lwd = 2, col = "firebrick",
     main = "Simulated Meixner Process (via AcceptReject)",
     xlab = "Time", ylab = expression(X[t]))
