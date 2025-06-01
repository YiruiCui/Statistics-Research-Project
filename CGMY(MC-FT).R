# ----------------------------------------------
# CGMY simulation via Fourier inversion (Algorithm A from Ballotta & Kyriakou 2014)
# ----------------------------------------------

library(pracma)

# --- CGMY Parameters (Set I in paper) ---
C <- 0.5
G <- 2.0
M <- 3.5
Y <- 0.5
Delta <- 1      # time increment

# --- Numerical parameters ---
N <- 2^10       # number of grid points (FFT resolution)
L <- 189         # Fourier domain truncation
D <- 13.9          # support of real space [-D/2, D/2]

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
for (l in 1:N) {
  sum = 0
  for (j in 1:N) {
    sum = sum + exp(-1i * (j-1) * (l-1) * h * delta) * exp(-1i * u[j] * x[1]) * phi_reg[j]
  }
  F_d[l] <- (h/(2*pi)) * exp(-1i * u[1] * (x[l] - x[1])) * sum
}

# Compute CDF approximations
F_cdf <- Re(F_d) + 0.5

# --- Draw samples via inverse transform sampling ---
set.seed(7914)
n_samples <- 100000
u_samples <- runif(n_samples)
x_samples <- approx(F_cdf, x, xout = u_samples, rule = 2)$y

# --- Plot ---
hist(x_samples, breaks = 200, freq = FALSE, col = "lightblue",
     main = "CGMY Increments via Algorithm A (MC-FT1)",
     xlab = "x")
lines(x, diff(c(0, F_cdf)) / delta, col = "red", lwd = 2)
legend("topright", legend = c("Simulated", "Numerical PDF"), col = c("lightblue", "red"), lwd = 2)
