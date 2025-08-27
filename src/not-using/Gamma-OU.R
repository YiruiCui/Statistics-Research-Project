# ----------------------------------------------
# Simulate Gamma–OU process using Schoutens (8.4.6)
# ----------------------------------------------

# --- Parameters ---
lambda <- 10       # mean reversion
a <- 10            # Poisson intensity parameter
b <- 100           # Exponential rate (i.e., Gamma(1, b))
y0 <- 0.08
T <- 1             # total time
dt <- 0.001        # fine resolution as in the book
time_grid <- seq(0, T, by = dt)
n_steps <- length(time_grid)

# --- Initialize ---
yt <- numeric(n_steps)
yt[1] <- y0
jump_times <- c()
jump_sizes <- c()

# --- Simulate Poisson process with intensity a * lambda ---
set.seed(7914)
increments <- rpois(n_steps - 1, a * lambda * dt)
poisson_process <- c(0, cumsum(increments))
N_total <- poisson_process[n_steps]
jump_sizes <- rexp(N_total, rate = b)

# --- Build process over time ---
jump_index <- 1
for (i in 2:n_steps) {

  # Decay from previous value
  yt[i] <- exp(-lambda * dt) * yt[i - 1]
  
  # Add new jumps in (t_prev, t_now]
  while (jump_index <= N_total && 
         jump_index <= poisson_process[i] && 
         jump_index > poisson_process[i-1]) {
    yt[i] <- yt[i] + jump_sizes[jump_index]
    jump_index <- jump_index + 1
  }
}

# --- Plot the Gamma–OU process path ---
plot(time_grid, yt, type = "l", col = "darkorange", lwd = 2,
     main = "Simulated Gamma–OU Process (from Schoutens)",
     xlab = "Time", ylab = expression(y[t]))
