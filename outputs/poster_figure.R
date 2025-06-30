library(grid)
library(png)
library(gridExtra)
library(here)
# Identify project location
here::i_am("outputs/Untitled.R")

# Load PNGs
img1 <- readPNG(here("outputs","GH_fit01.png"))
img2 <- readPNG(here("outputs","GH_QQplot.png"))
img3 <- readPNG(here("outputs","Meixner_fit01.png"))
img4 <- readPNG(here("outputs","Meixner_QQplot.png"))

# Convert to grobs
g1 <- rasterGrob(img1, interpolate = TRUE)
g2 <- rasterGrob(img2, interpolate = TRUE)
g3 <- rasterGrob(img3, interpolate = TRUE)
g4 <- rasterGrob(img4, interpolate = TRUE)

png(filename = here("outputs", "combined_plots.png"), width = 4000, height = 2400, res = 300)
# Set 2 rows × 2 columns layout
par(mfrow = c(2, 2))

# ---- Plot 1: Histogram with GH fit ----
hist(log_returns, breaks = 150, probability = TRUE,
     col = "lightblue", main = "S&P 500 Log-Returns vs GH Fit", xlab = "Log-Return")
curve(dghyp(x, gh_fit), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "GH Fit"), col = c("lightblue", "red"), lwd = 2, cex = 0.8)

# ---- Plot 2: Q–Q Plot for GH ----
plot(q_theoretical, log_returns_sorted,
     main = "Q–Q Plot: GH Fit vs Empirical Returns",
     xlab = "Theoretical Quantiles (GH)",
     ylab = "Empirical Quantiles",
     pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)

# ---- Plot 3: Histogram with Meixner fit ----
hist(log_returns, breaks = 150, probability = TRUE,
     col = "lightblue", main = "S&P 500 Log-Returns vs Meixner Fit", xlab = "Log-Return")
lines(x_vals, pdf_vals, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "Meixner Fit"), col = c("lightblue", "red"), lwd = 2, cex = 0.8)

# ---- Plot 4: Q–Q Plot for Meixner ----
plot(q_theoretical_Meixner, log_returns_sorted,
     main = "Q–Q Plot: Meixner Fit vs Empirical Returns",
     xlab = "Theoretical Quantiles (Meixner)",
     ylab = "Empirical Quantiles",
     pch = 16, col = "darkblue", cex = 0.6)
abline(0, 1, col = "red", lwd = 2)
# Arrange in 2x2 layout
dev.off()
# Optional: reset layout
par(mfrow = c(1, 1))


