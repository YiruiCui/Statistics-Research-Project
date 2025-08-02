library(here)  # Load {here} package for file path management

# Identify project location
here::i_am("analysis/SP500_SF_check.R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))