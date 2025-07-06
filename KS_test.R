library(KSgeneral)
library(twosamples)
# --- Accept-Reject Sampling with AcceptReject package ---
set.seed(7914)
meixner_samples <- AcceptReject::accept_reject(
  n = 10000L,
  f = meixner_pdf_m,
  continuous = TRUE,
  args_f = list(m=m, a=alpha, b=beta, d=delta),
  xlim = c(-100, 100)
)

# Computes p-value of two-sided KS test: log-returns vs fitted Meixner samples
KS2sample(log_returns, meixner_samples, conservative = T)

# Two-sample Kuiper test
kuiper_mexiner = kuiper_test(log_returns, meixner_samples)
summary(kuiper_mexiner)
#Kuiper2sample(log_returns, meixner_samples, conservative = T)

gh_samples <- rghyp(10000, 
                    object = ghyp(lambda = gh_params[1], 
                                  alpha = gh_params[2], 
                                  mu = gh_params[3], 
                                  sigma = gh_params[4], 
                                  gamma = gh_params[5])
                    )
# Computes p-value of two-sided KS test: log-returns vs fitted GH samples
KS2sample(log_returns, gh_samples, conservative = F)

# Two-sample Kuiper test
kuiper_gh = kuiper_test(log_returns, meixner_samples)
summary(kuiper_gh)
