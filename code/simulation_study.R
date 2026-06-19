# =============================================================================
# Simulation study: AIPW vs TMLE vs Residual-on-Residual regression
# under a partially linear data-generating model.
#
# =============================================================================

library(tmle)
library(lmtest)
library(sandwich)
library(here)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(future.apply)
library(progressr)

# ── Simulation parameters ───────────────────────────────────────────────────
n_sim     <- 10000
n_obs     <- 1000

# Scenario grid: sweep the true treatment effect and the outcome-error SD.
# Same seeds are reused across scenarios (common random numbers) so scenario
# contrasts share Monte Carlo noise.
psi_grid     <- c(-1, -0.5, 0, 0.5, 1)
sigma_Y_grid <- c(0.5, 1, 2)
scenarios    <- expand.grid(psi     = psi_grid,
                            sigma_Y = sigma_Y_grid,
                            KEEP.OUT.ATTRS = FALSE)

# Monte Carlo seeds: take the first `n_sim` values from the project's seed
# bank so runs are reproducible and seed assignment is stable across machines.
seed_bank <- read.csv(here("data", "random_seed_values.csv"))$random_seeds
if (length(seed_bank) < n_sim) {
  stop(sprintf("random_seed_values.csv has %d seeds; need %d.",
               length(seed_bank), n_sim))
}
sim_seeds <- as.integer(seed_bank[seq_len(n_sim)])

# ── Data-generating process ─────────────────────────────────────────────────
expit <- function(x) 1 / (1 + exp(-x))

generate_data <- function(n, psi, sigma_Y) {
  C1  <- runif(n, -1, 1)
  C2  <- runif(n, -1, 1)
  C3  <- runif(n, -1, 1)
  C4  <- runif(n, -1, 1)
  C5  <- rnorm(n)
  C6  <- rnorm(n)
  C7  <- rnorm(n)
  C8  <- rbinom(n, 1, 0.4)
  C9  <- rbinom(n, 1, 0.5)
  C10 <- rbinom(n, 1, 0.3)
  
  C <- cbind(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10)

  # propensity score: sinusoidal + curvilinear + interactions
  eta_A <- -0.4 +
    0.7 * sin(pi * C1) +
    0.5 * (C2^2) -
    0.6 * C3 * C4 +
    0.4 * C5 -
    0.5 * cos(C6) +
    0.3 * C7 * C5 +
    0.6 * C8 -
    0.4 * C9 +
    0.5 * C10 +
    0.3 * C1 * C8
  pi_C <- expit(eta_A)
  A    <- rbinom(n, 1, pi_C)

  # g(C): curvilinear + interactions
  gC <-  1.5 * sin(2 * C1) +
         1.0 * (C2^2) -
         0.8 * C3 * C4 +
         1.2 * C5 -
         0.9 * cos(C6) +
         0.6 * C7 * C5 +
         0.7 * C8 -
         0.5 * C9 +
         0.5 * C10 +
         0.4 * C1 * C8

  Y <- 125 + psi * A + gC + rnorm(n, 0, sigma_Y)

  data.frame(Y = Y, A = A,
             C1 = C1, C2 = C2, C3 = C3, C4 = C4,
             C5 = C5, C6 = C6, C7 = C7,
             C8 = C8, C9 = C9, C10 = C10)
}

# ── AIPW score ──────────────────────────────────────────────────────────────
aipw_func <- function(exposure, outcome, pscore, mu_hat, mu_hat0, mu_hat1) {
  ((2 * exposure - 1) * (outcome - mu_hat)) /
    ((2 * exposure - 1) * pscore + (1 - exposure)) +
    (mu_hat1 - mu_hat0)
}

# ── One Monte Carlo Replication ─────────────────────────────────────────────
# Nuisance models are fit with the *correct* (oracle) functional forms taken
# directly from generate_data(). Because the parametric models are correctly
# specified there's no overfitting risk, so neither SuperLearner nor cross-
# fitting is needed -- single full-data fits give unbiased nuisances. This is
# the oracle benchmark: if RonR/AIPW/TMLE don't behave well here, the issue is
# with the estimator, not the nuisance ML.
run_one_sim <- function(seed, psi, sigma_Y) {

  set.seed(seed)

  dat  <- generate_data(n_obs, psi, sigma_Y)
  covs <- dat[, paste0("C", 1:10)]

  # Propensity pi(C): matches eta_A in generate_data (note sin(pi*C1)).
  fit_pi <- glm(
    A ~ sin(pi * C1) + I(C2^2) + I(C3 * C4) + C5 + cos(C6) + I(C7 * C5)
        + C8 + C9 + C10 + I(C1 * C8),
    family = binomial(),
    data   = dat
  )
  pi_hat <- as.numeric(predict(fit_pi, type = "response"))
  pi_hat <- pmin(pmax(pi_hat, 0.01), 0.99)

  # Outcome mu(A, C): matches gC in generate_data (sin(2*C1)) plus linear A.
  fit_mu <- lm(
    Y ~ A + sin(2 * C1) + I(C2^2) + I(C3 * C4) + C5 + cos(C6) + I(C7 * C5)
        + C8 + C9 + C10 + I(C1 * C8),
    data = dat
  )
  mu_hat  <- as.numeric(predict(fit_mu))
  mu_hat0 <- as.numeric(predict(fit_mu, newdata = base::transform(dat, A = 0)))
  mu_hat1 <- as.numeric(predict(fit_mu, newdata = base::transform(dat, A = 1)))

  # Outcome nu(C) = E[Y | C]: uses the gC basis (no A). Strictly,
  # E[Y|C] = 125 + psi*pi(C) + g(C), so this absorbs psi*pi(C) into the
  # intercept + basis. Add sin(pi*C1) to the formula if you want it exact.
  fit_nu <- lm(
    Y ~ sin(2 * C1) + I(C2^2) + I(C3 * C4) + C5 + cos(C6) + I(C7 * C5)
        + C8 + C9 + C10 + I(C1 * C8),
    data = dat
  )
  nu_hat <- as.numeric(predict(fit_nu))

  uY <- dat$Y - nu_hat
  uA <- dat$A - pi_hat
  
  # ── Oracle OLS referent ──
  # fit_mu is the correctly-specified Y ~ A + gC-basis regression, so the
  # coefficient on A is the OLS estimate of psi and its model-based
  # (homoskedastic) SE is the textbook standard error.
  ols_coef <- summary(fit_mu)$coefficients
  ols_psi  <- unname(ols_coef["A", "Estimate"])
  ols_se   <- unname(ols_coef["A", "Std. Error"])

  # ── Residual-on-residual (no-intercept OLS, HC3 SE) ──
  fit_ror <- lm(uY ~ 0 + uA)
  ct      <- coeftest(fit_ror, vcov = vcovHC(fit_ror, type = "HC3"))
  ronr_psi <- unname(ct["uA", "Estimate"])
  ronr_se  <- unname(ct["uA", "Std. Error"])

  # ── AIPW ──
  score    <- aipw_func(dat$A, dat$Y, pi_hat, mu_hat, mu_hat0, mu_hat1)
  aipw_psi <- mean(score)
  aipw_se  <- sd(score) / sqrt(n_obs)

  # ── TMLE (continuous Y -> family = "gaussian") ──
  fit_tmle <- tryCatch(
    tmle(Y      = dat$Y,
         A      = dat$A,
         W      = covs,
         Q      = cbind(mu_hat0, mu_hat1),
         g1W    = pi_hat,
         Q.SL.library = NULL, g.SL.library = NULL,
         family = "gaussian"),
    error = function(e) NULL
  )
  if (!is.null(fit_tmle)) {
    tmle_psi <- fit_tmle$estimates$ATE$psi
    tmle_se  <- sqrt(fit_tmle$estimates$ATE$var.psi)
  } else {
    tmle_psi <- NA_real_
    tmle_se  <- NA_real_
  }

  result <- tibble(
    seed     = seed,
    psi_true = psi,
    sigma_Y  = sigma_Y,
    method   = c("OLS", "RonR", "AIPW", "TMLE"),
    psi      = c(ols_psi, ronr_psi, aipw_psi, tmle_psi),
    se       = c(ols_se,  ronr_se,  aipw_se,  tmle_se)
  )

  result
}

# ── Parallel driver ─────────────────────────────────────────────────────────

plan(multisession, workers = max(1, parallel::detectCores() - 1))

handlers(global = TRUE)
handlers("cli")

# Output path is set up before any sims run so a bad path fails fast, and so
# the periodic flushes below have somewhere to write.
out_dir  <- here("output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
csv_path <- file.path(out_dir, "sim_results.csv")

# Flush the accumulated results to CSV every `flush_every` replicates so a
# crash loses at most one chunk's worth of work. The CSV is overwritten in
# full each flush, which is cheap at this row count.
flush_every <- 500
chunks      <- split(sim_seeds, ceiling(seq_along(sim_seeds) / flush_every))

done <- tibble()
with_progress({
  p <- progressor(steps = nrow(scenarios) * length(sim_seeds))
  for (s in seq_len(nrow(scenarios))) {
    psi_s   <- scenarios$psi[s]
    sigma_s <- scenarios$sigma_Y[s]
    message(sprintf("[%s] scenario %d / %d: psi = %g, sigma_Y = %g",
                    Sys.time(), s, nrow(scenarios), psi_s, sigma_s))
    for (k in seq_along(chunks)) {
      t0 <- Sys.time()
      chunk_results <- future_lapply(chunks[[k]], function(seed) {
        out <- run_one_sim(seed, psi_s, sigma_s)
        p(sprintf("psi=%g sigma=%g seed %d", psi_s, sigma_s, seed))
        out
      }, future.seed = TRUE)
      done <- bind_rows(done, bind_rows(chunk_results))
      write.csv(done, csv_path, row.names = FALSE)
      message(sprintf("[%s]   chunk %d / %d flushed (%d rows total) in %.1f min",
                      Sys.time(), k, length(chunks),
                      nrow(done),
                      as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }
})

plan(sequential)

message(sprintf("Done. %s contains %d rows.", csv_path, nrow(done)))

# ── Performance summary ─────────────────────────────────────────────────────

z975 <- qnorm(0.975)

summary_tbl <- done %>%
  group_by(psi_true, sigma_Y, method) %>%
  summarise(
    n_reps   = n(),
    bias     = mean(psi - psi_true, na.rm = TRUE),
    emp_se   = sd(psi, na.rm = TRUE),
    avg_se   = mean(se, na.rm = TRUE),
    se_ratio = avg_se / emp_se,
    rmse     = sqrt(mean((psi - psi_true)^2, na.rm = TRUE)),
    cover_95 = mean(
      (psi - z975 * se) <= psi_true &
      (psi + z975 * se) >= psi_true,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(summary_tbl)

saveRDS(summary_tbl, file.path(out_dir, "sim_summary.rds"))

# ── Performance plots ───────────────────────────────────────────────────────
# One faceted figure: rows = metric (Bias, RMSE, 95% coverage), columns =
# sigma_Y, x-axis = psi, colour/shape = method. Dashed reference lines at
# bias = 0 and coverage = 0.95 (RMSE has no nominal target).

fig_dir <- file.path(out_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

plot_df <- summary_tbl %>%
  select(psi_true, sigma_Y, method, bias, rmse, cover_95) %>%
  mutate(cover_95 = 100 * cover_95) %>%
  pivot_longer(cols      = c(bias, rmse, cover_95),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(
    metric = factor(metric,
                    levels = c("bias", "rmse", "cover_95"),
                    labels = c("Bias", "RMSE", "95% coverage")),
    method = factor(method, levels = c("OLS", "RonR", "AIPW", "TMLE"))
  )

metric_levels <- c("Bias", "RMSE", "95% coverage")

hline_df <- tibble(
  metric = factor(c("Bias", "95% coverage"), levels = metric_levels),
  yint   = c(0, 95)
)

# Pin per-row y-limits via blank points keyed by `metric`. With
# scales = "free_y" this fixes the range exactly when the data fits inside;
# if any value strays outside, the range expands rather than silently
# clipping (so an out-of-range estimate is visible).
limits_df <- tibble(
  metric = factor(rep(metric_levels, each = 2), levels = metric_levels),
  value  = c(-0.05, 0.05,   0, 0.15,   90, 100)
)

p <- ggplot(plot_df, aes(x = psi_true, y = value,
                         colour = method, shape = method)) +
  geom_blank(data = limits_df, aes(y = value), inherit.aes = FALSE) +
  geom_hline(data = hline_df, aes(yintercept = yint),
             linetype = "dashed", colour = "grey50",
             inherit.aes = FALSE) +
  geom_line() +
  geom_point(size = 1.8) +
  facet_grid(metric ~ sigma_Y,
             scales   = "free_y",
             labeller = label_bquote(cols = sigma[Y] == .(sigma_Y))) +
  scale_x_continuous(breaks = psi_grid) +
  scale_colour_grey(start = 0.05, end = 0.7) +
  scale_shape_manual(values = c(OLS = 16, RonR = 17, AIPW = 15, TMLE = 18)) +
  labs(x = expression(psi),
       y = NULL,
       colour = "Method",
       shape  = "Method") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "sim_performance.pdf"),
       plot = p, width = 9, height = 7)
ggsave(file.path(fig_dir, "sim_performance.png"),
       plot = p, width = 9, height = 7, dpi = 200)

message(sprintf("Wrote %s and %s",
                file.path(fig_dir, "sim_performance.pdf"),
                file.path(fig_dir, "sim_performance.png")))
