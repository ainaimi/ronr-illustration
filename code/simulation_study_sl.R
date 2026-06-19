# =============================================================================
# Simulation study (SuperLearner GLM-library variant): AIPW vs TMLE vs RonR
# under a smooth-nonlinear DGP, with nuisances fit by a 4-member SuperLearner
# library of GLM variants.
#
# Library:
#   SL.glm              -- main effects only (matches the naive OLS comparator)
#   SL.glm.interaction  -- main effects + all two-way interactions
#   SL.glm.quad         -- main effects + quadratics for continuous columns
#   SL.glm.full         -- main effects + quadratics + two-way interactions
#                          (correctly specified for this DGP)
#
# All learners are GLM-based, so nuisance estimation is orders of magnitude
# faster than a ranger/xgboost-based library while still giving SL access to a
# correctly-specified working model via SL.glm.full.
#
# Contrasts with simulation_study.R (the oracle-parametric variant):
#   - DGP uses smooth polynomial + bilinear-interaction terms (no sin/cos
#     oscillation), so SL.glm.full recovers it via standard basis expansion.
#   - Nuisances pi, mu, nu are each fit by CV.SuperLearner sharing the same
#     outer folds, so cross-fitted predictions plug cleanly into AIPW/TMLE/RonR.
#   - The OLS comparator is the *naive* additive Y ~ A + C1 + ... + Cp
#     (intentionally misspecified) so it actually shows the parametric bias.
#   - Per-replicate SL ensemble weights are recorded for each nuisance so we
#     can audit which library member SL is leaning on.
# =============================================================================

library(tmle)
library(lmtest)
library(sandwich)
library(SuperLearner)
library(here)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(future.apply)
library(progressr)

# ── Custom SL learners (GLM variants) ───────────────────────────────────────
# Both learners decide column-by-column whether each variable in `X` is
# continuous (>2 unique values and not all in {0,1}) and store the resulting
# `cont_vars` on the fit so prediction-time augmentation is consistent with
# training-time augmentation. When X = cbind(A, covs) is passed (the mu fit),
# the A column is treated as binary and no A^2 or A:C^2 terms are added,
# although SL.glm.full does include linear A:C interactions.

SL.glm.quad <- function(Y, X, newX, family, obsWeights, ...) {
  if (is.matrix(X))    X    <- as.data.frame(X)
  if (is.matrix(newX)) newX <- as.data.frame(newX)

  is_cont <- function(x) {
    is.numeric(x) && length(unique(x)) > 2 && !all(x %in% c(0, 1))
  }
  cont_vars <- colnames(X)[vapply(X, is_cont, logical(1))]

  X_aug    <- X
  newX_aug <- newX
  for (nm in cont_vars) {
    X_aug[[paste0(nm, "_sq")]]    <- X[[nm]]^2
    newX_aug[[paste0(nm, "_sq")]] <- newX[[nm]]^2
  }

  fit.glm <- glm(Y ~ ., data = X_aug, family = family, weights = obsWeights)
  pred    <- predict(fit.glm, newdata = newX_aug, type = "response")
  fit     <- list(object = fit.glm, cont_vars = cont_vars)
  class(fit) <- "SL.glm.quad"
  list(pred = pred, fit = fit)
}

predict.SL.glm.quad <- function(object, newdata, ...) {
  if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
  newdata_aug <- newdata
  for (nm in object$cont_vars) {
    if (nm %in% colnames(newdata)) {
      newdata_aug[[paste0(nm, "_sq")]] <- newdata[[nm]]^2
    }
  }
  predict(object$object, newdata = newdata_aug, type = "response")
}

# SL.glm.full: quadratics for continuous columns + two-way interactions
# between the original input columns (no quadratic-by-anything interactions,
# which would burn d.f. on C^3 and C^2*C' terms the DGP doesn't have).
SL.glm.full <- function(Y, X, newX, family, obsWeights, ...) {
  if (is.matrix(X))    X    <- as.data.frame(X)
  if (is.matrix(newX)) newX <- as.data.frame(newX)

  is_cont <- function(x) {
    is.numeric(x) && length(unique(x)) > 2 && !all(x %in% c(0, 1))
  }
  cont_vars <- colnames(X)[vapply(X, is_cont, logical(1))]

  X_aug    <- X
  newX_aug <- newX
  for (nm in cont_vars) {
    X_aug[[paste0(nm, "_sq")]]    <- X[[nm]]^2
    newX_aug[[paste0(nm, "_sq")]] <- newX[[nm]]^2
  }

  orig        <- colnames(X)
  inter_terms <- combn(orig, 2, FUN = function(p) paste(p, collapse = ":"))
  rhs         <- paste(c(colnames(X_aug), inter_terms), collapse = " + ")
  fmla        <- as.formula(paste("Y ~", rhs))

  fit.glm <- glm(fmla, data = X_aug, family = family, weights = obsWeights)
  pred    <- predict(fit.glm, newdata = newX_aug, type = "response")
  fit     <- list(object = fit.glm, cont_vars = cont_vars)
  class(fit) <- "SL.glm.full"
  list(pred = pred, fit = fit)
}

predict.SL.glm.full <- function(object, newdata, ...) {
  if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
  newdata_aug <- newdata
  for (nm in object$cont_vars) {
    if (nm %in% colnames(newdata)) {
      newdata_aug[[paste0(nm, "_sq")]] <- newdata[[nm]]^2
    }
  }
  predict(object$object, newdata = newdata_aug, type = "response")
}

# Short labels for the four library members, used to name the recorded SL
# weight columns. Order matches the `sl_lib` vector in run_one_sim.
sl_lib_short <- c("SL.glm"             = "glm",
                  "SL.glm.interaction" = "int",
                  "SL.glm.quad"        = "quad",
                  "SL.glm.full"        = "full")

# ── Simulation parameters ───────────────────────────────────────────────────
# Each SL replicate is ~150 GLM fits per nuisance times 3 nuisances; with all
# learners GLM-based, per-replicate cost is dominated by glm() and is small
# relative to the oracle parametric script.
n_sim     <- 10000
n_obs     <- 500
num_folds <- 5     # outer CV folds (shared across the three SL nuisance fits)

# Scenario grid: sweep psi and outcome-error SD. Same seeds across scenarios
# (common random numbers) so scenario contrasts share Monte Carlo noise.
psi_grid     <- c(-1, -0.5, 0, 0.5, 1)
sigma_Y_grid <- c(0.5, 1, 2)
scenarios    <- expand.grid(psi     = psi_grid,
                            sigma_Y = sigma_Y_grid,
                            KEEP.OUT.ATTRS = FALSE)

seed_bank <- read.csv(here("data", "random_seed_values.csv"))$random_seeds
if (length(seed_bank) < n_sim) {
  stop(sprintf("random_seed_values.csv has %d seeds; need %d.",
               length(seed_bank), n_sim))
}
sim_seeds <- as.integer(seed_bank[seq_len(n_sim)])

# ── Data-generating process ─────────────────────────────────────────────────
# Six covariates: C1, C2 uniform; C3, C4 normal; C5, C6 binary. C6 is noise
# (enters neither pi nor mu); the rest are confounders or partial confounders.
# pi(C) and gC are quadratic in C1 plus a single bilinear interaction
# (C2*C3 in both) plus linear main effects. Smooth + low-order interactions:
# tractable for ranger, missed by linear-additive OLS.
expit <- function(x) 1 / (1 + exp(-x))

generate_data <- function(n, psi, sigma_Y) {
  C1 <- runif(n, -1, 1)
  C2 <- runif(n, -1, 1)
  C3 <- rnorm(n)
  C4 <- rnorm(n)
  C5 <- rbinom(n, 1, 0.4)
  C6 <- rbinom(n, 1, 0.5)

  eta_A <- -0.3 +
            0.9 * C1 +
            1.0 * (C1^2) -
            0.7 * C2 +
            0.8 * C2 * C3 -
            0.5 * C3 +
            0.6 * C5 -
            0.4 * C5 * C4
  pi_C <- expit(eta_A)
  A    <- rbinom(n, 1, pi_C)

  gC <- 1.2 * C1 +
        1.5 * (C1^2) -
        1.0 * C2 +
        1.4 * C2 * C3 -
        0.8 * C3 +
        0.7 * C4 +
        1.0 * C5 -
        0.6 * C5 * C4

  Y <- 125 + psi * A + gC + rnorm(n, 0, sigma_Y)

  data.frame(Y = Y, A = A,
             C1 = C1, C2 = C2, C3 = C3, C4 = C4, C5 = C5, C6 = C6)
}

# ── AIPW score ──────────────────────────────────────────────────────────────
aipw_func <- function(exposure, outcome, pscore, mu_hat, mu_hat0, mu_hat1) {
  ((2 * exposure - 1) * (outcome - mu_hat)) /
    ((2 * exposure - 1) * pscore + (1 - exposure)) +
    (mu_hat1 - mu_hat0)
}

# ── One Monte Carlo Replication ─────────────────────────────────────────────
# All three nuisances (pi, mu, nu) fit by CV.SuperLearner on the *same* outer
# folds via shared validRows, so cross-fitted predictions for AIPW/TMLE/RonR
# are consistently ordered. Counterfactual mu_hat0/mu_hat1 are assembled by
# positional assignment (NOT rbind in fold order — that bug silently mixes
# observations and biases AIPW/TMLE; see 2026-06-16 notes).
run_one_sim <- function(seed, psi, sigma_Y) {

  # Anchor the custom SL learners and their S3 predict methods so that
  # future.apply's globalsOf auto-detection picks them up and exports them to
  # parallel workers. Without this, the strings in `sl_lib` are not statically
  # detectable as globals and S3 dispatch on the class "SL.glm.quad" /
  # "SL.glm.full" would fail in the worker process.
  invisible(list(SL.glm.quad, predict.SL.glm.quad,
                 SL.glm.full, predict.SL.glm.full,
                 sl_lib_short))

  set.seed(seed)

  dat  <- generate_data(n_obs, psi, sigma_Y)
  covs <- dat[, paste0("C", 1:6)]

  # Shared outer folds for all three nuisance fits. `fold_index` is a list of
  # length `num_folds`; element k holds the validation-row indices for fold k.
  # Passed to each CV.SuperLearner as `cvControl$validRows`.
  fold_ids    <- sample(rep(seq_len(num_folds), length.out = n_obs))
  fold_index  <- split(seq_len(n_obs), fold_ids)

  sl_lib <- names(sl_lib_short)

  # Helper: collapse per-fold SL weights into one named vector in canonical
  # library order. CV.SuperLearner's `coef` is a V-by-library matrix whose
  # column names look like "SL.glm_All" (the "_All" suffix is the default
  # screen name); we strip the suffix before matching to sl_lib_short.
  pool_weights <- function(fit) {
    w <- colMeans(fit$coef)
    names(w) <- sub("_All$", "", names(w))
    w <- w[names(sl_lib_short)]
    setNames(unname(w), unname(sl_lib_short))
  }

  # Propensity pi(C) = P(A=1 | C).
  fit_pi <- CV.SuperLearner(
    Y              = dat$A,
    X              = covs,
    family         = binomial(),
    SL.library     = sl_lib,
    cvControl      = list(V = 5, validRows = fold_index),
    innerCvControl = list(list(V = 5)),
    control        = list(saveCVFitLibrary = FALSE),
    parallel       = "seq",
    verbose        = FALSE
  )
  pi_hat <- as.numeric(fit_pi$SL.predict)
  pi_hat <- pmin(pmax(pi_hat, 0.01), 0.99)
  pi_w   <- pool_weights(fit_pi)

  # Outcome mu(A, C) = E[Y | A, C]. A is a column of X.
  fit_mu <- CV.SuperLearner(
    Y              = dat$Y,
    X              = cbind(A = dat$A, covs),
    family         = gaussian(),
    SL.library     = sl_lib,
    cvControl      = list(V = 5, validRows = fold_index),
    innerCvControl = list(list(V = 5)),
    control        = list(saveCVFitLibrary = TRUE),
    parallel       = "seq",
    verbose        = FALSE
  )
  mu_hat <- as.numeric(fit_mu$SL.predict)
  mu_w   <- pool_weights(fit_mu)

  # Counterfactual predictions: for each outer fold k, predict on its
  # validation rows with A forced to 0 and 1, using the SL ensemble trained
  # on fold k's training data. Positional assignment keeps things in obs order.
  mu_hat0 <- numeric(n_obs)
  mu_hat1 <- numeric(n_obs)
  for (k in seq_len(num_folds)) {
    idx <- fold_index[[k]]
    nd0 <- cbind(A = 0L, covs[idx, , drop = FALSE])
    nd1 <- cbind(A = 1L, covs[idx, , drop = FALSE])
    mu_hat0[idx] <- predict(fit_mu$AllSL[[k]], newdata = nd0)$pred
    mu_hat1[idx] <- predict(fit_mu$AllSL[[k]], newdata = nd1)$pred
  }

  # Marginal outcome nu(C) = E[Y | C].
  fit_nu <- CV.SuperLearner(
    Y              = dat$Y,
    X              = covs,
    family         = gaussian(),
    SL.library     = sl_lib,
    cvControl      = list(V = 5, validRows = fold_index),
    innerCvControl = list(list(V = 5)),
    control        = list(saveCVFitLibrary = FALSE),
    parallel       = "seq",
    verbose        = FALSE
  )
  nu_hat <- as.numeric(fit_nu$SL.predict)
  nu_w   <- pool_weights(fit_nu)

  uY <- dat$Y - nu_hat
  uA <- dat$A - pi_hat

  # ── Naive linear OLS (misspecified comparator) ──
  # Linear-additive in C; misses the C1^2 and C2*C3 structure of the DGP, so
  # under nonlinear confounding the A coefficient absorbs bias. This is the
  # parametric strawman the SL-based estimators are meant to beat.
  fit_ols  <- lm(Y ~ A + C1 + C2 + C3 + C4 + C5 + C6, data = dat)
  ols_coef <- summary(fit_ols)$coefficients
  ols_psi  <- unname(ols_coef["A", "Estimate"])
  ols_se   <- unname(ols_coef["A", "Std. Error"])

  # ── Residual-on-residual (no-intercept OLS, HC3 SE) ──
  fit_ror  <- lm(uY ~ 0 + uA)
  ct       <- coeftest(fit_ror, vcov = vcovHC(fit_ror, type = "HC3"))
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

  out <- tibble(
    seed     = seed,
    psi_true = psi,
    sigma_Y  = sigma_Y,
    method   = c("OLS", "RonR", "AIPW", "TMLE"),
    psi      = c(ols_psi, ronr_psi, aipw_psi, tmle_psi),
    se       = c(ols_se,  ronr_se,  aipw_se,  tmle_se)
  )

  # Attach SL ensemble weights for each nuisance fit. Weights are per-replicate,
  # so the same value is repeated across the 4 method rows; downstream summaries
  # de-duplicate by filtering to one method.
  for (nm in names(pi_w)) out[[paste0("w_pi_", nm)]] <- pi_w[[nm]]
  for (nm in names(mu_w)) out[[paste0("w_mu_", nm)]] <- mu_w[[nm]]
  for (nm in names(nu_w)) out[[paste0("w_nu_", nm)]] <- nu_w[[nm]]

  out
}

# ── Parallel driver ─────────────────────────────────────────────────────────

plan(multisession, workers = max(1, parallel::detectCores() - 1))

handlers(global = TRUE)
handlers("cli")

out_dir  <- here("output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
csv_path <- file.path(out_dir, "sim_results_sl.csv")

flush_every <- 200
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
      chunk_results <- future_lapply(
        chunks[[k]],
        function(seed) {
          out <- run_one_sim(seed, psi_s, sigma_s)
          p(sprintf("psi=%g sigma=%g seed %d", psi_s, sigma_s, seed))
          out
        },
        future.seed     = TRUE,
        future.packages = c("SuperLearner", "tmle")
      )
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

saveRDS(summary_tbl, file.path(out_dir, "sim_summary_sl.rds"))
write.csv(summary_tbl, file.path(out_dir, "sim_summary_sl.csv"), row.names = FALSE)

# ── SL ensemble weights summary ─────────────────────────────────────────────
# Mean per-scenario SL weights for each of pi, mu, nu. Weights are recorded
# per-replicate and repeated across the 4 method rows, so we filter to one
# method before aggregating to avoid double counting.

weights_tbl <- done %>%
  filter(method == "AIPW") %>%
  group_by(psi_true, sigma_Y) %>%
  summarise(n_reps = n(),
            across(starts_with("w_"), mean),
            .groups = "drop")

print(weights_tbl)

saveRDS(weights_tbl, file.path(out_dir, "sim_summary_weights_sl.rds"))
write.csv(weights_tbl, file.path(out_dir, "sim_summary_weights_sl.csv"),
          row.names = FALSE)

# Quick diagnostic: mean weight on each library member, marginalised over all
# scenarios. If SL.glm.full is doing the heavy lifting, expect w_*_full close
# to 1 here. If the ensemble is splitting, SL is finding multiple usable fits.
overall_weights <- weights_tbl %>%
  summarise(across(starts_with("w_"), mean))
message("Overall mean SL weights (marginalised over scenarios):")
print(overall_weights)

# ── Performance plots ───────────────────────────────────────────────────────
# Rows = metric (Bias, RMSE, 95% coverage), columns = sigma_Y, x-axis = psi,
# colour/shape = method. Dashed reference lines at bias = 0 and coverage = 95.
# Bias y-range is widened relative to the oracle script to accommodate the
# misspecified-OLS bias.

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

limits_df <- tibble(
  metric = factor(rep(metric_levels, each = 2), levels = metric_levels),
  value  = c(-0.5, 0.5,   0, 0.6,   80, 100)
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

ggsave(file.path(fig_dir, "sim_performance_sl.pdf"),
       plot = p, width = 9, height = 7)
ggsave(file.path(fig_dir, "sim_performance_sl.png"),
       plot = p, width = 9, height = 7, dpi = 200)

message(sprintf("Wrote %s and %s",
                file.path(fig_dir, "sim_performance_sl.pdf"),
                file.path(fig_dir, "sim_performance_sl.png")))
