# =============================================================================
# POSITIVITY SWEEP (Story 1): AIPW vs TMLE vs RonR as overlap degrades.
#
# This is the positivity-stress variant of simulation_study_sl.R. The effect is
# held HOMOGENEOUS (a constant psi for every unit, no A-by-C interaction), so
# the population ATE is exactly psi regardless of overlap. We then sweep a
# single overlap-severity dial `zeta` that steepens the propensity linear
# predictor, pushing pi(C) toward the 0/1 boundaries. Because the effect is
# constant, RonR, AIPW, and TMLE all still TARGET the same psi; any divergence
# in their finite-sample behaviour as zeta grows is therefore instability, not
# a change of estimand. The claim under test (manuscript Discussion): RonR
# residualizes the exposure instead of inverting it, so it should stay stable
# while the inverse-propensity machinery in AIPW/TMLE destabilizes.
#
# What changes relative to simulation_study_sl.R:
#   - Scenario grid sweeps `zeta` (overlap severity) at a single fixed psi and
#     sigma_Y, rather than crossing psi x sigma_Y.
#   - generate_data() takes `zeta` and returns the true pi(C) so we can measure
#     the realized degree of positivity violation per replicate.
#   - Each replicate records overlap diagnostics (true-pi tail mass, estimated-
#     pi extremes, the largest realized IPW weight) and a TMLE-failure flag, so
#     the *magnitude* of the violation can be put on the same x-axis as the
#     instability it induces.
#   - The truncation bound applied to pi_hat is an explicit parameter (`pi_trunc`)
#     so its stabilizing role is visible and tunable rather than hard-wired.
#
# Library (unchanged): SL.glm, SL.glm.interaction, SL.glm.quad, SL.glm.full;
# only SL.glm.full is correctly specified for this DGP.
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

# Fixed effect and noise scale. Story 1 isolates positivity, so we hold these
# constant and sweep overlap only. psi is the true (homogeneous) ATE.
psi_fixed     <- 0.5
sigma_Y_fixed <- 1

# Truncation bound applied to estimated pi_hat before it enters AIPW/TMLE.
# pi_hat is clipped to [pi_trunc, 1 - pi_trunc]. This is the band-aid that
# stabilizes inverse-propensity weighting; we loosen it relative to the main
# simulation (which used 0.01) so the AIPW/TMLE instability the sweep is meant
# to expose is not masked by aggressive clipping. The largest realized IPW
# weight is therefore bounded by 1/pi_trunc.
pi_trunc <- 0.001

# Scenario grid: sweep the overlap-severity dial `zeta`. zeta = 1 reproduces
# the healthy-overlap propensity of the main simulation; larger zeta steepens
# the linear predictor and drives pi(C) toward 0/1. Same seeds across scenarios
# (common random numbers) so scenario contrasts share Monte Carlo noise.
zeta_grid    <- c(1, 2, 3, 4, 6)
scenarios    <- expand.grid(zeta = zeta_grid,
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
# (C2*C3 in both) plus linear main effects.
#
# Overlap dial: the propensity linear predictor is split into a fixed intercept
# (-0.3, holding marginal prevalence roughly centred) and a confounding core
# scaled by `zeta`. zeta = 1 is the healthy-overlap baseline; increasing zeta
# steepens the core, spreading pi(C) toward 0 and 1 and creating progressively
# severe practical positivity violations. The outcome model is UNCHANGED and
# the effect is constant (psi * A, no A-by-C term), so the ATE stays exactly
# psi at every zeta -- only overlap degrades.
expit <- function(x) 1 / (1 + exp(-x))

generate_data <- function(n, psi, sigma_Y, zeta = 1) {
  C1 <- runif(n, -1, 1)
  C2 <- runif(n, -1, 1)
  C3 <- rnorm(n)
  C4 <- rnorm(n)
  C5 <- rbinom(n, 1, 0.4)
  C6 <- rbinom(n, 1, 0.5)

  eta_core <-  0.9 * C1 +
               1.0 * (C1^2) -
        zeta * 0.7 * C2 +
               0.8 * C2 * C3 -
               0.5 * C3 +
               0.6 * C5 -
               0.4 * C5 * C4
  eta_A <- -0.3 + eta_core
  pi_C  <- expit(eta_A)
  A     <- rbinom(n, 1, pi_C)

  gC <- 1.2 * C1 +
        1.5 * (C1^2) -
        1.0 * C2 +
        1.4 * C2 * C3 -
        0.8 * C3 +
        0.7 * C4 +
        1.0 * C5 -
        0.6 * C5 * C4

  Y <- 125 + psi * A + gC + rnorm(n, 0, sigma_Y)

  # pi_true is carried out (not used by any estimator) so each replicate can
  # report the realized degree of positivity violation it actually drew.
  data.frame(Y = Y, A = A,
             C1 = C1, C2 = C2, C3 = C3, C4 = C4, C5 = C5, C6 = C6,
             pi_true = pi_C)
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
run_one_sim <- function(seed, zeta, psi, sigma_Y, pi_trunc) {

  # Anchor the custom SL learners and their S3 predict methods so that
  # future.apply's globalsOf auto-detection picks them up and exports them to
  # parallel workers. Without this, the strings in `sl_lib` are not statically
  # detectable as globals and S3 dispatch on the class "SL.glm.quad" /
  # "SL.glm.full" would fail in the worker process.
  invisible(list(SL.glm.quad, predict.SL.glm.quad,
                 SL.glm.full, predict.SL.glm.full,
                 sl_lib_short))

  set.seed(seed)

  dat  <- generate_data(n_obs, psi, sigma_Y, zeta)
  covs <- dat[, paste0("C", 1:6)]

  # Realized positivity violation for this replicate, measured on the TRUE
  # propensity pi_true (the ground-truth overlap we drew, independent of any
  # estimation error). These anchor the x-axis: they convert the abstract dial
  # `zeta` into the actual tail mass / boundary proximity it produced.
  pi_true            <- dat$pi_true
  true_pi_min        <- min(pi_true)
  true_pi_max        <- max(pi_true)
  prop_true_lt01     <- mean(pi_true < 0.01  | pi_true > 0.99)
  prop_true_out0595  <- mean(pi_true < 0.05  | pi_true > 0.95)

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
  pi_hat_raw <- as.numeric(fit_pi$SL.predict)
  # Estimated-pi extremes BEFORE truncation: how far the fitted propensity ran
  # toward the boundary (estimation instability, distinct from the true pi).
  pihat_raw_min <- min(pi_hat_raw)
  pihat_raw_max <- max(pi_hat_raw)
  pi_hat <- pmin(pmax(pi_hat_raw, pi_trunc), 1 - pi_trunc)
  n_trunc <- sum(pi_hat_raw < pi_trunc | pi_hat_raw > 1 - pi_trunc)
  # Largest inverse-probability weight actually used (treated weight 1/pi_hat,
  # control weight 1/(1 - pi_hat)). This is the quantity that destabilizes
  # AIPW/TMLE; bounded above by 1/pi_trunc by construction.
  max_ipw <- max(ifelse(dat$A == 1, 1 / pi_hat, 1 / (1 - pi_hat)))
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
    tmle_psi  <- fit_tmle$estimates$ATE$psi
    tmle_se   <- sqrt(fit_tmle$estimates$ATE$var.psi)
    tmle_fail <- FALSE
  } else {
    tmle_psi  <- NA_real_
    tmle_se   <- NA_real_
    tmle_fail <- TRUE   # tracked: na.rm in the summary would silently hide this
  }

  out <- tibble(
    seed     = seed,
    zeta     = zeta,
    psi_true = psi,
    sigma_Y  = sigma_Y,
    method   = c("OLS", "RonR", "AIPW", "TMLE"),
    psi      = c(ols_psi, ronr_psi, aipw_psi, tmle_psi),
    se       = c(ols_se,  ronr_se,  aipw_se,  tmle_se)
  )

  # Per-replicate positivity diagnostics. Constant across the 4 method rows
  # (downstream summaries de-duplicate by filtering to one method, as with the
  # SL weights). tmle_fail is method-specific but harmless to carry on all rows.
  out$true_pi_min       <- true_pi_min
  out$true_pi_max       <- true_pi_max
  out$prop_true_lt01    <- prop_true_lt01
  out$prop_true_out0595 <- prop_true_out0595
  out$pihat_raw_min     <- pihat_raw_min
  out$pihat_raw_max     <- pihat_raw_max
  out$n_trunc           <- n_trunc
  out$max_ipw           <- max_ipw
  out$tmle_fail         <- tmle_fail

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
csv_path <- file.path(out_dir, "sim_results_sl_positivity.csv")

flush_every <- 200
chunks      <- split(sim_seeds, ceiling(seq_along(sim_seeds) / flush_every))

done <- tibble()
with_progress({
  p <- progressor(steps = nrow(scenarios) * length(sim_seeds))
  for (s in seq_len(nrow(scenarios))) {
    zeta_s <- scenarios$zeta[s]
    message(sprintf("[%s] scenario %d / %d: zeta = %g (psi = %g, sigma_Y = %g)",
                    Sys.time(), s, nrow(scenarios), zeta_s,
                    psi_fixed, sigma_Y_fixed))
    for (k in seq_along(chunks)) {
      t0 <- Sys.time()
      chunk_results <- future_lapply(
        chunks[[k]],
        function(seed) {
          out <- run_one_sim(seed, zeta_s, psi_fixed, sigma_Y_fixed, pi_trunc)
          p(sprintf("zeta=%g seed %d", zeta_s, seed))
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
  group_by(zeta, method) %>%
  summarise(
    n_reps   = n(),
    n_valid  = sum(!is.na(psi)),       # reps that produced an estimate
    n_fail   = sum(is.na(psi)),        # estimator failures (TMLE under severe zeta)
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

print(summary_tbl, n = Inf)

saveRDS(summary_tbl, file.path(out_dir, "sim_summary_sl_positivity.rds"))
write.csv(summary_tbl, file.path(out_dir, "sim_summary_sl_positivity.csv"),
          row.names = FALSE)

# ── SL ensemble weights summary ─────────────────────────────────────────────
# Mean per-scenario SL weights for each of pi, mu, nu. Weights are recorded
# per-replicate and repeated across the 4 method rows, so we filter to one
# method before aggregating to avoid double counting.

weights_tbl <- done %>%
  filter(method == "AIPW") %>%
  group_by(zeta) %>%
  summarise(n_reps = n(),
            across(starts_with("w_"), mean),
            .groups = "drop")

print(weights_tbl)

saveRDS(weights_tbl, file.path(out_dir, "sim_summary_weights_sl_positivity.rds"))
write.csv(weights_tbl, file.path(out_dir, "sim_summary_weights_sl_positivity.csv"),
          row.names = FALSE)

# ── Positivity diagnostics summary ───────────────────────────────────────────
# Translate the abstract dial `zeta` into the realized magnitude of positivity
# violation it produced, averaged over replicates. This is the table that lets
# the instability in summary_tbl be read against an interpretable overlap
# metric (mean tail mass, mean worst-case true/estimated pi, mean largest IPW
# weight, mean number of truncated observations, and the TMLE failure rate).
# Diagnostics are per-replicate and repeated across the 4 method rows, so we
# filter to one method before aggregating.
posviol_tbl <- done %>%
  filter(method == "AIPW") %>%
  group_by(zeta) %>%
  summarise(
    n_reps             = n(),
    mean_true_pi_min   = mean(true_pi_min),
    mean_true_pi_max   = mean(true_pi_max),
    mean_prop_lt01     = mean(prop_true_lt01),     # P(pi<0.01 or >0.99)
    mean_prop_out0595  = mean(prop_true_out0595),  # P(pi<0.05 or >0.95)
    mean_pihat_raw_min = mean(pihat_raw_min),
    mean_pihat_raw_max = mean(pihat_raw_max),
    mean_n_trunc       = mean(n_trunc),
    mean_max_ipw       = mean(max_ipw),
    .groups = "drop"
  )

# TMLE failure rate by zeta (counted across all reps, not just AIPW rows).
tmle_fail_tbl <- done %>%
  filter(method == "TMLE") %>%
  group_by(zeta) %>%
  summarise(tmle_fail_rate = mean(tmle_fail), .groups = "drop")

posviol_tbl <- left_join(posviol_tbl, tmle_fail_tbl, by = "zeta")

message("Realized positivity violation by zeta:")
print(posviol_tbl, n = Inf)

saveRDS(posviol_tbl, file.path(out_dir, "sim_summary_posviol_sl_positivity.rds"))
write.csv(posviol_tbl,
          file.path(out_dir, "sim_summary_posviol_sl_positivity.csv"),
          row.names = FALSE)

# ── Performance plots ───────────────────────────────────────────────────────
# Story 1 figure: overlap severity on the x-axis, one panel per instability
# metric, colour/shape = method. Because the effect is homogeneous the ATE is
# psi at every zeta, so rising Emp.SE / RMSE and falling coverage for AIPW and
# TMLE -- against a flat RonR -- is the instability the sweep is built to show.

fig_dir <- file.path(out_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

metric_levels <- c("Bias", "RMSE", "SE ratio", "95% coverage")

plot_df <- summary_tbl %>%
  select(zeta, method, bias, rmse, se_ratio, cover_95) %>%
  mutate(cover_95 = 100 * cover_95) %>%
  pivot_longer(cols      = c(bias, rmse, se_ratio, cover_95),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(
    metric = factor(metric,
                    levels = c("bias", "rmse", "se_ratio", "cover_95"),
                    labels = metric_levels),
    method = factor(method, levels = c("OLS", "RonR", "AIPW", "TMLE"))
  )

# Reference lines: unbiasedness, SE ratio = 1, nominal coverage = 95.
hline_df <- tibble(
  metric = factor(c("Bias", "SE ratio", "95% coverage"), levels = metric_levels),
  yint   = c(0, 1, 95)
)

p <- ggplot(plot_df, aes(x = zeta, y = value,
                         colour = method, shape = method)) +
  geom_hline(data = hline_df, aes(yintercept = yint),
             linetype = "dashed", colour = "grey50",
             inherit.aes = FALSE) +
  geom_line() +
  geom_point(size = 1.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = zeta_grid) +
  scale_colour_grey(start = 0.05, end = 0.7) +
  scale_shape_manual(values = c(OLS = 16, RonR = 17, AIPW = 15, TMLE = 18)) +
  labs(x = expression("Overlap severity " * zeta * " (larger = worse positivity)"),
       y = NULL,
       colour = "Method",
       shape  = "Method") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "sim_performance_sl_positivity.pdf"),
       plot = p, width = 9, height = 8)
ggsave(file.path(fig_dir, "sim_performance_sl_positivity.png"),
       plot = p, width = 9, height = 8, dpi = 200)

# ── Instability vs. realized violation magnitude ─────────────────────────────
# The headline relationship the user asked for, made explicit: plot RMSE and
# coverage directly against an interpretable, realized overlap metric (the mean
# largest IPW weight, which scales with how close the worst pi got to 0/1)
# rather than against the abstract dial zeta. Joining posviol_tbl onto the
# performance summary gives every (zeta, method) row its realized x-coordinate.
relate_df <- summary_tbl %>%
  select(zeta, method, rmse, cover_95) %>%
  mutate(cover_95 = 100 * cover_95) %>%
  left_join(select(posviol_tbl, zeta, mean_max_ipw, mean_prop_out0595),
            by = "zeta") %>%
  pivot_longer(cols = c(rmse, cover_95),
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("rmse", "cover_95"),
                    labels = c("RMSE", "95% coverage")),
    method = factor(method, levels = c("OLS", "RonR", "AIPW", "TMLE"))
  )

p2 <- ggplot(relate_df, aes(x = mean_max_ipw, y = value,
                            colour = method, shape = method)) +
  geom_hline(data = tibble(metric = factor("95% coverage",
                                           levels = c("RMSE", "95% coverage")),
                           yint = 95),
             aes(yintercept = yint), linetype = "dashed", colour = "grey50",
             inherit.aes = FALSE) +
  geom_line() +
  geom_point(size = 1.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_x_log10() +
  scale_colour_grey(start = 0.05, end = 0.7) +
  scale_shape_manual(values = c(OLS = 16, RonR = 17, AIPW = 15, TMLE = 18)) +
  labs(x = "Mean largest inverse-probability weight (log scale)",
       y = NULL, colour = "Method", shape = "Method") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"))

ggsave(file.path(fig_dir, "sim_instability_vs_violation_sl_positivity.pdf"),
       plot = p2, width = 7, height = 7)
ggsave(file.path(fig_dir, "sim_instability_vs_violation_sl_positivity.png"),
       plot = p2, width = 7, height = 7, dpi = 200)

message(sprintf("Wrote performance and instability-vs-violation figures to %s",
                fig_dir))
