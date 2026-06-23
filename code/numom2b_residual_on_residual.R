library(haven)
library(dplyr)
library(fastDummies)
library(skimr)
library(tibble)
library(SuperLearner)
library(lmtest)
library(sandwich)
library(ggplot2)
library(here)
# remotes::install_github("tlverse/tmle3")
library(tmle3)
library(tmle)
library(xtable)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

# ── Data preparation ──────────────────────────────────────────────────────────

dat <- haven::read_dta(here("data","numom+hhs_diet_unimputed_v2.dta"))

names(dat)

dim(dat)

with(dat, table(
  exposure_missing = is.na(v_totdens_v1),
  outcome_missing  = is.na(pree_acog)
))

dat <- dat %>%
  filter(!is.na(pree_acog), !is.na(v_totdens_v1))

dat$v_dich <- ifelse(dat$v_totdens_v1 >= 1.25, 1, 0)

exposure   <- "v_dich"
outcomes   <- "pree_acog"
covariates <- c(
  # Sociodemographic
  "momage", "bmiprepreg", "smokerpre", "insurpub", "momeduc4", "momrace4",
  "gravidity", "married", "mombornus", "momenglishspeak", "accult",
  "dadeduc4", "povertypct",
  # Neighborhood / built environment
  "adi_nat_v1", "walk_nat_v1", "povperc_v1", "grocery1k_v1",
  # Medical history (pre-existing)
  "prehtn", "prediab", "apnea_highrisk", "medcat_thyroid",
  # Mental health / stress
  "epds_tot_v1", "stress_tot_v1", "anx_tot",
  # Physical activity / sleep
  "pa_totmetwk_new_v1", "v1_sleep_hrsnight_wkavg", "insomnia_tot",
  # Alcohol
  "alcduring", "v1_alctotwk_3mnthpriorpreg",
  # Site / pregnancy
  "publicsite", "pregplanned", "sleepsat", "puqe_tot",
  # Dietary components
  "hei2010_total_score_v1"
)

a_final <- dat %>%
  dplyr::select(all_of(outcomes), all_of(exposure), all_of(covariates))

factor_names <- c("smokerpre", "insurpub", "momeduc4", "momrace4",
                  "married", "mombornus", "momenglishspeak", "accult", "dadeduc4",
                  "prehtn", "prediab", "apnea_highrisk", "medcat_thyroid", "alcduring")
a_final[, factor_names] <- lapply(a_final[, factor_names], factor)

a_final <- dummy_cols(a_final, ignore_na = TRUE,
                      remove_selected_columns = TRUE, remove_first_dummy = TRUE)
skim(a_final)

# ── Missing data handling ─────────────────────────────────────────────────────

nodes_n <- list(
  W = setdiff(names(a_final), c("pree_acog", "v_dich")),
  A = "v_dich",
  Y = "pree_acog"
)

processed_n   <- process_missing(a_final, nodes_n, complete_nodes = "A")
processed_dat <- processed_n$data
processed_dat$id <- 1:nrow(processed_dat)

# ── Covariate matrix ──────────────────────────────────────────────────────────

covs <- processed_dat %>%
  dplyr::select(-id, -pree_acog, -v_dich)

# ── Cross-validation folds ────────────────────────────────────────────────────

set.seed(123)
n         <- nrow(processed_dat)
num.folds <- 10
folds     <- sample(rep(1:num.folds, length.out = n))
fold_dat   <- tibble(id = seq_len(n), folds)
fold_index <- split(fold_dat$id, fold_dat$folds)

# ── Custom glmnet screen ──────────────────────────────────────────────────────

.SL.require <- function(package,
                        message = paste("loading required package (", package, ") failed", sep = "")) {
  if (!requireNamespace(package, quietly = FALSE)) stop(message, call. = FALSE)
  invisible(TRUE)
}

screen.glmnet1 <- function(Y, X, family, alpha = 1, minscreen = 10,
                           nfolds = 10, nlambda = 100, ...) {
  .SL.require("glmnet")
  if (!is.matrix(X)) X <- model.matrix(~-1 + ., X)
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, lambda = NULL,
                             type.measure = "deviance", nfolds = nfolds,
                             family = family$family, alpha = alpha, nlambda = nlambda)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, increased lambda")
    sumCoef       <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, function(x) sum(x != 0))
    newCut        <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0)
  }
  return(whichVariable)
}

sl.lib <- list(
  c("SL.mean",   "All"),
  c("SL.ranger", "All"),
  c("SL.glmnet", "All"),
  c("SL.glm",    "All"),
  c("SL.earth",  "All")
)

# ── Super Learner: outcome model for R on R ───────────────────────────────────

fit_nu <- CV.SuperLearner(
  Y          = processed_dat$pree_acog,
  X          = covs,
  method     = "method.NNLS",
  family     = binomial(),
  SL.library = sl.lib,
  cvControl  = list(V = 10, validRows = fold_index),
  control    = list(saveCVFitLibrary = FALSE),
  parallel   = "seq",
  verbose    = FALSE
)

nu_hat <- fit_nu$SL.predict
uY     <- processed_dat$pree_acog - nu_hat

# ── Super Learner: outcome model for AIPW & TMLE ──────────────────────────────
covariates_aug <- cbind(covs, v_dich = processed_dat$v_dich)
fit_mu <- CV.SuperLearner(
  Y          = processed_dat$pree_acog,
  X          = covariates_aug,
  method     = "method.NNLS",
  family     = binomial(),
  SL.library = sl.lib,
  cvControl  = list(V = 10, validRows = fold_index),
  control    = list(saveCVFitLibrary = TRUE),
  parallel   = "seq",
  verbose    = FALSE
)

# Natural-A predictions: $SL.predict is cross-fitted and already in original
# observation order — no loop needed.
mu_hat <- as.numeric(fit_mu$SL.predict)

# Counterfactual predictions: pre-allocate and assign POSITIONALLY by fold so
# predictions stay aligned to each individual's row (rbind would reorder to
# fold order and silently mis-pair mu with A/Y/pi downstream).
mu_hat1 <- numeric(n)
mu_hat0 <- numeric(n)
for (i in 1:num.folds) {
  idx <- fold_index[[i]]
  mu_hat1[idx] <- predict(fit_mu$AllSL[[i]],
                          newdata = base::transform(covariates_aug[idx, ], v_dich = 1),
                          onlySL = TRUE)$pred
  mu_hat0[idx] <- predict(fit_mu$AllSL[[i]],
                          newdata = base::transform(covariates_aug[idx, ], v_dich = 0),
                          onlySL = TRUE)$pred
}

# ── Super Learner: exposure model ─────────────────────────────────────────────

fit_pi <- CV.SuperLearner(
  Y          = processed_dat$v_dich,
  X          = covs,
  method     = "method.NNLS",
  family     = binomial(),
  SL.library = sl.lib,
  cvControl  = list(V = 10, validRows = fold_index),
  control    = list(saveCVFitLibrary = FALSE),
  parallel   = "seq",
  verbose    = FALSE
)

summary(fit_pi)
coef(fit_pi)

pi <- fit_pi$SL.predict
uA <- processed_dat$v_dich - pi

summarize(pi)

# ── Propensity score overlap plot ─────────────────────────────────────────────

ps_dat <- data.frame(
  ps      = pi,
  exposed = factor(processed_dat$v_dich, levels = c(0, 1),
                   labels = c("Low veg intake", "High veg intake"))
)

ggplot(ps_dat, aes(x = ps, fill = exposed, color = exposed)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  scale_fill_manual(values  = c("#2166ac", "#d6604d")) +
  scale_color_manual(values = c("#2166ac", "#d6604d")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x    = "Estimated propensity score",
    y    = "Density",
    fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "inside", legend.position.inside = c(0.82, 0.85))

ggsave(here("figures", "ps_overlap.pdf"), width = 6, height = 4)

# ── Standardized mean differences (love plot) ─────────────────────────────────

exposed   <- covs[processed_dat$v_dich == 1, ]
unexposed <- covs[processed_dat$v_dich == 0, ]

smd <- function(x1, x0) {
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
    sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
}

smd_df <- data.frame(
  covariate = names(covs),
  smd       = sapply(names(covs), function(v) smd(exposed[[v]], unexposed[[v]]))
)

ggplot(smd_df, aes(x = smd, y = reorder(covariate, abs(smd)))) +
  geom_vline(xintercept = 0,           color = "gray40", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.1, 0.1), color = "gray60", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 3, color = "#2166ac") +
  labs(x = "Standardized Mean Difference", y = NULL) +
  theme_bw(base_size = 13)

ggsave(here("figures", "smd_love_plot.pdf"), width = 6, height = 5)

# ── Residual-on-residual regression ──────────────────────────────────────────

dat_ror <- data.frame(uY = uY, uA = uA)
fit_ror <- lm(uY ~ uA, data = dat_ror)
summary(fit_ror)

ct <- coeftest(fit_ror, vcov = vcovHC(fit_ror, type = "HC3"))
print(ct)

beta_hat <- ct["uA", "Estimate"]
se_hat   <- ct["uA", "Std. Error"]

c(estimate = beta_hat,
  lower_95  = beta_hat - qnorm(0.975) * se_hat,
  upper_95  = beta_hat + qnorm(0.975) * se_hat)


## AIPW
aipw_func <- function(exposure, outcome,
                      pscore, mu_hat, mu_hat0, mu_hat1) {
  aipw_score <- ((2 * exposure- 1) * (outcome-
                                        mu_hat))/((2 * exposure- 1) * pscore +
                                                    (1- exposure)) + (mu_hat1- mu_hat0)
  return(aipw_score)
}

aipw_score <- aipw_func(exposure = processed_dat$v_dich, 
                        outcome = processed_dat$pree_acog,
                        pscore = pi, 
                        mu_hat = mu_hat, 
                        mu_hat0 = mu_hat0, 
                        mu_hat1 = mu_hat1)

cbind(as.matrix(mu_hat), mu_hat0, mu_hat1)

colnames(aipw_score) <- NULL
## what is aipw_score?
aipw_psi <- mean(aipw_score)
aipw_se <- sd(aipw_score)/sqrt(nrow(processed_dat))
aipw_ate <- c(aipw_psi, aipw_se, aipw_lcl = aipw_psi - 1.96*aipw_se, aipw_Ucl = aipw_psi + 1.96*aipw_se)


fit_tmle <- tmle(Y = processed_dat$pree_acog, A = processed_dat$v_dich,
                 W = covs, Q = cbind(mu_hat0, mu_hat1),
                 g1W = pi, family = "binomial")

# ── Results table ─────────────────────────────────────────────────────────────

tmle_est <- fit_tmle$estimates$ATE$psi
tmle_se  <- sqrt(fit_tmle$estimates$ATE$var.psi)
tmle_ci  <- fit_tmle$estimates$ATE$CI

results <- data.frame(
  Method                          = c("Residual-on-Residual", "AIPW", "TMLE"),
  `Risk Difference (per 100)`     = c(beta_hat, aipw_psi, tmle_est) * 100,
  `Lower 95\\% CL`                = c(beta_hat - qnorm(0.975) * se_hat,
                                      aipw_psi  - qnorm(0.975) * aipw_se,
                                      tmle_ci[1]) * 100,
  `Upper 95\\% CL`                = c(beta_hat + qnorm(0.975) * se_hat,
                                      aipw_psi  + qnorm(0.975) * aipw_se,
                                      tmle_ci[2]) * 100,
  check.names = FALSE
)

results[, -1] <- lapply(results[, -1], round, digits = 2)

table_caption <- paste(
  "Adjusted risk difference (per 100 participants) for the association",
  "between periconceptional high vegetable intake ($\\geq$ 1.25 cups per 1{,}000 kcal)",
  "and preeclampsia, estimated by residual-on-residual regression, augmented",
  "inverse probability weighting (AIPW), and targeted maximum likelihood",
  "estimation (TMLE) in the Nulliparous Pregnancy Outcomes Study: Monitoring",
  "Mothers-to-be Heart Health Study (nuMoM2b-HHS; $n = 7{,}923$). All estimators were fit with a 10-fold",
  "cross-fitted Super Learner ensemble, adjusting for sociodemographic,",
  "behavioral, medical, neighborhood, and dietary confounders.",
  "CL, confidence limit."
)

print(
  xtable(results,
         caption = table_caption,
         label   = "tab:numom_results",
         align   = c("l", "l", "r", "r", "r")),
  include.rownames     = FALSE,
  caption.placement    = "top",
  sanitize.text.function = identity,
  booktabs             = TRUE
)
