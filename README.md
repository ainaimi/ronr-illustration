# Residual-on-Residual Regression as a Tool for Analyzing Observational Data

Replication code for the manuscript *Residual-on-Residual Regression as a Tool
for Analyzing Observational Data* (Naimi, Jin, Yu, Parisi, Bodnar).

The paper compares residual-on-residual (RonR) regression to AIPW and TMLE for
estimating a confounder-adjusted average treatment effect, using a simulation
study and an application to the Nulliparous Pregnancy Outcomes Study:
Monitoring Mothers-to-Be Heart Health Study (nuMoM2b-HHS) data on
periconceptional vegetable intake and preeclampsia.

## Repository layout

```
code/
  simulation_study.R               Oracle parametric-nuisance simulation (sanity check
                                   that the AIPW/TMLE/RonR implementations are correct).
  simulation_study_sl.R            Main simulation study (Super Learner nuisances,
                                   4-member GLM library). Produces sim_results_sl.csv
                                   and the figure used in the manuscript.
  simulation_study_sl_positivity.R
                                   Positivity stress test: reuses the Super Learner
                                   nuisance setup of simulation_study_sl.R but holds
                                   the effect homogeneous and sweeps a single overlap-
                                   severity dial, so AIPW/TMLE/RonR all target the same
                                   ATE while overlap degrades. Produces
                                   sim_results_sl_positivity.csv plus the
                                   positivity-sweep summaries and figures.
  numom2b_residual_on_residual.R   Application analysis on the nuMoM2b-HHS data.
  seed_gen.R                       Generates the reproducible seed bank.

data/
  random_seed_values.csv   Reproducible Monte Carlo seed bank used by the
                           simulation scripts. The nuMoM2b-HHS application data
                           are not redistributed here (see "Data access" below).

output/
  sim_summary_sl.{csv,rds}         Per-scenario performance summaries used in
                                   the manuscript text.
  sim_summary_weights_sl.{csv,rds} Per-scenario mean Super Learner ensemble
                                   weights for the propensity, conditional
                                   outcome, and marginal outcome nuisance fits.
  sim_summary_sl_positivity.{csv,rds}
                                   Positivity stress-test performance summaries, one
                                   row per (overlap severity, method), including the
                                   per-severity TMLE failure count.
  sim_summary_weights_sl_positivity.{csv,rds}
                                   Mean Super Learner ensemble weights by overlap
                                   severity for the positivity stress test.
  sim_summary_posviol_sl_positivity.{csv,rds}
                                   Realized positivity-violation diagnostics by overlap
                                   severity (true- and estimated-propensity extremes,
                                   tail mass near 0/1, largest inverse-probability
                                   weight, TMLE failure rate).
  figures/sim_performance_sl.{pdf,png}
                                   Figure 1 of the manuscript.
  figures/sim_performance_sl_positivity.{pdf,png}
                                   Positivity stress-test performance (bias, RMSE,
                                   SE ratio, 95% coverage) versus overlap severity.
  figures/sim_instability_vs_violation_sl_positivity.{pdf,png}
                                   RMSE and 95% coverage versus the realized largest
                                   inverse-probability weight.

figures/                  Figures used directly in the manuscript .tex.

manuscript/               LaTeX source for the manuscript (with appendix) and
                          its bibliography.
```

The large raw per-replicate results files (`output/sim_results.csv`,
`output/sim_results_sl.csv`, and `output/sim_results_sl_positivity.csv`) are not
committed; they are regenerable by running the corresponding scripts under
`code/`.

## Reproducing the simulation study

The main simulation results in the manuscript come from
`code/simulation_study_sl.R`. From an R session at the repository root:

```r
source("code/simulation_study_sl.R")
```

This will (i) regenerate `output/sim_results_sl.csv` from the seeded scenario
grid, (ii) re-derive `output/sim_summary_sl.{csv,rds}` and
`output/sim_summary_weights_sl.{csv,rds}`, and (iii) re-render
`output/figures/sim_performance_sl.{pdf,png}`. Wall-clock cost on an 11-core
laptop is roughly 1.5--2 hours; the script parallelises across replicates via
`future.apply::multisession`.

The oracle benchmark in `code/simulation_study.R` uses correctly-specified
parametric nuisances and provides a sanity check that the AIPW, TMLE, and RonR
implementations recover bias near zero and 95\% coverage when the nuisance
functional forms are oracle-known.

### Positivity stress test

`code/simulation_study_sl_positivity.R` re-runs the Super Learner comparison as
overlap (positivity) degrades. The treatment effect is homogeneous (a constant
`psi` for every unit, with no treatment-by-covariate interaction), so the
population ATE is exactly `psi` no matter how poor the overlap. A single
overlap-severity dial `zeta` then steepens the propensity score and pushes
`pi(C)` toward the 0/1 boundaries; `zeta = 1` reproduces the healthy overlap of
the main simulation. Because the estimand never changes, AIPW, TMLE, and RonR
all still target the same `psi`, so any divergence in their finite-sample
behaviour as `zeta` grows reflects instability rather than a change of estimand.
The propensity-truncation bound is also loosened (`pi_trunc = 0.001`, versus
`0.01` in the main simulation) so that this instability is not masked by
aggressive clipping.

```r
source("code/simulation_study_sl_positivity.R")
```

The grid sweeps `zeta` over {1, 2, 3, 4, 6} at fixed `psi = 0.5` and
`sigma_Y = 1`, again with 10,000 replicates per scenario drawn from the same
seed bank. Each replicate also records the realized degree of positivity
violation it drew -- true- and estimated-propensity extremes, the tail mass near
0/1, the largest inverse-probability weight, and a TMLE-failure flag -- so the
resulting instability can be read against an interpretable overlap metric rather
than the abstract dial. Summaries go to
`output/sim_summary_*_sl_positivity.{csv,rds}` (performance, ensemble weights,
and positivity diagnostics) and two figures to `output/figures/`. With five
overlap scenarios rather than the main study's fifteen (psi x sigma_Y) cells,
wall-clock cost is roughly a third of the main simulation.

Complete implementation details (data-generating mechanism, nuisance library,
cross-fitting scheme, estimator definitions) are also written out in the
manuscript's appendix.

### Required R packages

```
SuperLearner, tmle, sandwich, lmtest, future.apply, progressr,
here, tibble, dplyr, tidyr, ggplot2
```

## Data access

The nuMoM2b-HHS data used in the application section of the manuscript are not
redistributed in this repository. Access to the underlying study data is
governed by the nuMoM2b-HHS Data and Specimen Distribution Subcommittee; see
the study consortium website for application procedures.

`code/numom2b_residual_on_residual.R` is provided so that the analytic pipeline applied
to the application data can be inspected and replicated by authorised users
who have obtained the underlying data through the proper channels.

## Citation

Forthcoming.
