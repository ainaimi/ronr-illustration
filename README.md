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
  figures/sim_performance_sl.{pdf,png}
                                   Figure 1 of the manuscript.

figures/                  Figures used directly in the manuscript .tex.

manuscript/               LaTeX source for the manuscript (with appendix) and
                          its bibliography.
```

The large raw per-replicate results files (`output/sim_results.csv` and
`output/sim_results_sl.csv`) are not committed; they are regenerable by running
the corresponding scripts under `code/`.

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
