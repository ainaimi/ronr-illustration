---
project: Residual_on_Residual_Illustration
type: paper
status: in-progress
owner: ain
collaborators: [Qianhui Jin, Lisa Bodnar]
journal: American Journal of Epidemiology (Practice of Epidemiology)
deadline: 
last_updated: 2026-06-18
---

## Next Actions
- [ ] Decide on `fit_nu` form (gC basis vs + sin(πC1) vs tower-property identity)
- [ ] Run scenario-grid pilot (100 reps × 15 scenarios) before committing to full sweep
- [ ] Update summary block to `group_by(method, psi_true, sigma_Y)` and use `psi_true` instead of removed global
- [ ] Manuscript TODOs: AIPW closed-form expression, two `CITE` placeholders, `[Qianhui, what does v_dich represent?]`, abstract
- [ ] Submit to AJE Practice of Epidemiology

## Notes
Applied paper demonstrating residual-on-residual regression on the numom+hhs diet data (`data/numom+hhs_diet_unimputed_v2.dta`). Manuscript now built in LaTeX (`manuscript/2026_06_15-ResidualOnResidual.tex/.pdf`) with proof sketch (`ronr_proof_sketch.tex`); earlier `.docx` drafts retained for history (QJ-suffixed versions are prior revisions from Qianhui).

Simulation study (`code/simulation_study.R`) compares OLS / RonR / AIPW / TMLE across a 5 × 3 grid of (psi, sigma_Y) using oracle parametric nuisances — SuperLearner machinery was removed in the 2026-06-18 session. Results land in `output/sim_results.csv`. Detailed session context in `_notes/`.

Literature corpus in `_corpus/` (Robinson 1988, Chernozhukov 2018, four Bodnar et al. applied papers).
