# leader-transport

Doubly-robust **transportability** of the cardiovascular effect of liraglutide from
the **LEADER** trial to real-world **Veterans Affairs (VA)** target cohorts.

This repository holds the R code for the manuscript *"Real-world cardiovascular
effects of liraglutide: a transportability analysis of the LEADER trial"*
(under review, *American Journal of Epidemiology*, AJE-00137-2026). It has two parts:

1. the **applied transportability analysis** (LEADER → VA cohorts A–E), and
2. a **simulation study** characterizing the estimator under misspecification and
   positivity violations.

Individual-level data are **not** included — see [Data](#data).

## Methods in brief

- **Goal.** Transport the effect of liraglutide vs. placebo, estimated in the LEADER
  RCT, to VA target cohorts that progressively relax the trial's eligibility criteria.
- **Estimand.** Target-population difference in survival probability,
  `psi(t) = E[ S(t; A=1) - S(t; A=0) | S = 0 ]`, at landmark times `t`, for four
  endpoints: composite MACE, MI, stroke, and all-cause death. (`S` indexes trial vs.
  target membership; result files carry an `rmst` prefix for historical reasons.)
- **Estimator.** Augmented IPW (doubly robust), combining two nuisances:
  approximate balancing weights (`optweight`; Zubizarreta, 2015) that balance trial
  covariates to the target moments, and a pseudo-observation survival outcome model
  (`eventglm`) learned with `SuperLearner` (5-fold cross-fitting in the applied analysis).
- **Inference.** Closed-form influence-function (EIC) variance and Wald intervals.

## Repository structure

```
leader-transport/
├── transport_helpers.R          # AIPW estimator: temp(), temp_m(), make_samp()
├── Trans_LEADER_primary.R           # primary analysis, cohorts A–E    -> Results/Primary/
├── Trans_LEADER_sensitivity.R   # 6 sensitivity analyses (cohort A) -> Results/Sensitivity/
├── make_figures.R                   # real-data figures, Results/Primary/ -> Figures/
├── Simulation/
│   ├── 01_dgp.R                     # data-generating process, 8 scenarios
│   ├── 02_estimators.R             # dr_balance(): DR transportability estimator
│   ├── 03_run_simulation.R         # driver: runs scenarios, aggregates, plots
│   ├── simulation_table.csv        # bias & coverage by scenario × time
│   ├── simulation.png
│   └── simulation_weight_distribution.png
├── Results/
│   ├── Primary/                     # survival-difference results + weight diagnostics, cohorts A–E
│   └── Sensitivity/                 # cohort-A sensitivity outputs + overlap/love plots
└── Old Results/                     # superseded outputs from an earlier draft
```

Source files:
[`transport_helpers.R`](transport_helpers.R) ·
[`Trans_LEADER_primary.R`](Trans_LEADER_primary.R) ·
[`Trans_LEADER_sensitivity.R`](Trans_LEADER_sensitivity.R) ·
[`make_figures.R`](make_figures.R) ·
[`Simulation/01_dgp.R`](Simulation/01_dgp.R) ·
[`Simulation/02_estimators.R`](Simulation/02_estimators.R) ·
[`Simulation/03_run_simulation.R`](Simulation/03_run_simulation.R)

### Core functions ([`transport_helpers.R`](transport_helpers.R))

| Function | Purpose |
|---|---|
| `temp()` | AIPW transport estimator for one endpoint; pseudo-observation outcome regression with K-fold cross-fitting; returns the survival-difference ATE plus `EY0`/`EY1`, each with EIC variances and Wald CIs at every landmark time. |
| `temp_m()` | Runs `temp()` over the four endpoints (Composite, MI, Stroke, Death). |
| `make_samp()` | Aligns trial balancing weights to the stacked (target, trial) row order. |

### Applied target cohorts (A–E)

Nested VA cohorts that successively relax the LEADER inclusion criteria:

| Cohort | Definition |
|---|---|
| A | A1C ≥ 7%, age ≥ 50, and (cardiac disease or CKD stage 3–4) |
| B | drops the CVD/CKD requirement |
| C | drops the age requirement |
| D | drops the A1C requirement |
| E | full VA cohort (no inclusion criteria) |

## Requirements

R (≥ 4.1) and:

```r
install.packages(c(
  "WeightIt", "optweight", "osqp", "cobalt", "MatchIt", "optmatch",
  "SuperLearner", "ranger", "earth", "glmnet", "gam",
  "survival", "eventglm", "senseweight",
  "data.table", "dplyr", "tidyr", "readr", "table.express", "tableone",
  "ggplot2", "cowplot"
))
```

The SuperLearner libraries use `SL.glmnet`, `SL.earth`, `SL.ranger`, and `SL.gam`
(hence `glmnet`, `earth`, `ranger`, `gam`); `senseweight` for
sensitivity analysis.

## Reproducing the simulation (no data required)

```bash
Rscript Simulation/03_run_simulation.R
```

`03_run_simulation.R` sources `01_dgp.R` and `02_estimators.R`, runs all 8 scenarios
at `n_sim = 1000`, and writes [`Simulation/simulation_table.csv`](Simulation/simulation_table.csv)
plus the two figures. It sets `out_dir <- ~/Documents/LEADER/Simulation` and sources
the other scripts from there, so either place the repo at `~/Documents/LEADER` or edit
`out_dir` at the top.

Scenarios:

| Scenario (code) | Label | Construction |
|---|---|---|
| `baseline` | Baseline | linear selection & outcome |
| `omitted_selection` | Incorrect Selection Model | unmeasured `U` in selection only |
| `omitted_outcome` | Incorrect Outcome Model | `U` (and `U×A`) in outcome only |
| `omitted_both` | Exchangeability Violation | `U` in both nuisances |
| `both_nonlinear` | Nonlinear Misspecification | quadratic/interaction terms in both |
| `positivity_practical` | Practical Positivity Violation | strong selection slopes |
| `positivity_structural` | Structural Positivity Violation | trial restricted to `x1 ≤ −0.5` |
| `outcome_extrapolation` | Outcome Extrapolation | target shifted off trial support |

Outputs: [`Simulation/simulation_table.csv`](Simulation/simulation_table.csv) (bias
and 95% coverage by scenario and follow-up time), `simulation.png`, and
`simulation_weight_distribution.png`.

## Reproducing the applied analysis (restricted VA data required)

Run order matters — the sensitivity script consumes a `.rda` written by the primary:

```bash
Rscript Trans_LEADER_primary.R          # -> Results/Primary/ and cohort_a_setup.rda
Rscript Trans_LEADER_sensitivity.R  # consumes cohort_a_setup.rda -> Results/Sensitivity/
```

- [`Trans_LEADER_primary.R`](Trans_LEADER_primary.R): for each cohort A–E, fits
  `optweight` balancing weights (trial → target), then writes weight diagnostics and
  survival-difference results for the four endpoints across landmark months
  (6, 12, …, 54) to [`Results/Primary/`](Results/Primary).
- [`Trans_LEADER_sensitivity.R`](Trans_LEADER_sensitivity.R) (cohort A):
  (0) cross-fitting K = 5 vs K = 1; (1) higher-order (quadratic + interaction) balance;
  (2) sex-stratified; (3) balance-tolerance grid; (4) single-algorithm outcome learners;
  (5) Huang (2024) variance-based sensitivity for an unmeasured moderator
  (`senseweight`); plus overlap/love-plot diagnostics. Outputs go to
  [`Results/Sensitivity/`](Results/Sensitivity).

**Before you run**, the applied scripts reference a secure server via absolute paths:

1. `result_path` and the `load()`/`save()` paths target restricted VA storage
   (`P:/ORD_Raghavan_…`); edit them to your environment.
2. `Trans_LEADER_primary.R` sources `transport_helpers.R` from the working directory,
   whereas `Trans_LEADER_sensitivity.R` sources it by absolute path — adjust if needed.

### Figures

```bash
Rscript make_figures.R
```

[`make_figures.R`](make_figures.R) reads the cohort result files from
`Results/Primary/` and writes `Figures/Figure2.png`, `Figure3.png`, and `Figure4.png`:

- Figure 2 — transported survival probabilities (left) and risk differences over
  follow-up (right) for VA target cohort A, by outcome.
- Figure 3 — transported risk differences at the landmark month (`end_time`, default
  36) for cohort A, by outcome.
- Figure 4 — inclusion criteria (panel A) and transported risk differences across VA
  cohorts A–E (panel B).

The result files hold transported estimates only. The LEADER trial series is overlaid
when a file `Results/Primary/LD_trans_rmst_rslt_LEADER.csv` with the same columns and
`cohort_name == "LEADER"` is present; otherwise only the transported series is drawn.
Rows whose estimates fall outside their valid range (survival probability in `[0, 1]`,
risk difference in `[-1, 1]`) are dropped as degenerate weight solutions; a cohort with
no valid rows is excluded.

## Data

Not included. The applied scripts expect three restricted R objects:

- `data1_gen.rda` → `data1_gen` — VA target-population records (cohorts A–E are subsets)
- `data1_LD.rda` → `data1_LD` — LEADER trial records, formatted to stack with the target
- `LD_data.rda` → `LD_data` — LEADER records with endpoints (`Time_comp`/`out_comp`,
  `time_MI`/`MI`, `time_Stroke`/`Stroke`, `Time_Death`/`Death`), `ARM`, `HBA1CBL`,
  `EGFREPB`, `AGE`, `SUBJID`, `SEX`

Covariates balanced and modeled: `White`, `age`, `smoker`, `HTN`, `A1C`, `Med_cat`,
`egfr`, `priorCHF`, `priorstroke`, `prior_MI`, `priorPCI_CABG`, `CKD_1`,
`cardiac_disease`, `Hyperlipidemia`, `LiverDisease`.

LEADER individual-participant data are governed by the trial sponsor; the VA cohort is
governed by VA data-use rules (VINCI/ORD). Neither can be redistributed here.

## Citation

If you use this code, please cite the manuscript (update on publication):

> [Authors, incl. K. P. Josey and A. J. Spieker]. *Real-world cardiovascular effects of
> liraglutide: a transportability analysis of the LEADER trial.* American Journal of
> Epidemiology (under review), 2026. Manuscript AJE-00137-2026.

```bibtex
@article{leader_transport_2026,
  author  = {Josey, Kevin P. and others and Spieker, Andrew J.},
  title   = {Real-world cardiovascular effects of liraglutide: a transportability analysis of the LEADER trial},
  journal = {American Journal of Epidemiology (under review)},
  year    = {2026},
  note    = {Manuscript AJE-00137-2026}
}
```

## Contact

Questions and bug reports: open an issue on
[github.com/kjosey/leader-transport](https://github.com/kjosey/leader-transport),
or contact Kevin Josey ([@kjosey](https://github.com/kjosey)).
