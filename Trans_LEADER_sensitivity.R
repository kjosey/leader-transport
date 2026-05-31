####-----------------------------------------------------------------------------
## Sensitivity analyses (cohort A, peer-review companion to the primary
## analysis in Trans_LEADER_primary.R).
##
##   Sensitivity 0 : cross-fitting on/off (K = 5 primary vs K = 1 in-sample).
##   Sensitivity 1 : balance on quadratic and interaction moments.
##   Sensitivity 2 : sex-stratified transportability.
##   Sensitivity 3 : balance-tolerance grid.
##   Sensitivity 4 : single-algorithm outcome regressions.
##   Sensitivity 5 : Huang (2024) variance-based sensitivity for an omitted
##                   moderator, via the senseweight package.
####-----------------------------------------------------------------------------

library(data.table)
library(WeightIt)
library(optweight)

library(SuperLearner)
library(ranger)
library(earth)
library(glmnet)

library(survival)
library(eventglm)

library(ggplot2)
library(cowplot)
library(cobalt)
library(senseweight)

source("P:/ORD_Raghavan_201905055D/Grace/transportability/LEADER_reviewer/transport_helpers.R")

result_path <- "P:/ORD_Raghavan_201905055D/Grace/transportability/LEADER_reviewer"
load(file.path(result_path, "cohort_a_setup.rda"))

## Primary balancing weights (reused by Sensitivity 0 and Sensitivity 4).
trt_wts_primary <- opt_wt_a$weights[data_a$LEADER == 1]
names(trt_wts_primary) <- data_a$SUBJID[data_a$LEADER == 1]
samp_primary <- make_samp(cohort_a, LD_data1, trt_wts_primary)

###--- Diagnostics (cohort A) -------------------------------------------------

p_love <- love.plot(opt_wt_a, binary = "std", thresholds = c(m = 0.1))
ggsave(file.path(result_path, "loveplot_cohortA.png"),
       p_love, width = 7, height = 5, dpi = 300)

trt_wts <- opt_wt_a$weights[data_a$LEADER == 1]

p_wt <- ggplot(data.frame(w = trt_wts), aes(x = w)) +
  geom_histogram(bins = 50, fill = "grey40", color = "white") +
  labs(x = "Approximate balancing weight (trial)", y = "Count") + theme_bw()

ggsave(file.path(result_path, "weight_histogram_cohortA.png"),
       p_wt, width = 6, height = 4, dpi = 300)

overlap_vars <- c("age", "A1C", "egfr")
overlap_long <- rbindlist(lapply(overlap_vars, function(v) {
  data.frame(variable = v,
             value = data_a[[v]],
             population = ifelse(data_a$LEADER == 1, "LEADER (trial)", "VA Cohort A (target)"),
             weight = ifelse(data_a$LEADER == 1, opt_wt_a$weights, 1))
  }))

p_overlap <- ggplot(overlap_long, aes(x = value, weight = weight, fill = population, color = population)) +
  geom_density(alpha = 0.35, adjust = 1.1) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(x = NULL, y = "Density (trial weighted to target)", fill = NULL, color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(result_path, "covariate_overlap_cohortA.png"),
       p_overlap, width = 9, height = 3.5, dpi = 300)

formula_1 <- as.formula(paste("LEADER ~", paste(var, collapse = "+")))
ps_fit <- glm(formula_1, data = data_a, family = binomial())
ps_long <- data.frame(ps = predict(ps_fit, type = "response"),
                      population = ifelse(data_a$LEADER == 1, "LEADER (trial)", "VA Cohort A (target)"))

p_ps <- ggplot(ps_long, aes(x = ps, fill = population, color = population)) +
  geom_density(alpha = 0.35) +
  labs(x = "Estimated P(S = 1 | X)", y = "Density", fill = NULL, color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

ggsave(file.path(result_path, "propensity_overlap_cohortA.png"),
       p_ps, width = 6, height = 4, dpi = 300)

###--- Sensitivity 0: cross-fitting on/off (K = 5 vs K = 1) -------------------

set.seed(42)

cf_k5 <- temp_m(LD_data1 = LD_data1, cohort = cohort_a,
                cohort_name = "K = 5 (cross-fit)",
                vars = var_trt, cuts = psd_cut, samp = samp_primary, K = 5L)

cf_k1 <- temp_m(LD_data1 = LD_data1, cohort = cohort_a,
                cohort_name = "K = 1 (no cross-fit)",
                vars = var_trt, cuts = psd_cut, samp = samp_primary, K = 1L)

write.csv(rbind(cf_k5, cf_k1),
          file.path(result_path, "LD_trans_rmst_rslt_cohortA_crossfit_sens.csv"),
          row.names = FALSE)

###--- Sensitivity 1: balance on interactions and squared moments ------------

sq_terms <- paste0("I(", c("age", "A1C", "egfr"), "^2)")
int_terms <- c("age:A1C", "age:egfr", "A1C:egfr",
               "age:cardiac_disease", "A1C:Med_cat")

formula_hi <- as.formula(paste("LEADER ~", paste(c(var, sq_terms, int_terms), collapse = " + ")))

opt_wt_a_hi <- weightit(formula_hi, data = data_a, focal = 0,
                        method = "optweight", estimand = "ATT", tols = 0.05)
trt_wts_hi <- opt_wt_a_hi$weights[data_a$LEADER == 1]
names(trt_wts_hi) <- data_a$SUBJID[data_a$LEADER == 1]

samp_hi <- make_samp(cohort_a, LD_data1, trt_wts_hi)
tbl_hi <- temp_m(LD_data1 = LD_data1, cohort = cohort_a,
                 cohort_name = "higher-order balance",
                 vars = var_trt, cuts = psd_cut, samp = samp_hi)

write.csv(tbl_hi, file.path(result_path, "LD_trans_rmst_rslt_cohortA_higher_order.csv"), row.names = FALSE)

###--- Sensitivity 2: sex-stratified transportability ------------------------

var_strat <- setdiff(var, "SEX")
var_trt_strat <- setdiff(var_trt, "SEX")
formula_strat <- as.formula(paste("LEADER ~", paste(var_strat, collapse = " + ")))

transport_by_sex <- function(sex_level, tag) {

  trial_s <- data1_LD[SEX == sex_level]
  cohort_s <- cohort_a[SEX == sex_level]
  LD_s <- LD_data1[SEX == sex_level]

  if (nrow(trial_s) < 50 || nrow(cohort_s) < 50 || nrow(LD_s) < 50) return(NULL)

  d_s <- rbind(trial_s, cohort_s)
  wf_s <- weightit(formula_strat, data = d_s, focal = 0,
                   method = "optweight", estimand = "ATT", tols = 0.05)
  w_s <- wf_s$weights[d_s$LEADER == 1]
  names(w_s) <- d_s$SUBJID[d_s$LEADER == 1]

  temp_m(LD_data1 = LD_s, cohort = cohort_s, cohort_name = tag,
         vars = var_trt_strat, cuts = psd_cut,
         samp = make_samp(cohort_s, LD_s, w_s))

}

write.csv(rbind(transport_by_sex("M", "cohort_a men"),
                transport_by_sex("F", "cohort_a women (exploratory)")),
          file.path(result_path, "LD_trans_cohortA_sex_stratified.csv"),
          row.names = FALSE)


###--- Sensitivity 3: balance-tolerance grid ---------------------------------

tol_results <- lapply(c(0.005, 0.01, 0.05, 0.1), function(tol) {

  wf <- weightit(formula_1, data = data_a, focal = 0,
                 method = "optweight", estimand = "ATT", tols = tol)

  w_t <- wf$weights[data_a$LEADER == 1]
  names(w_t) <- data_a$SUBJID[data_a$LEADER == 1]

  temp_m(LD_data1 = LD_data1, cohort = cohort_a,
         cohort_name = sprintf("tol = %.3f", tol),
         vars = var_trt, cuts = psd_cut,
         samp = make_samp(cohort_a, LD_data1, w_t))

})

write.csv(do.call(rbind, tol_results), file.path(result_path, "LD_trans_rmst_rslt_cohortA_tolerance_sens.csv"), row.names = FALSE)

###--- Sensitivity 4: single-algorithm outcome regressions -------------------

single_lib_grid <- list(glmnet = "SL.glmnet", earth = "SL.earth",
                        ranger = "SL.ranger", gam = "SL.gam")

lib_results <- lapply(names(single_lib_grid), function(tag) {
  temp_m(LD_data1 = LD_data1, cohort = cohort_a,
         cohort_name = paste0("outcome = ", tag),
         vars = var_trt, cuts = psd_cut, samp = samp_primary,
         sl.lib = single_lib_grid[[tag]])
})

write.csv(do.call(rbind, lib_results), file.path(result_path, "LD_trans_rmst_rslt_cohortA_outcome_library_sens.csv"), row.names = FALSE)

###--- Sensitivity 5: Huang (2024) variance-based sensitivity ----------------

w_named <- opt_wt_a$weights[data_a$LEADER == 1]
names(w_named) <- data_a$SUBJID[data_a$LEADER == 1]
w_trial <- w_named[match(LD_data1$SUBJID, names(w_named))]
Z_trial <- as.numeric(LD_data1$liraglutide)

pseudo_at <- function(evt, t) {

  o <- switch(evt,
              Composite = list(t = LD_data1$Time_comp, e = LD_data1$out_comp),
              MI = list(t = LD_data1$time_MI, e = LD_data1$MI),
              Stroke = list(t = LD_data1$time_Stroke, e = LD_data1$Stroke),
              Death = list(t = LD_data1$Time_Death, e = LD_data1$Death))

  eventglm::pseudo_infjack(
    Surv(t_, e_) ~ 1, time = t, cause = 1, type = "survival",
    data = data.frame(t_ = o$t, e_ = as.numeric(o$e))
  )

}

huang_cell <- function(evt, t) {

  Y <- pseudo_at(evt, t)
  k <- is.finite(w_trial) & is.finite(Y) & is.finite(Z_trial)
  s2 <- var(Y[k & Z_trial == 1]) + var(Y[k & Z_trial == 0])
  cbind(evt_name = evt, psd_cut = t,
        summarize_sensitivity(weights = w_trial[k], Y = Y[k], Z = Z_trial[k],
                              sigma2 = s2, estimand = "PATE"))

}

ate_rows <- subset(as.data.frame(tbl_rslt_a), value == "ATE", select = c("evt_name", "psd_cut"))
rv_grid <- do.call(rbind, Map(huang_cell, ate_rows$evt_name, ate_rows$psd_cut))
write.csv(rv_grid, file.path(result_path, "huang_rv_grid.csv"), row.names = FALSE)

## Focal estimand: composite at 36 months (benchmarking and contour) -------

Y0 <- pseudo_at("Composite", 36)
k0 <- is.finite(w_trial) & is.finite(Y0) & is.finite(Z_trial)
w0 <- w_trial[k0]
Y0 <- Y0[k0]
Z0 <- Z_trial[k0]
s2 <- var(Y0[Z0 == 1]) + var(Y0[Z0 == 0])
sens <- summarize_sensitivity(weights = w0, Y = Y0, Z = Z0,
                              sigma2 = s2, estimand = "PATE")

LD_k0 <- as.data.frame(LD_data1)[k0, var]

combined <- rbind(
  data.frame(LD_k0,
             S = 1L, Z = Z0, Y = Y0),
  data.frame(as.data.frame(cohort_a[, var, with = FALSE]),
             S = 0L, Z = NA_real_, Y = NA_real_)
)

combined[] <- lapply(combined, function(c) if (is.character(c)) factor(c) else c)

bench <- run_benchmarking(weighting_vars = var, data = combined,
                          treatment = "Z", outcome = "Y", selection = "S",
                          estimate = sens$Estimate, RV = sens$RV,
                          sigma2 = s2, estimand = "PATE")
write.csv(bench, file.path(result_path, "huang_benchmark_composite_36mo.csv"),
          row.names = FALSE)

ggsave(file.path(result_path, "huang_contour_composite_36mo.png"),
       contour_plot(var(w0), s2, sens$Estimate, bench,
                    benchmark = TRUE, shade = TRUE,
                    shade_var = c("age", "A1C"), label_size = 4),
       width = 7, height = 5, dpi = 300)