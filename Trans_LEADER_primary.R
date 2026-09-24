####-----------------------------------------------------------------------------
## Primary transportability analysis: cohorts A through E.
##
## Cohort definitions:
##   A: A1C >= 7 AND age >= 50 AND (cardiac_disease | CKD stage 3-4)
##   B: drop A1C requirement
##   C: drop age requirement
##   D: drop CVD/CKD requirement
##   E: drop all inclusion criteria
####-----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(data.table)
library(table.express)
library(optmatch)
library(tableone)
library(WeightIt)
library(MatchIt)
library(cobalt)
library(osqp)
library(optweight)

library(SuperLearner)
library(ranger)
library(earth)
library(glmnet)

library(survival)
library(eventglm)

source("transport_helpers.R")

result_path <- "P:/ORD_Raghavan_201905055D/Grace/transportability/result/CDS"

###--- Data --------------------------------------------------------------------
load(file = "P:/ORD_Raghavan_201905055d/Grace/transportability/data/data1_gen.rda")
data2_gen <- data1_gen[smoker != "", ]

cohort_a <- data2_gen[A1C >= 7 & age >= 50 & (cardiac_disease == 1 | CKD %in% c("stage 3", "stage 4")), ]
cohort_b <- data2_gen[age >= 50 & (cardiac_disease == 1 | CKD %in% c("stage 3", "stage 4")), ]
cohort_c <- data2_gen[A1C >= 7 & (cardiac_disease == 1 | CKD %in% c("stage 3", "stage 4")), ]
cohort_d <- data2_gen[A1C >= 7 & age >= 50, ]
cohort_e <- data2_gen

load(file = "P:/ORD_Raghavan_201905055d/Grace/transportability/data/data1_LD.rda")
load(file = "P:/ORD_Raghavan_201905055d/Grace/transportability/data/LD_data.rda")

LD_data1 <- LD_data[(!is.na(EGFREPB) & !is.na(HBA1CBL)) & !is.na(AGE), ]
LD_data1[CKD == "Stage 1", CKD_1 := "stage 1"]
LD_data1[CKD == "Stage 2", CKD_1 := "stage 2"]
LD_data1[is.na(CKD_1), CKD_1 := "stage 3/4"]
LD_data1[race_1 == "WHITE", White := 1]
LD_data1[is.na(White), White := 0]
LD_data1[ARM == "Liraglutide", liraglutide := 1]
LD_data1[ARM == "Placebo", liraglutide := 0]

###--- Balancing and outcome-model covariates ---------------------------------
var <- c("White", "age", "smoker", "HTN", "A1C", "Med_cat", "egfr",
         "priorCHF", "priorstroke", "prior_MI", "priorPCI_CABG", "CKD_1",
         "cardiac_disease", "Hyperlipidemia", "LiverDisease")
var_trt <- var

psd_cut <- c(6, 12, 18, 24, 30, 36, 42, 48, 54)

###--- Cohort-level pipeline --------------------------------------------------
run_primary <- function(cohort, cohort_name, seed = 42) {

  set.seed(seed)
  data_combined <- rbind(data1_LD, cohort)
  formula_1 <- as.formula(paste("LEADER ~", paste(var, collapse = "+")))

  opt_wt <- weightit(formula_1, data = data_combined, focal = 0,
                     method = "optweight", estimand = "ATT", tols = 0.05)

  trt_wts <- opt_wt$weights[data_combined$LEADER == 1]
  names(trt_wts) <- data_combined$SUBJID[data_combined$LEADER == 1]

  wt_diag <- data.frame(
    cohort = cohort_name,
    n_trial = length(trt_wts),
    ess = sum(trt_wts)^2 / sum(trt_wts^2),
    mean = mean(trt_wts),
    median = median(trt_wts),
    sd = sd(trt_wts),
    min = min(trt_wts),
    max = max(trt_wts),
    p99 = quantile(trt_wts, 0.99, names = FALSE)
  )

  write.csv(wt_diag, file.path(result_path, sprintf("weight_diagnostics_%s.csv", cohort_name)), row.names = FALSE)

  samp <- make_samp(cohort, LD_data1, trt_wts)
  tbl_result <- temp_m(LD_data1 = LD_data1, cohort = cohort,
                       cohort_name = cohort_name,
                       vars = var_trt, cuts = psd_cut, samp = samp)

  write.csv(tbl_result, file.path(result_path, sprintf("LD_trans_rmst_rslt_%s.csv", cohort_name)), row.names = FALSE)

  list(data = data_combined, opt_wt = opt_wt, tbl_rslt = tbl_result)

}

###--- Run cohorts A through E ------------------------------------------------
res_a <- run_primary(cohort_a, "cohort_a")

## For sensitivity-analysis companion script
data_a <- res_a$data
opt_wt_a <- res_a$opt_wt
tbl_rslt_a <- res_a$tbl_rslt

save(data_a, opt_wt_a, tbl_rslt_a, LD_data1, cohort_a, data1_LD, var, var_trt, psd_cut,
     file = file.path(result_path, "cohort_a_setup.rda"))

res_b <- run_primary(cohort_b, "cohort_b")
res_c <- run_primary(cohort_c, "cohort_c")
res_d <- run_primary(cohort_d, "cohort_d")
res_e <- run_primary(cohort_e, "cohort_e")

results <- list(cohort_a = res_a$tbl_rslt,
                cohort_b = res_b$tbl_rslt,
                cohort_c = res_c$tbl_rslt,
                cohort_d = res_d$tbl_rslt,
                cohort_e = res_e$tbl_rslt)

write.csv(do.call(rbind, results), file.path(result_path, "LD_trans_rmst_rslt_cohorts_AE.csv"), row.names = FALSE)