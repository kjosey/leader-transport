####-----------------------------------------------------------------------------
## Adverse-event transportability analysis (cohort A).
##
## Binary outcomes:
##   1. Any AE leading to medication discontinuation
##   2. Composite of abdominal symptoms leading to discontinuation
##   3. Acute/chronic pancreatitis (composite)
##   4. Acute gallstone disease — composite, cholecystitis acute, cholelithiasis
##   5. Cancer leading to discontinuation
##
## Estimator: cross-fit AIPW on a single binary outcome (temp_binary()), using
## the same balancing weights and covariate set as the primary analysis.
####-----------------------------------------------------------------------------

library(data.table)
library(WeightIt)
library(optweight)

library(SuperLearner)
library(ranger)
library(earth)
library(glmnet)

library(survival)
library(ggplot2)
library(cowplot)
library(ggsci)

source("transport_helpers.R")

result_path <- "P:/ORD_Raghavan_201905055D/Grace/transportability/result/CDS"
load(file.path(result_path, "cohort_a_setup.rda"))

###--- Build trial frame with AE columns --------------------------------------
LD_adv <- fread("P:/ORD_Raghavan_201905055D/Grace/transportability/data/data_LEADER_adv.csv",
                header = TRUE, sep = ",")

LD_adv1 <- LD_adv[, .(SUBJID, ARM, adverse_discont, abdominal_symptoms,
                      adv_cancer, gall_comp, Cholelithiasis,
                      Cholecystitis_acute, pancr)]

ae_cols <- setdiff(colnames(LD_adv1), c("SUBJID", "ARM"))
setnafill(LD_adv1, cols = ae_cols, fill = 0)

LD_ae <- copy(data1_LD)
LD_ae[LD_adv1, on = c("SUBJID"), (ae_cols) := mget(ae_cols)]
LD_ae[ARM == "Liraglutide", liraglutide := 1]
LD_ae[ARM == "Placebo", liraglutide := 0]

###--- Weights aligned to (target, trial) order -------------------------------
trt_wts <- opt_wt_a$weights[data_a$LEADER == 1]
names(trt_wts) <- data_a$SUBJID[data_a$LEADER == 1]
samp <- make_samp(cohort_a, LD_ae, trt_wts)

###--- Design matrix and S indicator ------------------------------------------
fmla <- as.formula(paste("~", paste(var_trt, collapse = "+")))
X0 <- model.frame(fmla, data = cohort_a)
X1 <- model.frame(fmla, data = LD_ae)
X <- rbind(X0, X1)
S <- rep(c(0, 1), c(nrow(X0), nrow(X1)))
Z1 <- LD_ae$liraglutide

###--- Loop over AE outcomes --------------------------------------------------
ae_events <- list(
  "Adv discontinue Med" = "adverse_discont",
  "Abdominal Symptoms" = "abdominal_symptoms",
  "Cancer discontinue Med" = "adv_cancer",
  "Composite acute Gall Stone" = "gall_comp",
  "Cholelithiasis acute" = "Cholelithiasis",
  "Cholecystitis acute" = "Cholecystitis_acute",
  "Pancreatitis" = "pancr"
)

set.seed(42)
tbl_list <- lapply(names(ae_events), function(nm) {
  tryCatch(
    temp_binary(S = S, X = X, Z1 = Z1,
                Y1 = LD_ae[[ae_events[[nm]]]],
                evt_name = nm, samp = samp),
    error = function(e) NULL
  )
})

tbl_adv <- as.data.table(do.call(rbind, tbl_list))
tbl_adv[, cohort_name := "cohort_a"]

write.csv(tbl_adv[value == "ATE"],
          file.path(result_path, "LD_trans_adverse_ATE_cohortA.csv"),
          row.names = FALSE)

###--- Naive (unweighted) LEADER ATE for comparison ---------------------------
LD_adv2 <- data.frame(LD_adv1)
LD_adv2$ARM <- factor(LD_adv2$ARM, levels = c("Liraglutide", "Placebo"))

naive_rows <- lapply(ae_cols, function(col) {
  tt <- t.test(LD_adv2[[col]] ~ ARM, data = LD_adv2, var.equal = FALSE)
  data.frame(
    evt_name = col,
    mn_Liraglutide = tt$estimate[["mean in group Liraglutide"]],
    mn_placebo = tt$estimate[["mean in group Placebo"]],
    estimate = tt$estimate[["mean in group Liraglutide"]] -
               tt$estimate[["mean in group Placebo"]],
    se = tt$stderr,
    lower = tt$conf.int[1],
    upper = tt$conf.int[2]
  )
})

adv_LD_plot <- do.call(rbind, naive_rows)

save(tbl_adv, adv_LD_plot,
     file = file.path(result_path, "LD_trans_adverse_rslt_cohortA.rda"))

###--- Side-by-side forest plot (LEADER vs VA-weighted LEADER) ----------------
## House style matched to Figures 2 and 4: LEADER blue, VA-weighted orange,
## risk differences on the percent scale, clinically grouped event order.
ae_clean <- c(
  "Adv discontinue Med" = "Any adverse event\n(discontinuation)",
  "Abdominal Symptoms" = "Abdominal symptoms\n(discontinuation)",
  "Composite acute Gall Stone" = "Acute gallstone disease\n(composite)",
  "Cholecystitis acute" = "Cholecystitis",
  "Cholelithiasis acute" = "Cholelithiasis",
  "Pancreatitis" = "Pancreatitis",
  "Cancer discontinue Med" = "Cancer\n(discontinuation)"
)
ld_label <- setNames(names(ae_events), unlist(ae_events))

leader_df <- data.frame(
  estimate = adv_LD_plot$estimate,
  lower = adv_LD_plot$lower,
  upper = adv_LD_plot$upper,
  Adv_event = ld_label[adv_LD_plot$evt_name],
  cohort = "LEADER"
)

cohort_df <- as.data.frame(tbl_adv[value == "ATE",
                                   .(estimate, lower, upper,
                                     Adv_event = evt_name,
                                     cohort = "VA-weighted LEADER")])

p_adv1 <- rbind(leader_df, cohort_df)
p_adv1[c("estimate", "lower", "upper")] <-
  lapply(p_adv1[c("estimate", "lower", "upper")], function(x) 100 * as.numeric(x))
p_adv1$Adv_event <- factor(ae_clean[p_adv1$Adv_event], levels = rev(ae_clean))
p_adv1$cohort <- factor(p_adv1$cohort, levels = c("LEADER", "VA-weighted LEADER"))

src_cols <- c("LEADER" = "#2C7FB8", "VA-weighted LEADER" = "#E0A11A")
dodge <- position_dodge(width = 0.5)

p <- ggplot(p_adv1, aes(x = estimate, y = Adv_event, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, lwd = 0.3) +
  geom_errorbar(aes(xmax = upper, xmin = lower), position = dodge, width = 0) +
  geom_point(position = dodge, size = 3) +
  labs(x = "Adverse event risk difference, liraglutide - placebo (%)",
       y = NULL, color = NULL) +
  scale_color_manual(values = src_cols) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(result_path, "adverse_event_forest_cohortA.png"),
       p, width = 9, height = 6.4, dpi = 400)
ggsave(file.path(result_path, "adverse_event_forest_cohortA.tiff"),
       p, width = 9, height = 6.4, dpi = 400, compression = "lzw")
