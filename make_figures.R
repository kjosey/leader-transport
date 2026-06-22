####-----------------------------------------------------------------------------
## Final manuscript figures for the LEADER transportability analysis.
##
## One script for all three cross-population figures, each a faithful re-creation
## of the originally submitted layout (NOT the simplified make_figures.R version),
## regenerated from the FINAL tolerance-0.05 estimates incl. the now-estimable
## Cohort D:
##
##   Figure 2  : per outcome, transported survival curves (LEADER vs VA-weighted
##               LEADER, liraglutide vs placebo) beside the risk difference in
##               event-free survival over follow-up; four outcomes stacked (4x2).
##   Figure 4  : (A) UpSet-style inclusion-criteria matrix for LEADER and the VA
##               target cohorts A-E; (B) transported risk difference over
##               follow-up, one panel per outcome, LEADER and VA cohorts A-E.
##   eFigure 4 : effective sample size by cohort (bars) over the inclusion-criteria
##               matrix (UpSet layout).
##   eFigure 2 : sex-stratified transported risk differences over follow-up
##               (LEADER men/women to the male and exploratory female VA targets).
##   eFigure 3 : higher-order/interaction vs first-moment balancing, transported
##               risk differences over follow-up (Cohort A).
##
## Sources (final-data folders):
##   LEADER trial survival by arm ... Results/Trial_Results/LEADER_sp_results_by_arm.csv
##   LEADER trial RD series ......... Results/Trial_Results/LEADER_spdiff_results.csv
##   transported cohort A ........... Results/Tolerance_5_cohortA/
##   transported cohorts B-E ........ Results/Tolerance_5_cohortB_to_cohortE/
##   weight diagnostics (ESS) ....... the same Tolerance_5 folders
####-----------------------------------------------------------------------------

library(ggplot2)
library(cowplot)

base <- "~/Documents/LEADER"
res <- file.path(base, "Results")
fig_dir <- file.path(base, "Figures")
dir.create(fig_dir, showWarnings = FALSE)

t5a <- file.path(res, "Tolerance_5_cohortA")
t5be <- file.path(res, "Tolerance_5_cohortB_to_cohortE")
trial <- file.path(res, "Trial_Results")

outcome_levels <- c("Composite", "MI", "Stroke", "Death")
outcome_titles <- c(Composite = "Composite", MI = "MI", Stroke = "Stroke",
                    Death = "All-cause Mortality")

###--- Shared helpers ----------------------------------------------------------

## Transported risk differences (ATE rows). A transported risk difference must
## lie in [-1, 1]; rows outside that range signal a degenerate weight solution
## and are dropped (none occur with tolerance 0.05, including Cohort D).
read_ate <- function(path, label) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$value == "ATE", ]
  d <- d[is.finite(d$estimate) & abs(d$estimate) <= 1, ]
  data.frame(cohort = label, evt_name = d$evt_name, psd_cut = d$psd_cut,
             estimate = d$estimate, lower = d$lower, upper = d$upper,
             stringsAsFactors = FALSE)
}

## UpSet-style criteria matrix. `required` is a named list (one logical vector
## per column, aligned to crit_levels listed top-to-bottom). A filled dot marks a
## satisfied criterion; a connector runs from the topmost to the bottommost
## satisfied criterion within a column.
upset_panel <- function(crit_levels, required, x_labels, dot_size = 4,
                        base_size = 12) {
  ny <- length(crit_levels)
  cols <- names(required)
  nx <- length(cols)
  grid <- do.call(rbind, lapply(seq_len(nx), function(j) {
    data.frame(x = j, y = (ny + 1) - seq_len(ny), req = required[[cols[j]]])
  }))
  seg <- do.call(rbind, lapply(seq_len(nx), function(j) {
    ys <- (ny + 1) - which(required[[cols[j]]])
    if (length(ys) < 1) return(NULL)
    data.frame(x = j, ymin = min(ys), ymax = max(ys))
  }))
  ggplot() +
    geom_segment(data = seg, aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.6, color = "grey15") +
    geom_point(data = grid[grid$req, ], aes(x, y), size = dot_size,
               color = "grey15") +
    scale_x_continuous(breaks = seq_len(nx), labels = x_labels,
                       limits = c(0.5, nx + 0.5)) +
    scale_y_continuous(breaks = seq_len(ny), labels = rev(crit_levels),
                       limits = c(0.5, ny + 0.5)) +
    labs(x = "Cohort", y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.y = element_text(face = "bold"))
}

###--- Shared data: transported and trial risk differences ---------------------

cohort_levels <- c("LEADER", "VA A", "VA B", "VA C", "VA D", "VA E")
## LEADER in dark grey; cohorts A-E in a qualitative palette approximating the
## originally submitted figure (orange, cyan, red, green, purple).
cohort_cols <- c("LEADER" = "#3B4252", "VA A" = "#E69F00", "VA B" = "#29B6D8",
                 "VA C" = "#E15759", "VA D" = "#59A14F", "VA E" = "#8E6FB0")

trans <- rbind(
  read_ate(file.path(t5a, "LD_trans_rslt_cohort_a.csv"), "VA A"),
  read_ate(file.path(t5be, "LD_trans_rslt_cohort_b.csv"), "VA B"),
  read_ate(file.path(t5be, "LD_trans_rslt_cohort_c.csv"), "VA C"),
  read_ate(file.path(t5be, "LD_trans_rslt_cohort_d.csv"), "VA D"),
  read_ate(file.path(t5be, "LD_trans_rslt_cohort_e.csv"), "VA E")
)

ld_rd <- read.csv(file.path(trial, "LEADER_spdiff_results.csv"), stringsAsFactors = FALSE)
ld_rd <- data.frame(cohort = "LEADER", evt_name = ld_rd$Event, psd_cut = ld_rd$Month,
                    estimate = ld_rd$estimate, lower = ld_rd$lower, upper = ld_rd$upper,
                    stringsAsFactors = FALSE)

rd <- rbind(ld_rd, trans)
rd$cohort <- factor(rd$cohort, levels = cohort_levels)
rd$evt_name <- factor(rd$evt_name, levels = outcome_levels)

####=============================================================================
## Figure 2 : Cohort A survival curves and risk difference over follow-up
####=============================================================================

src_leader <- "LEADER"
src_target <- "VA-weighted LEADER"
src_levels <- c(src_leader, src_target)
src_cols <- c("#2C7FB8", "#E0A11A")
names(src_cols) <- src_levels

## Survival (EY1 = liraglutide, EY0 = placebo) for the transported Cohort A.
ca <- read.csv(file.path(t5a, "LD_trans_rslt_cohort_a.csv"), stringsAsFactors = FALSE)
ca <- ca[ca$value %in% c("EY0", "EY1"), ]
surv_target <- data.frame(
  psd_cut = ca$psd_cut, arm = ifelse(ca$value == "EY1", "Liraglutide", "Placebo"),
  source = src_target, evt_name = ca$evt_name, estimate = ca$estimate,
  stringsAsFactors = FALSE)

ls <- read.csv(file.path(trial, "LEADER_sp_results_by_arm.csv"), stringsAsFactors = FALSE)
surv_leader <- data.frame(
  psd_cut = ls$Month, arm = ls$Arm, source = src_leader, evt_name = ls$Event,
  estimate = ls$Prob_surv, stringsAsFactors = FALSE)

surv <- rbind(surv_target, surv_leader)
## anchor each curve at S(0) = 1
anchor <- expand.grid(psd_cut = 0, arm = c("Liraglutide", "Placebo"),
                      source = src_levels, evt_name = outcome_levels,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
anchor$estimate <- 1
surv <- rbind(surv, anchor)
surv$source <- factor(surv$source, levels = src_levels)
surv$evt_name <- factor(surv$evt_name, levels = outcome_levels)

## Risk difference for Cohort A: LEADER vs VA-weighted LEADER (reuse `rd`).
rd2 <- rd[rd$cohort %in% c("LEADER", "VA A"), ]
rd2$source <- factor(ifelse(rd2$cohort == "LEADER", src_leader, src_target),
                     levels = src_levels)

surv_panel <- function(evt) {
  d <- surv[surv$evt_name == evt, ]
  ggplot(d, aes(psd_cut, estimate, color = source, linetype = arm)) +
    geom_line(linewidth = 0.7) +
    geom_point(aes(shape = source), size = 1.3) +
    scale_color_manual(values = src_cols, name = "Cohort", drop = FALSE) +
    scale_linetype_manual(values = c(Liraglutide = "solid", Placebo = "dashed"),
                          name = "Arm") +
    scale_shape_manual(values = c(16, 18), guide = "none") +
    coord_cartesian(ylim = c(0.75, 1)) +
    scale_y_continuous(breaks = seq(0.75, 1, 0.05)) +
    labs(title = outcome_titles[[evt]], x = "Month", y = "Mean Survival probability") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

rd2_panel <- function(evt) {
  d <- rd2[rd2$evt_name == evt, ]
  ggplot(d, aes(psd_cut, estimate, color = source, fill = source)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.3) +
    scale_color_manual(values = src_cols, name = "Cohort", drop = FALSE) +
    scale_fill_manual(values = src_cols, name = "Cohort", drop = FALSE) +
    labs(title = outcome_titles[[evt]], x = "Month",
         y = "Difference in survival probability") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

fig2_list <- list()
for (evt in outcome_levels) {
  fig2_list[[length(fig2_list) + 1]] <- surv_panel(evt)
  fig2_list[[length(fig2_list) + 1]] <- rd2_panel(evt)
}
fig2 <- plot_grid(plotlist = fig2_list, ncol = 2)
ggsave(file.path(fig_dir, "Figure2.png"), fig2, width = 12, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "Figure2.tiff"), fig2, width = 12, height = 14,
       dpi = 300, compression = "lzw")

####=============================================================================
## Figure 4 : inclusion-criteria matrix (A) and RD over follow-up, cohorts A-E (B)
####=============================================================================

crit6 <- c("High ASCVD Risk", "Age >= 50", "A1c >= 7%", "eGFR > 15",
           "Oral med and/or long-acting insulin", "Type 2 diabetes")
req6 <- list(
  "LEADER" = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
  "VA A" = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  "VA B" = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
  "VA C" = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE),
  "VA D" = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  "VA E" = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
)
fig4a <- upset_panel(crit6, req6,
                     c("LEADER", "VA A", "VA B", "VA C", "VA D", "VA E"))

panel_rd <- function(evt) {
  d <- rd[rd$evt_name == evt, ]
  ggplot(d, aes(psd_cut, estimate, color = cohort, fill = cohort)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    scale_color_manual(values = cohort_cols, name = "Cohort", drop = FALSE) +
    scale_fill_manual(values = cohort_cols, name = "Cohort", drop = FALSE) +
    labs(title = outcome_titles[[evt]], x = "Month",
         y = "Difference in survival probability") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

fig4b <- plot_grid(plotlist = lapply(outcome_levels, panel_rd), ncol = 2)

fig4 <- plot_grid(fig4a, fig4b, ncol = 1, rel_heights = c(1, 2.2),
                  labels = c("A", "B"))
ggsave(file.path(fig_dir, "Figure4.png"), fig4, width = 9, height = 11, dpi = 300)
ggsave(file.path(fig_dir, "Figure4.tiff"), fig4, width = 9, height = 11,
       dpi = 300, compression = "lzw")

####=============================================================================
## eFigure 4 : effective sample size bars over the inclusion-criteria matrix
####=============================================================================

ess_of <- function(path) read.csv(path, stringsAsFactors = FALSE)$ess
ess <- data.frame(
  cohort = factor(c("LEADER", "A", "B", "C", "D", "E"),
                  levels = c("LEADER", "A", "B", "C", "D", "E")),
  ess = c(9336,
          ess_of(file.path(t5a, "weight_diagnostics_cohort_a.csv")),
          ess_of(file.path(t5be, "weight_diagnostics_cohort_b.csv")),
          ess_of(file.path(t5be, "weight_diagnostics_cohort_c.csv")),
          ess_of(file.path(t5be, "weight_diagnostics_cohort_d.csv")),
          ess_of(file.path(t5be, "weight_diagnostics_cohort_e.csv")))
)

p_ess <- ggplot(ess, aes(cohort, ess)) +
  geom_col(fill = "grey15", width = 0.8) +
  scale_y_continuous(limits = c(0, 10000),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Effective Sample Size") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_blank())

crit3 <- c("High ASCVD Risk", "Age >= 50", "A1c >= 7%")
req3 <- list(
  "LEADER" = c(TRUE, TRUE, TRUE),
  "A" = c(TRUE, TRUE, TRUE),
  "B" = c(TRUE, TRUE, FALSE),
  "C" = c(TRUE, FALSE, TRUE),
  "D" = c(FALSE, TRUE, TRUE),
  "E" = c(FALSE, FALSE, FALSE)
)
p_mat <- upset_panel(crit3, req3, c("LEADER", "A", "B", "C", "D", "E"))

efig4 <- plot_grid(p_ess, p_mat, ncol = 1, align = "v", axis = "lr",
                   rel_heights = c(2, 1.4))
ggsave(file.path(fig_dir, "eFigure4.png"), efig4, width = 7, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "eFigure4.tiff"), efig4, width = 7, height = 8, dpi = 300, compression = "lzw")

####=============================================================================
## Shared panel for transported risk differences over follow-up, one series per
## colour (used by eFigures 2 and 3). LEADER is carried for reference.
####=============================================================================

read_ate_series <- function(path, label, by_cohort = FALSE, recode = NULL) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$value == "ATE", ]
  d <- d[is.finite(d$estimate) & abs(d$estimate) <= 1, ]
  series <- if (by_cohort) recode[d$cohort_name] else label
  data.frame(series = series, evt_name = d$evt_name, psd_cut = d$psd_cut,
             estimate = d$estimate, lower = d$lower, upper = d$upper,
             stringsAsFactors = FALSE)
}

leader_ref <- data.frame(series = "LEADER", evt_name = ld_rd$evt_name,
                         psd_cut = ld_rd$psd_cut, estimate = ld_rd$estimate,
                         lower = ld_rd$lower, upper = ld_rd$upper,
                         stringsAsFactors = FALSE)

rd_series_panel <- function(d, evt, cols, title) {
  dd <- d[d$evt_name == evt, ]
  ggplot(dd, aes(psd_cut, estimate, color = series, fill = series)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.3) +
    scale_color_manual(values = cols, name = NULL, drop = FALSE) +
    scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
    labs(title = title, x = "Month", y = "Difference in survival probability") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom")
}

efig_grid <- function(dat, levels, cols, file) {
  dat$series <- factor(dat$series, levels = levels)
  dat$evt_name <- factor(dat$evt_name, levels = outcome_levels)
  p <- plot_grid(plotlist = lapply(outcome_levels, function(e)
    rd_series_panel(dat, e, cols, outcome_titles[[e]])), ncol = 2)
  ggsave(file.path(fig_dir, paste0(file, ".png")), p, width = 11, height = 9, dpi = 300)
  ggsave(file.path(fig_dir, paste0(file, ".tiff")), p, width = 11, height = 9,
         dpi = 300, compression = "lzw")
  p
}

####=============================================================================
## eFigure 2 : sex-stratified transported risk differences over follow-up.
##   LEADER men transported to the male VA target and LEADER women to the small
##   (exploratory) female VA subgroup; LEADER overall shown for reference.
####=============================================================================

sex_recode <- c("cohort_a men" = "Men",
                "cohort_a women (exploratory)" = "Women (exploratory)")
sex_dat <- rbind(
  leader_ref,
  read_ate_series(file.path(t5a, "LD_trans_cohortA_sex_stratified.csv"),
                  by_cohort = TRUE, recode = sex_recode)
)
sex_levels <- c("LEADER", "Men", "Women (exploratory)")
sex_cols <- c("LEADER" = "#3B4252", "Men" = "#2C7FB8",
              "Women (exploratory)" = "#E0A11A")
efig2 <- efig_grid(sex_dat, sex_levels, sex_cols, "eFigure2")

####=============================================================================
## eFigure 3 : higher-order / interaction balancing vs first-moment balancing.
##   Transported risk differences over follow-up (Cohort A); LEADER for reference.
####=============================================================================

ho_dat <- rbind(
  leader_ref,
  read_ate_series(file.path(t5a, "LD_trans_rslt_cohort_a.csv"), "First-moment balance"),
  read_ate_series(file.path(t5a, "LD_trans_rmst_rslt_cohortA_higher_order.csv"),
                  "Higher-order balance")
)
ho_levels <- c("LEADER", "First-moment balance", "Higher-order balance")
ho_cols <- c("LEADER" = "#3B4252", "First-moment balance" = "#2C7FB8",
             "Higher-order balance" = "#E0A11A")
efig3 <- efig_grid(ho_dat, ho_levels, ho_cols, "eFigure3")
