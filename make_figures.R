####-----------------------------------------------------------------------------
## Real-data figures for the LEADER transportability analysis.
##
## Reads the primary result files written by Trans_LEADER_primary.R and builds:
##   Figure 2 : transported survival curves (left) and risk differences over
##              follow-up (right), VA target cohort A.
##   Figure 3 : risk differences at the landmark time, VA target cohort A.
##   Figure 4 : inclusion criteria (A) and transported effects across the
##              VA target cohorts A-E (B).
####-----------------------------------------------------------------------------

library(ggplot2)
library(cowplot)

###--- Configuration ----------------------------------------------------------

results_dir <- "~/Documents/LEADER/Results/Primary"
fig_dir <- "Figures"
leader_file <- file.path(results_dir, "LD_trans_rmst_rslt_LEADER.csv")

end_time <- 36

cohort_levels <- c("cohort_a", "cohort_b", "cohort_c", "cohort_d", "cohort_e")
cohort_labels <- c(cohort_a = "A", cohort_b = "B", cohort_c = "C",
                   cohort_d = "D", cohort_e = "E")

outcome_levels <- c("Composite", "MI", "Stroke", "Death")
outcome_labels <- c(Composite = "Composite MACE", MI = "Non-fatal MI",
                    Stroke = "Non-fatal stroke", Death = "All-cause mortality")

src_leader <- "LEADER trial"
src_target <- "VA target (transported)"
src_cols <- c("#2C7FB8", "#E0A11A")
names(src_cols) <- c(src_leader, src_target)

dir.create(fig_dir, showWarnings = FALSE)

###--- Load results -----------------------------------------------------------

valid_rows <- function(d) {
  finite <- is.finite(d$estimate) & is.finite(d$lower) & is.finite(d$upper)
  surv <- d$value %in% c("EY0", "EY1")
  ate <- d$value == "ATE"
  finite & ((surv & d$estimate >= 0 & d$estimate <= 1) |
            (ate & abs(d$estimate) <= 1))
}

read_result <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  keep <- valid_rows(d)
  if (any(!keep)) {
    message(sprintf("  %s: dropped %d of %d rows with out-of-range estimates ",
                    basename(path), sum(!keep), nrow(d)),
            "(degenerate weights).")
  }
  d[keep, , drop = FALSE]
}

primary <- do.call(rbind, lapply(cohort_levels, function(cn) {
  path <- file.path(results_dir, sprintf("LD_trans_rmst_rslt_%s.csv", cn))
  if (!file.exists(path)) return(NULL)
  read_result(path)
}))

dropped <- setdiff(cohort_levels, unique(primary$cohort_name))
if (length(dropped)) {
  message("No finite estimates for: ", paste(dropped, collapse = ", "),
          " (excluded from the figures).")
}

leader <- if (file.exists(leader_file)) read_result(leader_file) else NULL
if (is.null(leader)) {
  message("LEADER trial estimates not found at ", leader_file,
          " - drawing the transported series only. Add a file with the same ",
          "columns and cohort_name == \"LEADER\" to overlay the trial series.")
}

## Label helper applied to every plotting frame.
prep <- function(d, source_label) {
  d$source <- source_label
  d$evt_name <- factor(d$evt_name, levels = outcome_levels)
  d
}

###--- Figure 2: survival curves (left) and RD over follow-up (right) ----------
## VA target cohort A; LEADER overlaid when available.

target_a <- prep(primary[primary$cohort_name == "cohort_a", ], src_target)
leader_a <- if (!is.null(leader)) prep(leader, src_leader) else NULL

surv <- rbind(target_a, leader_a)
surv <- surv[surv$value %in% c("EY0", "EY1"), ]
surv$arm <- ifelse(surv$value == "EY1", "Liraglutide", "Placebo")

## Anchor the curves at S(0) = 1.
anchor <- expand.grid(psd_cut = 0, arm = c("Liraglutide", "Placebo"),
                      source = unique(surv$source),
                      evt_name = factor(outcome_levels, levels = outcome_levels),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
anchor$estimate <- 1
surv <- rbind(surv[, c("psd_cut", "arm", "source", "evt_name", "estimate")], anchor)

p_surv <- ggplot(surv, aes(psd_cut, estimate, color = source, linetype = arm)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ evt_name, labeller = as_labeller(outcome_labels)) +
  scale_color_manual(values = src_cols, name = NULL) +
  scale_linetype_manual(values = c(Liraglutide = "solid", Placebo = "dashed"),
                        name = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(x = "Months since randomization", y = "Survival probability") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

rd <- rbind(target_a, leader_a)
rd <- rd[rd$value == "ATE", ]

p_rd <- ggplot(rd, aes(psd_cut, estimate, color = source, fill = source)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ evt_name, labeller = as_labeller(outcome_labels)) +
  scale_color_manual(values = src_cols, name = NULL) +
  scale_fill_manual(values = src_cols, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  labs(x = "Months since randomization",
       y = "Risk difference in event-free survival\n(liraglutide - placebo)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

fig2 <- plot_grid(p_surv, p_rd, nrow = 1, labels = c("", ""))
ggsave(file.path(fig_dir, "Figure2.png"), fig2, width = 14, height = 6, dpi = 300)

###--- Figure 3: risk differences at the landmark time (cohort A) --------------

fig3_df <- rd[rd$psd_cut == end_time, ]

p_fig3 <- ggplot(fig3_df,
                 aes(estimate, evt_name, color = source)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = src_cols, name = NULL) +
  scale_y_discrete(limits = rev(outcome_levels), labels = rev(outcome_labels)) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x, 1), "%")) +
  labs(x = sprintf("Risk difference at %d months (liraglutide - placebo)", end_time),
       y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Figure3.png"), p_fig3, width = 8, height = 5, dpi = 300)

###--- Figure 4A: inclusion criteria across target populations -----------------
## Criteria encoded from the cohort definitions in Trans_LEADER_primary.R.

crit_levels <- c("HbA1c >= 7%", "Age >= 50", "CVD or CKD 3-4")
criteria <- data.frame(
  pop = rep(c("LEADER", "A", "B", "C", "D", "E"), each = 3),
  criterion = rep(crit_levels, times = 6),
  required = c(TRUE,  TRUE,  TRUE,    # LEADER
               TRUE,  TRUE,  TRUE,    # A
               TRUE,  TRUE,  FALSE,   # B (drop CVD/CKD)
               TRUE,  FALSE, TRUE,    # C (drop age)
               FALSE, TRUE,  TRUE,    # D (drop A1C)
               FALSE, FALSE, FALSE),  # E (no criteria)
  stringsAsFactors = FALSE
)
criteria$pop <- factor(criteria$pop, levels = rev(c("LEADER", "A", "B", "C", "D", "E")))
criteria$criterion <- factor(criteria$criterion, levels = crit_levels)

p_fig4a <- ggplot(criteria, aes(criterion, pop)) +
  geom_tile(aes(fill = required), color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(required, "✓", "")), size = 5) +
  scale_fill_manual(values = c("TRUE" = "#BDD7E7", "FALSE" = "grey92"),
                    guide = "none") +
  labs(x = NULL, y = "Target population", title = "Inclusion criteria") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid = element_blank())

###--- Figure 4B: transported effects across target cohorts --------------------

across <- prep(primary[primary$value == "ATE" & primary$psd_cut == end_time, ],
               src_target)
across$pop <- factor(cohort_labels[across$cohort_name],
                     levels = rev(c("A", "B", "C", "D", "E")))

p_fig4b <- ggplot(across, aes(estimate, pop)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2,
                 color = src_cols[[src_target]]) +
  geom_point(size = 2.4, color = src_cols[[src_target]]) +
  facet_wrap(~ evt_name, nrow = 1, labeller = as_labeller(outcome_labels)) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(x = sprintf("Transported risk difference at %d months", end_time),
       y = "VA target population",
       title = "Transported treatment effect") +
  theme_bw(base_size = 12)

fig4 <- plot_grid(p_fig4a, p_fig4b, nrow = 1, rel_widths = c(1, 2.4),
                  labels = c("A", "B"))
ggsave(file.path(fig_dir, "Figure4.png"), fig4, width = 15, height = 5, dpi = 300)
