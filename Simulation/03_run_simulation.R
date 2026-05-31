###=============================================================================
### Run the simulation, aggregate, and plot
###=============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

out_dir <- file.path("~", "Documents", "LEADER", "Simulation")
out_dir <- normalizePath(out_dir, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(out_dir, "01_dgp.R"))
source(file.path(out_dir, "02_estimators.R"))

scenarios <- c("baseline",
               "omitted_selection",
               "omitted_outcome",
               "omitted_both",
               "both_nonlinear",
               "positivity_practical",
               "positivity_structural",
               "outcome_extrapolation")

n_sim <- 1000
n_total <- 1000
cuts <- c(1, 2, 3, 4, 5, 6)

set.seed(42)

run_it <- function(scenario_name) {
  
  estimates <- matrix(NA_real_, n_sim, length(cuts))
  variances <- matrix(NA_real_, n_sim, length(cuts))
  truths <- matrix(NA_real_, n_sim, length(cuts))
  weights_list <- vector("list", n_sim)
  
  for (i in seq_len(n_sim)) {
    
    sim_data <- generate_data(n_total = n_total, scenario = scenario_name,
                              cuts = cuts)
    df <- sim_data$df
    S <- df$S
    W <- data.frame(x1 = df$x1, x2 = df$x2)
    Z1 <- df$A[S == 1]
    t1 <- df$time[S == 1]
    e1 <- df$status[S == 1]
    
    fit <- tryCatch(
      dr_balance(S = S, X = W, Z1 = Z1, time1 = t1, event1 = e1,
                 cuts = cuts, method = "optweight", tols = 0.05),
      error = function(e) list(estimate = rep(NA, length(cuts)),
                               variance = rep(NA, length(cuts)),
                               wfit = NULL)
    )
    
    estimates[i, ] <- fit$estimate
    variances[i, ] <- fit$variance
    truths[i, ] <- sim_data$truth$diff
    
    if (!is.null(fit$wfit)) {
      weights_list[[i]] <- fit$wfit$weights[S == 1]
    }
    
  }
  
  # Flag extreme outliers (numerical blow-ups) as NA
  outlier <- abs(estimates) > 10 | variances > 100
  estimates[outlier] <- NA
  variances[outlier] <- NA
  
  avg_truth <- colMeans(truths, na.rm = TRUE)
  bias <- colMeans(estimates, na.rm = TRUE) - avg_truth
  rmse <- sqrt(colMeans((estimates - truths)^2, na.rm = TRUE))
  
  lo <- estimates - 1.96 * sqrt(variances)
  hi <- estimates + 1.96 * sqrt(variances)
  truth_mat <- matrix(rep(avg_truth, n_sim), byrow = TRUE, nrow = n_sim)
  coverage <- colMeans(truth_mat >= lo & truth_mat <= hi, na.rm = TRUE)
  
  summary_df <- data.frame(scenario = scenario_name, time = cuts,
                           truth = avg_truth, bias = bias, rmse = rmse,
                           coverage = coverage)
  
  rep_idx <- which(!vapply(weights_list, is.null, logical(1)))[1]
  
  weight_df <- if (length(rep_idx) && !is.na(rep_idx)) {
    data.frame(scenario = scenario_name, w = weights_list[[rep_idx]])
  } else {
    data.frame(scenario = scenario_name, w = numeric(0))
  }
  
  list(summary = summary_df, weights = weight_df)
  
}

results_by_scenario <- lapply(scenarios, run_it)
final_results <- do.call(rbind, lapply(results_by_scenario, `[[`, "summary"))
weights_long <- do.call(rbind, lapply(results_by_scenario, `[[`, "weights"))

final_results$scenario <- factor(final_results$scenario,
                                 levels = scenarios,
                                 labels = c("Baseline",
                                            "Incorrect Selection Model",
                                            "Incorrect Outcome Model",
                                            "Exchangeability Violation",
                                            "Nonlinear Misspecification",
                                            "Practical Positivity Violation",
                                            "Structural Positivity Violation",
                                            "Outcome Extrapolation"))

results_long <- final_results %>%
  select(scenario, time, bias, coverage) %>%
  pivot_longer(cols = c("bias", "coverage"),
               names_to = "metric", 
               values_to = "value") %>%
  mutate(metric = case_when(metric == "bias" ~ "Bias",
                            metric == "coverage" ~ "Coverage Probability"))

write_csv(results_long, file.path(out_dir, "simulation_table.csv"))

ref_lines <- data.frame(metric = c("Bias", "Coverage Probability"),
                        hline = c(0, 0.95))

faceted_plot <- ggplot(results_long, aes(x = time, y = value, color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.2) +
  facet_wrap(~ metric, scales = "free_y") +
  geom_hline(data = ref_lines,
             aes(yintercept = hline),
             linetype = "dashed",
             color = "black") +
  labs(title = "Estimator Performance by Scenario",
       x = "Follow-up Time (t)",
       y = "Value",
       color = "Scenario") +
  scale_x_continuous(breaks = cuts) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 21, face = "bold"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 15),
        axis.title = element_text(size = 19),
        axis.text = element_text(size = 15),
        strip.text = element_text(size = 19, face = "bold")) +
  guides(color = guide_legend(nrow = 2, override.aes = list(linewidth = 1.6, size = 4)))

ggsave(file.path(out_dir, "simulation.png"), faceted_plot,
       width = 13, height = 7.5, dpi = 300)


###=============================================================================
### Weight-density panel (Reviewer 2 major comment 4)
###=============================================================================

weights_long$scenario <- factor(weights_long$scenario,
                                levels = scenarios,
                                labels = c("Baseline",
                                           "Incorrect Selection Model",
                                           "Incorrect Outcome Model",
                                           "Exchangeability Violation",
                                           "Nonlinear Misspecification",
                                           "Practical Positivity Violation",
                                           "Structural Positivity Violation",
                                           "Outcome Extrapolation"))

weight_plot <- ggplot(subset(weights_long, w > 0),
                      aes(x = log10(w))) +
  geom_density(fill = "grey40", 
               color = "grey20", 
               alpha = 0.5) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             color = "black") +
  facet_wrap(~ scenario, 
             scales = "free_y",
             ncol = 4) +
  labs(x = expression(log[10] * "(trial weight)"),
       y = "Density") +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(size = 10, face = "bold"))

ggsave(file.path(out_dir, "simulation_weight_distribution.png"), 
       weight_plot, width = 12, height = 5, dpi = 300)
