####-----------------------------------------------------------------------------
## Shared helpers for the LEADER transportability analysis.
##
##   temp()         AIPW transport estimator; cross-fit outcome regression
##   temp_m()       multi-outcome wrapper (composite, MI, stroke, death)
##   temp_binary()  cross-fit AIPW for a single binary outcome (adverse events)
##   make_samp()    align trial-side weights to the (target, trial) row order
##                  expected by temp_m(), temp(), and temp_binary()
##
####-----------------------------------------------------------------------------

temp <- function(S, X, Z1, time1, event1, evt_name, cuts, samp, K = 5L,
                 sl.lib = c("SL.mean", "SL.glmnet", "SL.earth",
                            "SL.gam", "SL.ranger")) {
  
  event1 <- as.numeric(event1)
  n <- nrow(X)
  n0 <- sum(S == 0)
  n1 <- sum(S == 1)
  idx_S1 <- which(S == 1)
  idx_S0 <- which(S == 0)
  X1 <- X[idx_S1, , drop = FALSE]
  X0 <- X[idx_S0, , drop = FALSE]
  
  ps_dat <- data.frame(t_ = time1, e_ = event1)
  
  psi1.mat <- vapply(cuts, function(tt) {
    eventglm::pseudo_infjack(
      formula = Surv(t_, e_) ~ 1,
      time = tt,
      cause = 1,
      type = "survival",
      data = ps_dat
    )
  }, numeric(length(time1)))
  
  fold_id <- integer(n1)
  
  for (z_val in c(0, 1)) {
    z_idx <- which(Z1 == z_val)
    fold_id[z_idx] <- sample(rep(seq_len(K), length.out = length(z_idx)))
  }
  
  pseudo.out <- sapply(seq_len(ncol(psi1.mat)), function(i) {
    
    psi1 <- psi1.mat[, i]
    
    mu0_S1 <- numeric(n1)
    mu1_S1 <- numeric(n1)
    mu0_S0_mat <- matrix(0, n0, K)
    mu1_S0_mat <- matrix(0, n0, K)
    
    for (k in seq_len(K)) {
      
      if (K == 1L) {
        tr <- te <- seq_len(n1)
      } else {
        tr <- which(fold_id != k)
        te <- which(fold_id == k)
      }
      
      fit_k <- SuperLearner(Y = psi1[tr],
                            X = data.frame(X1[tr, , drop = FALSE], Z = Z1[tr]),
                            SL.library = sl.lib, family = gaussian())
      
      mu0_S1[te] <- predict(fit_k, newdata = data.frame(X1[te, , drop = FALSE], Z = 0))$pred
      mu1_S1[te] <- predict(fit_k, newdata = data.frame(X1[te, , drop = FALSE], Z = 1))$pred
      
      mu0_S0_mat[, k] <- predict(fit_k, newdata = data.frame(X0, Z = 0))$pred
      mu1_S0_mat[, k] <- predict(fit_k, newdata = data.frame(X0, Z = 1))$pred
      
    }
    
    mu0_S0 <- rowMeans(mu0_S0_mat)
    mu1_S0 <- rowMeans(mu1_S0_mat)
    
    psi <- numeric(n)
    psi[idx_S1] <- psi1
    Z <- numeric(n)
    Z[idx_S1] <- Z1
    
    mu_obs <- numeric(n)
    mu_obs[idx_S1] <- Z1 * mu1_S1 + (1 - Z1) * mu0_S1
    mu_0 <- numeric(n)
    mu_0[idx_S1] <- mu0_S1
    mu_0[idx_S0] <- mu0_S0
    mu_1 <- numeric(n)
    mu_1[idx_S1] <- mu1_S1
    mu_1[idx_S0] <- mu1_S0
    
    w0 <- (1 - Z) * S * samp / mean(1 - Z1)
    w1 <- Z * S * samp / mean(Z1)
    
    theta0 <- sum(S * w0 * (psi - mu_obs)) / n1 + mean(mu_0[idx_S0])
    theta0_eic <- (S * w0 * (psi - mu_obs)) / mean(S == 1) +
      ((S == 0) * mu_0) / mean(S == 0) - theta0
    
    theta1 <- sum(S * w1 * (psi - mu_obs)) / n1 + mean(mu_1[idx_S0])
    theta1_eic <- (S * w1 * (psi - mu_obs)) / mean(S == 1) +
      ((S == 0) * mu_1) / mean(S == 0) - theta1
    
    tau <- sum(S * (w1 - w0) * (psi - mu_obs)) / n1 + mean(mu_1[idx_S0] - mu_0[idx_S0])
    tau_eic <- (S * (w1 - w0) * (psi - mu_obs)) / mean(S == 1) +
      ((S == 0) * (mu_1 - mu_0)) / mean(S == 0) - tau
    
    c(ate = tau, ate_var = var(tau_eic) / n,
      EY0 = theta0, EY0_var = var(theta0_eic) / n,
      EY1 = theta1, EY1_var = var(theta1_eic) / n)
    
  })
  
  ate_est <- pseudo.out[1, ]; ate_var <- pseudo.out[2, ]
  EY0_est <- pseudo.out[3, ]; EY0_var <- pseudo.out[4, ]
  EY1_est <- pseudo.out[5, ]; EY1_var <- pseudo.out[6, ]
  
  rbind(
    data.frame(value = "ATE", estimate = ate_est, variance = ate_var,
               lower = ate_est - 1.96 * sqrt(ate_var),
               upper = ate_est + 1.96 * sqrt(ate_var),
               psd_cut = cuts, evt_name = evt_name),
    data.frame(value = "EY0", estimate = EY0_est, variance = EY0_var,
               lower = EY0_est - 1.96 * sqrt(EY0_var),
               upper = EY0_est + 1.96 * sqrt(EY0_var),
               psd_cut = cuts, evt_name = evt_name),
    data.frame(value = "EY1", estimate = EY1_est, variance = EY1_var,
               lower = EY1_est - 1.96 * sqrt(EY1_var),
               upper = EY1_est + 1.96 * sqrt(EY1_var),
               psd_cut = cuts, evt_name = evt_name)
  )
  
}

temp_m <- function(LD_data1, cohort, cohort_name, cuts, vars, samp, K = 5L,
                   sl.lib = c("SL.mean", "SL.glmnet", "SL.earth",
                              "SL.gam", "SL.ranger")) {
  
  fmla <- as.formula(paste("~", paste(vars, collapse = "+")))
  X0 <- model.frame(fmla, data = cohort)
  X1 <- model.frame(fmla, data = LD_data1)
  X <- rbind(X0, X1)
  S <- rep(c(0, 1), c(nrow(X0), nrow(X1)))
  
  Z1 <- LD_data1$liraglutide
  
  outcomes <- list(
    Composite = list(t = LD_data1$Time_comp, e = LD_data1$out_comp),
    MI = list(t = LD_data1$time_MI, e = LD_data1$MI),
    Stroke = list(t = LD_data1$time_Stroke, e = LD_data1$Stroke),
    Death = list(t = LD_data1$Time_Death, e = LD_data1$Death)
  )
  
  tbl_list <- lapply(names(outcomes), function(nm) {
    o <- outcomes[[nm]]
    temp(S = S, X = X, Z1 = Z1, time1 = o$t, event1 = as.numeric(o$e),
         cuts = cuts, evt_name = nm, samp = samp, K = K, sl.lib = sl.lib)
  })
  
  cbind(do.call(rbind, tbl_list), cohort_name)
  
}

make_samp <- function(cohort_df, trial_df, trial_wts_by_subjid) {
  w_trial <- trial_wts_by_subjid[match(trial_df$SUBJID, names(trial_wts_by_subjid))]
  c(rep(1, nrow(cohort_df)), w_trial)
}

temp_binary <- function(S, X, Z1, Y1, evt_name, samp, K = 5L,
                        sl.lib = c("SL.mean", "SL.glmnet", "SL.earth",
                                   "SL.gam", "SL.ranger")) {

  Y1 <- as.numeric(Y1)
  n <- nrow(X)
  n0 <- sum(S == 0)
  n1 <- sum(S == 1)
  idx_S1 <- which(S == 1)
  idx_S0 <- which(S == 0)
  X1 <- X[idx_S1, , drop = FALSE]
  X0 <- X[idx_S0, , drop = FALSE]

  fold_id <- integer(n1)
  for (z_val in c(0, 1)) {
    z_idx <- which(Z1 == z_val)
    fold_id[z_idx] <- sample(rep(seq_len(K), length.out = length(z_idx)))
  }

  mu0_S1 <- numeric(n1)
  mu1_S1 <- numeric(n1)
  mu0_S0_mat <- matrix(0, n0, K)
  mu1_S0_mat <- matrix(0, n0, K)

  for (k in seq_len(K)) {

    if (K == 1L) {
      tr <- te <- seq_len(n1)
    } else {
      tr <- which(fold_id != k)
      te <- which(fold_id == k)
    }

    fit_k <- SuperLearner(Y = Y1[tr],
                          X = data.frame(X1[tr, , drop = FALSE], Z = Z1[tr]),
                          SL.library = sl.lib, family = binomial())

    mu0_S1[te] <- predict(fit_k, newdata = data.frame(X1[te, , drop = FALSE], Z = 0))$pred
    mu1_S1[te] <- predict(fit_k, newdata = data.frame(X1[te, , drop = FALSE], Z = 1))$pred
    mu0_S0_mat[, k] <- predict(fit_k, newdata = data.frame(X0, Z = 0))$pred
    mu1_S0_mat[, k] <- predict(fit_k, newdata = data.frame(X0, Z = 1))$pred

  }

  mu0_S0 <- rowMeans(mu0_S0_mat)
  mu1_S0 <- rowMeans(mu1_S0_mat)

  psi <- numeric(n)
  psi[idx_S1] <- Y1
  Z <- numeric(n)
  Z[idx_S1] <- Z1

  mu_obs <- numeric(n)
  mu_obs[idx_S1] <- Z1 * mu1_S1 + (1 - Z1) * mu0_S1
  mu_0 <- numeric(n)
  mu_0[idx_S1] <- mu0_S1
  mu_0[idx_S0] <- mu0_S0
  mu_1 <- numeric(n)
  mu_1[idx_S1] <- mu1_S1
  mu_1[idx_S0] <- mu1_S0

  w0 <- (1 - Z) * S * samp / mean(1 - Z1)
  w1 <- Z * S * samp / mean(Z1)

  theta0 <- sum(S * w0 * (psi - mu_obs)) / n1 + mean(mu_0[idx_S0])
  theta0_eic <- (S * w0 * (psi - mu_obs)) / mean(S == 1) +
    ((S == 0) * mu_0) / mean(S == 0) - theta0

  theta1 <- sum(S * w1 * (psi - mu_obs)) / n1 + mean(mu_1[idx_S0])
  theta1_eic <- (S * w1 * (psi - mu_obs)) / mean(S == 1) +
    ((S == 0) * mu_1) / mean(S == 0) - theta1

  tau <- sum(S * (w1 - w0) * (psi - mu_obs)) / n1 + mean(mu_1[idx_S0] - mu_0[idx_S0])
  tau_eic <- (S * (w1 - w0) * (psi - mu_obs)) / mean(S == 1) +
    ((S == 0) * (mu_1 - mu_0)) / mean(S == 0) - tau

  rbind(
    data.frame(value = "ATE", estimate = tau, variance = var(tau_eic) / n,
               lower = tau - 1.96 * sqrt(var(tau_eic) / n),
               upper = tau + 1.96 * sqrt(var(tau_eic) / n),
               evt_name = evt_name),
    data.frame(value = "EY0", estimate = theta0, variance = var(theta0_eic) / n,
               lower = theta0 - 1.96 * sqrt(var(theta0_eic) / n),
               upper = theta0 + 1.96 * sqrt(var(theta0_eic) / n),
               evt_name = evt_name),
    data.frame(value = "EY1", estimate = theta1, variance = var(theta1_eic) / n,
               lower = theta1 - 1.96 * sqrt(var(theta1_eic) / n),
               upper = theta1 + 1.96 * sqrt(var(theta1_eic) / n),
               evt_name = evt_name)
  )

}