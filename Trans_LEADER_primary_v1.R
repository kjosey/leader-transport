
####-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## This file runs final models for sensitivity analyses
  
##-----------------------------------------------------------------------------------------------------------------------------------------------------------------

library(dplyr)
library(readr)
library(cobalt)
library(MatchIt)
library ("data.table")
library(table.express)
library(optmatch)
library(randomForest)
library(dbarts)
library(tableone)
library(WeightIt)
library(osqp)
library(optweight)

# Super Learner Libraries
library(SuperLearner)
library(ranger)
library(earth)
library(glmnet)
  
library(survival)
library(survminer) 
library(reshape)
library(ggplot2)
library(cowplot)
  
###-------------------------------------------------------------------------------------------------------------------------------------------------------
  #  Secondary/sensitivity analyses:
  #  a.	Running the transportability of LEADER to cohort A only for the survival probability, risk difference, and HR estimates, including sex as a balancing variable.
###---------------------------------------------------------------------------------------------------------------------------------------------------------
  
  load(file = "P:/ORD_Raghavan_201905055d/Grace/transportability/data/data1_gen.rda")
  
  data2_gen <- data1_gen[(smoker != ''),]
  cohort_a <- data2_gen[((A1C >= 7) & (age >= 50) & (cardiac_disease==1 | CKD %in% c('stage 3', 'stage 4'))),]
  
  colSums(is.na(cohort_a))
  
  load(file = "P:/ORD_Raghavan_201905055d/Grace/transportability/data/data1_LD.rda")
  data_a <- rbind(data1_LD, cohort_a)
  
  var <- c( "White" , "age" ,"smoker", "HTN", "A1C", "Med_cat", "egfr", "priorCHF" , "priorstroke", "prior_MI", "priorPCI_CABG" ,"CKD_1", 
            "cardiac_disease", "Hyperlipidemia", "LiverDisease","SEX")  ##"COPD","Afib", "Dementia", "num_med", "race_1",   "CAD"", "Hypertension", "Cancer",
  
  
  ### sampling wt -------------------------
  formula_1 <- as.formula(paste("LEADER~", paste(var, collapse="+")))
  
  opt_wt_a <- weightit(formula_1, data = data_a, focal = 0,  method = "optweight", estimand = "ATT", tols = 0.05)
  summary(opt_wt_a)
  love.plot(opt_wt_a, binary = "std", thresholds = c(m = .1))
  
  data_a[, ":="(wt_a_opt = opt_wt_a[["weights"]])]
  
  data_a1 <- data_a[LEADER==1, .(SUBJID, wt_a_opt),]
 
  
  ###load LD outcome data --------------------
  load(file = "P:/ORD_Raghavan_201905055d/Grace/transportability/data/LD_data.rda")
  LD_data1 <- LD_data[(!is.na(EGFREPB) & !is.na(HBA1CBL)) & !is.na(AGE),]  ### trial data with outcome
  
  LD_data1[CKD =='Stage 1', ":=" (CKD_1 = 'stage 1')]
  LD_data1[CKD =='Stage 2', ":=" (CKD_1 = 'stage 2')]
  LD_data1[is.na(CKD_1), ":=" (CKD_1 = 'stage 3/4') ]
  
  LD_data1[race_1 == 'WHITE', ":=" (White = 1)]
  LD_data1[is.na(White), ":="(White = 0)]
  
  LD_data1[data_a1, on = c('SUBJID'), smpl_wt_a_opt:= i.wt_a_opt]
  
  
  ### outcome model-----------------------------------
  ## RMST --------------------
  
  temp <- function(S, X, Z1, time1, event1, evt_name, cuts, samp) {
    
    event1 <- as.numeric(event1)
    
    as.data.frame(X)
    as.numeric(Z1)
    as.numeric(S)
    X1 <- subset(X, S == 1)
    
    n <- nrow(X)
    n0 <- sum(S == 0)
    n1 <- sum(S == 1)
    
    # pseudo obersvations for RMST
    lfit <- surv_fit(Surv(time, event) ~ 1, data = data.frame(event = event1, time = time1))
    psi1.mat <- as.matrix(pseudo(lfit, times = cuts, type = "surv"))
    
    pseudo.out <- sapply(1:ncol(psi1.mat), function(i, ...) {
      
      # outcome model
      psi1 <- psi1.mat[,i]
      out <- SuperLearner(Y = psi1, X = data.frame(X1, Z = Z1),
                          SL.library = c("SL.mean", "SL.glmnet", "SL.earth", "SL.gam", "SL.ranger"), family = gaussian())
      
      # prediction over every observation
      mu_tmp <- cbind(c(predict(out, newdata = data.frame(X, Z = 0))$pred),
                      c(predict(out, newdata = data.frame(X, Z = 1))$pred))
      
      
      psi <- rep(0, n)
      psi[S == 1] <- psi1
      Z <- rep(0, n)
      Z[S == 1] <- Z1
      
      # weights multipliers
      
      w0 <- (1 - Z) * S * samp / mean(1 - Z1)
      w1 <- Z * S * samp / mean(Z1)
      
      mu <- cbind(Z*mu_tmp[,2] + (1 - Z)*mu_tmp[,1], mu_tmp)
      
      # Get EIF-based Estimators
      theta0 <- sum(I(S==1) * w0 * (psi - mu[,1]))/n1 + mean(mu[S == 0,2])
      theta0_eic <- c(I(S==1) * w0 * (psi - mu[,1]))/mean(I(S == 1)) + c(I(S == 0) * mu[,2])/mean(I(S == 0)) - theta0
      
      theta1 <- sum(I(S==1) * w1 * (psi - mu[,1]))/n1 + mean(mu[S == 0,3])
      theta1_eic <- c(I(S==1) * w1 * (psi - mu[,1]))/mean(I(S == 1)) + c(I(S == 0) * mu[,3])/mean(I(S == 0)) - theta1
      
      tau <- sum(I(S==1) * (w1 - w0) * (psi - mu[,1]))/n1 + mean(mu[S == 0,3] - mu[S == 0,2])
      tau_eic <- c(I(S==1) * (w1 - w0) * (psi - mu[,1]))/mean(I(S == 1)) + c(I(S == 0) * (mu[,3] - mu[,2]))/mean(I(S == 0)) - tau
      
      # Point Estimates
      
      EY0 <- theta0
      EY0_var <- var(theta0_eic)/n
      
      EY1 <- theta1
      EY1_var <- var(theta1_eic)/n
      
      ate <- tau
      ate_var <- var(tau_eic)/n
      
      # output to next level
      return(c(ate = ate, ate_var = ate_var, EY0 = EY0, EY0_var = EY0_var, EY1 = EY1, EY1_var = EY1_var) )
      
    })
    
    ate_estimate <- pseudo.out[1,]
    ate_variance <- pseudo.out[2,]
    ate_lower <- ate_estimate - 1.96*sqrt(ate_variance)
    ate_upper <- ate_estimate + 1.96*sqrt(ate_variance)
    
    EY0_estimate <- pseudo.out[3,]
    EY0_variance <- pseudo.out[4,]
    EY0_lower <- EY0_estimate - 1.96*sqrt(EY0_variance)
    EY0_upper <- EY0_estimate + 1.96*sqrt(EY0_variance)
    
    EY1_estimate <- pseudo.out[5,]
    EY1_variance <- pseudo.out[6,]
    EY1_lower <- EY1_estimate - 1.96*sqrt(EY1_variance)
    EY1_upper <- EY1_estimate + 1.96*sqrt(EY1_variance)
    
    ## tbl <- cbind(ate_estimate, ate_variance, ate_lower, ate_upper, psd_cut, evt_name)
    
    tbl <- rbind(data.frame(value = "ATE", estimate = ate_estimate, variance = ate_variance, 
                            lower = ate_lower, upper = ate_upper, psd_cut = psd_cut, evt_name = evt_name),
                 data.frame(value = "EY0", estimate = EY0_estimate, variance = EY0_variance, 
                            lower = EY0_lower, upper = EY0_upper, psd_cut = psd_cut, evt_name = evt_name),
                 data.frame(value = "EY1", estimate = EY1_estimate, variance = EY1_variance, 
                            lower = EY1_lower, upper = EY1_upper, psd_cut = psd_cut, evt_name = evt_name))
    
    return (tbl)
    
  }
  ##--------------------------------------------------------
  LD_data1[ARM == "Liraglutide", ":=" (Liralutide = 1)]
  LD_data1[ARM == "Placebo", ":=" (Liralutide = 0)]
  
  #LD_data1[, ":="(Time_comp_day = Time_comp*(365/12))]
  
  psd_cut <- c( 6, 12, 18, 24, 30, 36, 42, 48, 54)
  ##var_trt <- c("SEX", "age")
  
  
  var_trt <- c( "White" , "age" ,"smoker", "HTN", "A1C", "Med_cat", "egfr", "priorCHF" , "priorstroke", "prior_MI", "priorPCI_CABG" ,"CKD_1", 
                "cardiac_disease", "Hyperlipidemia", "LiverDisease", "SEX")
  
  temp_m <- function(LD_data1, cohort, cohort_name,  cuts, vars, samp) {
    
    fmla <- as.formula(paste("~ ", paste(vars, collapse="+")))
    X0 <- model.frame(fmla, data = cohort)
    X1 <- model.frame(fmla, data = LD_data1)
    X <- rbind(X0, X1)
    S <- rep(c(0,1), c(nrow(X0), nrow(X1)))
    
    # Trial Specific
    Z1 <- LD_data1$Liralutide 
    
    event1_comp <- as.numeric(LD_data1$out_comp)
    time1_comp <- LD_data1$Time_comp
    event1_MI <- as.numeric(LD_data1$MI)
    time1_MI <- LD_data1$time_MI
    event1_strk <- as.numeric(LD_data1$Stroke)
    time1_strk <- LD_data1$time_Stroke
    event1_death <- as.numeric(LD_data1$Death)
    time1_death <- LD_data1$Time_Death
    
    comp <- temp(S = S, X = X, Z1 = Z1, time1 = time1_comp, event1 = event1_comp, cuts = cuts, evt_name = "Composite", samp = samp)
    MI <- temp(S = S, X = X, Z1 = Z1, time1 = time1_MI, event1 = event1_MI, cuts = cuts, evt_name = "MI" , samp = samp)
    strk <- temp(S = S, X = X, Z1 = Z1, time1 = time1_strk, event1 = event1_strk, cuts = cuts, evt_name = "Stroke" , samp = samp)
    death <- temp(S = S, X = X, Z1 = Z1, time1 = time1_death, event1 = event1_death, cuts = cuts, evt_name = "Death" , samp = samp)
    
    tbl_est <- cbind(rbind(comp, MI, strk, death), cohort_name)
    
    
    return(tbl_est)
    
  }
  
  
  tbl_rslt_a <- temp_m(LD_data1 = LD_data1, cohort = cohort_a, cohort_name = "cohort_a with sex",  vars = var_trt, cuts = psd_cut, samp = data_a$wt_a_opt)
  
  
  tbl_est_a <- print(tbl_rslt_a)
  write.csv( tbl_est_a, "P:/ORD_Raghavan_201905055D/Grace/transportability/result/CDS/LD_trans_rmst_rslt_cohortA_sex.csv")
  
  