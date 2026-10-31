## Util functions

library(devtools)
library(tidyverse)
library(pbapply)
# devtools::install_github("LauBok/PStrata")
library(PStrata) 

# Sample -----------------------------
## Calculate stratum probability
calc_stratum_probability <- function(strata_coef_matrix, data_matrix) {
  t(apply(
    exp(data_matrix %*% strata_coef_matrix), 1, function(x) x / sum(x)
  ))
}

## Calculate survival probability
calc_surv_probability_AFT <- function(outcome_coef_matrix, sigma_vector, data_matrix, time_sequence) {
  mu <- data_matrix %*% outcome_coef_matrix
  lapply(1:ncol(mu), function(j) pnorm((rep(1, nrow(mu)) %*% t(log(time_sequence)) - mu[, j]) / sigma_vector[j], lower.tail = FALSE))
}

calc_surv_probability <- function(outcome_coef_matrix, phi_vector, data_matrix, time_sequence) {
  exp.mu <- exp(data_matrix %*% outcome_coef_matrix)
  lapply(1:ncol(exp.mu), function(j) exp(- (exp.mu[, j]) %*% t(time_sequence^phi_vector[j]) / phi_vector[j]))
}

calc_RACE <- function(outcome_coef_matrix, phi_vector, data_matrix, time_sequence) {
  exp.mu <- exp(data_matrix %*% outcome_coef_matrix)
  c = t(t(exp.mu) / phi_vector)
  lapply(1:ncol(exp.mu), function(j) (phi_vector[j] * c[, j]^(1/phi_vector[j]))^(-1) * pgamma(c[, j] %*% t(time_sequence^phi_vector[j]), 1/phi_vector[j]) * gamma(1/phi_vector[j]))
}

## generate survival time from Weibull Cox model
get_survival_cox <- function(n, mu, phi) {
  # inverse CDF sampling
  u <- runif(n, 0, 1)
  return ((-log(u) * phi * exp(-mu))^(1 / phi))
}


## generate survival time from normal AFT model
get_survival_AFT <- function(n, mu, sigma) {
  return (exp(rnorm(n, mu, sigma)))
}

## sample S from multinomial model
# sample_S <- function(strata_coef_matrix, data_matrix_S, strata) {
#   strata_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
#   apply(strata_prob, 1, function(x)
#     sample(c(0, 1, 3), size = 1, replace = TRUE, prob = x)
#   )
# }

## sample S from multinomial model
sample_S <- function(strata_coef_matrix, data_matrix_S, strata) {
  strata_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
  apply(strata_prob, 1, function(x)
    sample(strata, size = 1, replace = TRUE, prob = x)
  )
}

## sample Y from Logistic model
sample_Y_binary <- function(outcome_coef_matrix, data_matrix) {
  outcome_matrix_linear <- data_matrix %*% outcome_coef_matrix
  outcome_matrix_prob <- plogis(outcome_matrix_linear)
  # Generate binary data with the same dimensions as outcome_matrix_prob
  # outcome_matrix <- matrix(rbinom(n = length(outcome_matrix_prob), 
  #                                 size = 1, 
  #                                 prob = outcome_matrix_prob), 
  #                          nrow = nrow(outcome_matrix_prob), 
  #                          ncol = ncol(outcome_matrix_prob))
  # 
  outcome_matrix <- matrix(NA, nrow = nrow(outcome_matrix_prob), ncol = ncol(outcome_matrix_prob))
  for (i in 1:ncol(outcome_matrix_prob)) {
    outcome_matrix[, i] <- rbinom(n = nrow(outcome_matrix_prob), size = 1, prob = outcome_matrix_prob[, i])
  }
  return (list(outcome_matrix, outcome_matrix_prob))
}

# Sample binary ICE from logistic model
sample_D <- function(ice_coef, data_matrix, Z){
  vec_linear <- data_matrix %*% ice_coef + Z
  vec_prob <- plogis(vec_linear)
  #Generate binary data with the same dimensions as outcome_matrix_prob
  D <- rbinom(n = length(vec_prob), size = 1, prob = vec_prob)
}

## sample Y from Weibull Cox model
sample_Y_cox <- function(outcome_coef_matrix, phi_vector, data_matrix_Y) {
  outcome_matrix <- data_matrix_Y %*% outcome_coef_matrix
  do.call(
    cbind, 
    map(
      1:length(phi_vector),
      ~get_survival_cox(nrow(data_matrix_Y), outcome_matrix[, .x], phi_vector[.x])
    )
  )
}

## sample Y from Weibull Cox model
sample_Y_AFT <- function(outcome_coef_matrix, sigma_vector, data_matrix_Y) {
  outcome_matrix <- data_matrix_Y %*% outcome_coef_matrix
  do.call(
    cbind, 
    map(
      1:length(sigma_vector),
      ~get_survival_AFT(nrow(data_matrix_Y), outcome_matrix[, .x], sigma_vector[.x])
    )
  )
}

# Analysis ---------------------------

stratum_surv_probability <- function(strata_coef_matrix, outcome_coef_matrix, phi_vector,
                                     data_matrix_S, data_matrix_G, time_sequence, strata_list) {
  stratum_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
  surv_prob <- calc_surv_probability(outcome_coef_matrix, phi_vector, data_matrix_G, time_sequence)
  sapply(1:length(surv_prob), 
         function(j) t(surv_prob[[j]]) %*% stratum_prob[, strata_list[j]] / sum(stratum_prob[, strata_list[j]])
  )
}

stratum_surv_probability_AFT <- function(strata_coef_matrix, outcome_coef_matrix, sigma_vector,
                                     data_matrix_S, data_matrix_G, time_sequence, strata_list) {
  stratum_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
  surv_prob <- calc_surv_probability_AFT(outcome_coef_matrix, sigma_vector, data_matrix_G, time_sequence)
  sapply(1:length(surv_prob), 
         function(j) t(surv_prob[[j]]) %*% stratum_prob[, strata_list[j]] / sum(stratum_prob[, strata_list[j]])
  )
}

stratum_RACE <- function(strata_coef_matrix, outcome_coef_matrix, phi_vector,
                                     data_matrix_S, data_matrix_G, time_sequence, strata_list) {
  stratum_prob <- calc_stratum_probability(strata_coef_matrix, data_matrix_S)
  RACE <- calc_RACE(outcome_coef_matrix, phi_vector, data_matrix_G, time_sequence)
  sapply(1:length(RACE), 
         function(j) t(RACE[[j]]) %*% stratum_prob[, strata_list[j]] / sum(stratum_prob[, strata_list[j]])
  )
}

## Reparameterization from (theta, beta_G) given in PStrata to (phi, alpha) in the paper.
## phi = exp(theta), alpha = beta_G + theta
get_surv_prob_matrices <- function(result, time_sequence, strata_list) {
  beta_S <- rstan::extract(result$post_samples, par = "beta_S")[[1]]
  beta_G <- rstan::extract(result$post_samples, par = "beta_G")[[1]]
  theta <- rstan::extract(result$post_samples, par = "theta")[[1]]
  tmp_adjust <- t(c(1, rep(0, ncol(beta_G[1, , ]) - 1)))
  data_matrix_S = result$PSobject$S.formula$fixed_eff_matrix
  data_matrix_G = result$PSobject$Y.formula$fixed_eff_matrix
  results <- pblapply(1:nrow(theta), function(i) {
    stratum_surv_probability(
      strata_coef_matrix = t(rbind(0, beta_S[i, ,])),
      outcome_coef_matrix = t(beta_G[i, ,] + theta[i, ] %*% tmp_adjust),
      phi_vector = exp(theta[i, ]),
      data_matrix_S = data_matrix_S,
      data_matrix_G = data_matrix_G,
      time_sequence = time_sequence,
      strata_list = strata_list
    )
  })
  return (results)
}

get_RACE_matrices <- function(result, time_sequence, strata_list) {
  beta_S <- rstan::extract(result$post_samples, par = "beta_S")[[1]]
  beta_G <- rstan::extract(result$post_samples, par = "beta_G")[[1]]
  theta <- rstan::extract(result$post_samples, par = "theta")[[1]]
  tmp_adjust <- t(c(1, rep(0, ncol(beta_G[1, , ]) - 1)))
  data_matrix_S = result$PSobject$S.formula$fixed_eff_matrix
  data_matrix_G = result$PSobject$Y.formula$fixed_eff_matrix
  results <- pblapply(1:nrow(theta), function(i) {
    stratum_RACE(
      strata_coef_matrix = t(rbind(0, beta_S[i, ,])),
      outcome_coef_matrix = t(beta_G[i, ,] + theta[i, ] %*% tmp_adjust),
      phi_vector = exp(theta[i, ]),
      data_matrix_S = data_matrix_S,
      data_matrix_G = data_matrix_G,
      time_sequence = time_sequence,
      strata_list = strata_list
    )
  })
  return (results)
}

## phi = exp(theta), alpha = beta_G + theta
get_surv_prob_matrices_AFT <- function(result, time_sequence, strata_list) {
  beta_S <- rstan::extract(result$post_samples, par = "beta_S")[[1]]
  beta_G <- rstan::extract(result$post_samples, par = "beta_G")[[1]]
  sigma <- rstan::extract(result$post_samples, par = "sigma")[[1]]
  data_matrix_S = result$PSobject$S.formula$fixed_eff_matrix
  data_matrix_G = result$PSobject$Y.formula$fixed_eff_matrix
  results <- pblapply(1:nrow(sigma), function(i) {
    stratum_surv_probability_AFT(
      strata_coef_matrix = t(rbind(0, beta_S[i, ,])),
      outcome_coef_matrix = t(beta_G[i, ,]),
      sigma_vector = sigma[i, ],
      data_matrix_S = data_matrix_S,
      data_matrix_G = data_matrix_G,
      time_sequence = time_sequence,
      strata_list = strata_list
    )
  })
  return (results)
}

summarize_surv_prob <- function(surv_prob_matrices) {
  num_group <- ncol(surv_prob_matrices[[1]])
  lapply(1:num_group, function(g) {
    apply(
      sapply(1:length(surv_prob_matrices), function(b) surv_prob_matrices[[b]][, g]),
      1,
      function(x) c(mean = mean(x[!is.nan(x)]), lwr = unname(quantile(x, 0.025, na.rm = T)), upr = unname(quantile(x, 0.975, na.rm = T)))
    )
  })
}

combine_fitted_true_surv_prob <- function(true_surv_prob, summary_surv_prob, 
                                          true_id, summary_id) {
  lapply(1:6, function(g) {
    true_surv_prob_g = true_surv_prob[, true_id[g]]
    summarize_surv_prob_g = summary_surv_prob[[summary_id[g]]]
    return (rbind(true = c(true_surv_prob_g), summarize_surv_prob_g))
  })
}

make_plot_data <- function(summary_data, time_sequence) {
  # syntax: summary_data[[group]][type, time]
  # output: data.frame(S, Z, time, true, mean, lwr, upr)
  S = c("never-taker", "never-taker", "complier", "complier", "always-taker", "always-taker")
  Z = c("control", "treatment", "control", "treatment", "control", "treatment")
  lapply(1:length(summary_data), function(g) {
    data.frame(S = S[g], Z = Z[g], time = time_sequence, 
               true = summary_data[[g]]["true", ], mean = summary_data[[g]]["mean", ],
               lwr = summary_data[[g]]["lwr", ], upr = summary_data[[g]]['upr', ])
  }) %>% {do.call(bind_rows, .)}
}

plot_surv_prob <- function(plot_data) {
  ggplot(plot_data) +
    geom_line(aes(x = time, y = true, linetype = "true")) +
    geom_line(aes(x = time, y = mean, linetype = "estimate")) +
    geom_ribbon(aes(x = time, ymin = lwr, ymax = upr), alpha = 0.2) +
    facet_grid(S ~ Z) + 
    xlab("Time") + ylab("Survival probability") + theme_bw()
}

add_to <- function(data, result, ER_type) {
  strata_treatments <- c(
    "never-taker (control arm)", 
    "never-taker (treatment arm)", 
    "complier (control arm)", 
    "complier (treatment arm)", 
    "always-taker (control arm)", 
    "always-taker (treatment arm)" 
  )
  for (i in 1:6) {
    new_df <- data.frame(
      true = result$combined_data[[i]]["true", ],
      mean = result$combined_data[[i]]["mean", ],
      lwr = result$combined_data[[i]]["lwr", ],
      upr = result$combined_data[[i]]["upr", ],
      time = result$time_sequence,
      ER_type = ER_type,
      strata_treatment = strata_treatments[i]
    )
    data <- bind_rows(data, new_df)
  }
  return (data)
}

# simulation <- function(
#     sim_data, truth_ER, fitted_ER, misspecification = FALSE,
#     fitted_type = 'Cox', estimand = 'prob'
# ) {
#   set.seed(0)
#   if (truth_ER) {
#     true_strata_list = c(1, 2, 2, 3)
#     true_model_id = c(1, 1, 2, 3, 4, 4)
#   }
#   else {
#     true_strata_list = c(1, 1, 2, 2, 3, 3)
#     true_model_id = c(1, 2, 3, 4, 5, 6)
#   }
#   if (fitted_ER) {
#     fitted_strata_list = c(1, 2, 2, 3)
#     fitted_model_id = c(1, 1, 2, 3, 4, 4)
#     strata = c(n = "00*", c = "01", a = "11*")
#   }
#   else {
#     fitted_strata_list = c(1, 1, 2, 2, 3, 3)
#     fitted_model_id = c(1, 2, 3, 4, 5, 6)
#     strata = c(n = "00", c = "01", a = "11")
#   }
#   
#   if (misspecification) {
#     S.formula = Z + D ~ X1
#     Y.formula = Y + delta ~ X1
#   }
#   else {
#     S.formula = Z + D ~ X1 + X2
#   Y.formula = Y + delta ~ X1 + X2
#   }
#   
#   time_sequence <- seq(0.01, 15, length.out = 100)
#   
#   result <- PStrata::PStrata(
#     S.formula = S.formula,
#     Y.formula = Y.formula,
#     Y.family = survival(fitted_type),
#     data = sim_data$data,
#     strata = strata,
#     warmup = 500, iter = 1000,
#     chains = 6, core = 6,
#     prior_intercept = prior_normal(0, 10),
#     survival.time.points = seq(0.01, 15, length.out = 100)
#   )
#   
#   data_matrix <- cbind(1, sim_data$data$X1, sim_data$data$X2, 0)
#   
#   if (estimand == 'prob') {
#     true_surv_prob <- stratum_surv_probability(
#       sim_data$strata_coef_matrix, sim_data$outcome_coef_matrix, sim_data$phi,
#       data_matrix, data_matrix, time_sequence, true_strata_list
#     )
#     if (fitted_type == "Cox") {
#       fitted_surv_prob <- get_surv_prob_matrices(result, time_sequence, fitted_strata_list)
#     }
#     else {
#       fitted_surv_prob <- get_surv_prob_matrices_AFT(result, time_sequence, fitted_strata_list)
#     }
#     summary_surv_prob <- summarize_surv_prob(fitted_surv_prob)
#     combined_data <- combine_fitted_true_surv_prob(
#       true_surv_prob, summary_surv_prob, true_model_id, fitted_model_id
#     )
#   }
#   else {
#     true_RACE <- stratum_RACE(
#       sim_data$strata_coef_matrix, sim_data$outcome_coef_matrix, sim_data$phi,
#       data_matrix, data_matrix, time_sequence, fitted_strata_list
#     )
#     if (fitted_type == "Cox") {
#       fitted_RACE <- get_RACE_matrices(result, time_sequence, fitted_strata_list)
#     }
#     else {
#       fitted_RACE <- get_RACE_matrices_AFT(result, time_sequence, fitted_strata_list)
#     }
#     summary_RACE <- summarize_surv_prob(fitted_RACE)
#     combined_data <- combine_fitted_true_surv_prob(
#       true_RACE, summary_RACE, true_model_id, fitted_model_id
#     )
#   }
#   
#   gPlot <- plot_surv_prob(make_plot_data(combined_data, time_sequence))
#   return (list(
#     time_sequence = time_sequence, data = sim_data, PStrata_result = result, combined_data = combined_data, plot = gPlot
#   ))
# }


#================================PI function==========================
# The initial value of EM algorithm is set to 0
# Covariates are recommended to be scaled 
# "times" are a vector of the times we are interested in 
# epsz are the sensitivity parameters (z=0,1), we let epsz(X)=exp(-X %*% etaz), where etaz are pre-determined
####################### MR Estimators ############################
est.MR.01.1 = function(times,U,delta,z,s,Xtp,Xps,XC,XT,eps1){
  u.order <- order(U)
  U <- U[u.order]
  Xtp <- Xtp[u.order, ]
  Xps <- Xps[u.order, ]
  XC <- XC[u.order, ]
  XT <- XT[u.order, ]
  delta <- delta[u.order]
  z <- z[u.order]
  s <- s[u.order]
  n = length(z)
  fit.tp <- glm(z~Xtp, family=binomial(link = "logit"))
  g0_est <- fit.tp$coef[1]
  g1_est <- fit.tp$coef[-1]
  tp_est <- as.vector(exp(g0_est+ Xtp %*% g1_est)/(1 + exp(g0_est+ Xtp %*% g1_est)))  #estimation of tp
  Xpsint = cbind(1,Xps)
  set00 = subset(1:n,(s==0)&(z==0))
  set01 = subset(1:n,(s==1)&(z==0))
  set10 = subset(1:n,(s==0)&(z==1))
  set11 = subset(1:n,(s==1)&(z==1))
  n00 = length(set00)
  n01 = length(set01)
  n10 = length(set10)
  n11 = length(set11)
  p = length(Xps[1,])
  t1 = t2 = numeric(p+1)
  EM = function(t1,t2){
    theta0 = theta3 = theta4 = numeric(p+1)
    theta1 = t1
    theta2 = t2
    pr1 = pr2 = numeric(length(set11))
    pr3 = pr4 = numeric(length(set00))
    for(j in 1:3){                         ###EM?㷨??theta
      pr1 = 1/(1+exp(Xpsint[set11,]%*%theta1))
      pr2 = 1-pr1
      pr3 = 1/(1+exp(Xpsint[set00,]%*%theta2))
      pr4 = 1-pr3
      V = c(rep(0,length(set11)),rep(1,length(set11)),rep(2,length(set10)),rep(1,length(set01)),
            rep(0,length(set00)),rep(2,length(set00)))
      tmpX<-rbind(Xpsint[set11,],Xpsint[set11,],Xpsint[set10,],Xpsint[set01,],Xpsint[set00,],Xpsint[set00,])
      weight = c(pr1,pr2,rep(1,length(set10) + length(set01)),pr3,pr4)
      fit<-multinom(V~tmpX[,2:(p+1)],data = as.data.frame(tmpX),weights = weight,trace=F)
      theta1<-coefficients(fit)[1,]
      theta2<-coefficients(fit)[2,]
    }
    while (sum(abs(theta1-theta3))/sum(abs(theta1))>10^-5 || sum(abs(theta2-theta4))/sum(abs(theta2))>10^-5){
      theta3 = theta1
      theta4 = theta2
      pr1 = 1/(1+exp(Xpsint[set11,]%*%theta1))
      pr2 = 1-pr1
      pr3 = 1/(1+exp(Xpsint[set00,]%*%theta2))
      pr4 = 1-pr3
      V = c(rep(0,length(set11)),rep(1,length(set11)),rep(2,length(set10)),rep(1,length(set01)),
            rep(0,length(set00)),rep(2,length(set00)))
      tmpX<-rbind(Xpsint[set11,],Xpsint[set11,],Xpsint[set10,],Xpsint[set01,],Xpsint[set00,],Xpsint[set00,])
      weight = c(pr1,pr2,rep(1,length(set10) + length(set01)),pr3,pr4)
      fit<-multinom(V~tmpX[,2:(p+1)],data = as.data.frame(tmpX),weights = weight,trace=F)
      theta1<-coefficients(fit)[1,]
      theta2<-coefficients(fit)[2,]
    }
    return(c(theta1,theta2))
  }
  tmp = EM(t1,t2)
  ps01_est = 1/(1+exp(Xpsint%*%tmp[1:(p+1)])+exp(Xpsint%*%tmp[(p+2):(2*p+2)]))
  ps11_est = exp(Xpsint%*%tmp[1:(p+1)])/(1+exp(Xpsint%*%tmp[1:(p+1)])+exp(Xpsint%*%tmp[(p+2):(2*p+2)]))
  ps00_est = exp(Xpsint%*%tmp[(p+2):(2*p+2)])/(1+exp(Xpsint%*%tmp[1:(p+1)])+exp(Xpsint%*%tmp[(p+2):(2*p+2)]))
  p1_est = 1-ps00_est
  p0_est = ps11_est
  p1.mean = mean(z*(s-p1_est)/tp_est+p1_est)
  p0.mean = mean((1-z)*(s-p0_est)/(1-tp_est)+p0_est)
  w01_est = (p1_est-p0_est)/(p1.mean-p0.mean)
  ######## Sensitivity Parameters ##################
  w1.01.est = (eps1 * p1_est)/(eps1 * (p1_est-p0_est)+p0_est)
  ################## PH Model for T and C #########################
  n11 = length(set11)
  U11 = U[set11]
  delta11 = delta[set11]
  cdelta11 = 1-delta11
  XT11 = XT[set11,]
  XC11 = XC[set11,]
  Tdata = data.frame(U11,delta11,XT11)
  Cdata = data.frame(U11,cdelta11,XC11)
  fit.C11<-coxph(Surv(U11, cdelta11) ~ XC11,data=Cdata)
  fit.T11<-coxph(Surv(U11,delta11) ~ XT11,data=Tdata)
  coeC = fit.C11$coefficients
  coeT = fit.T11$coefficients
  
  CumHazT = CumHazC = matrix(0,n,n)            #cumulative hazard of T and C
  basehazC = basehaz(fit.C11,centered = F)
  basehazT = basehaz(fit.T11,centered = F)
  basehazT1 = basehazC1 = numeric(n)           #basehaz for T and C in n time points
  for(i in 1:n){
    basehazC1[i] = basehazC[max1(which(basehazC$time<=U[i])),1]
    basehazT1[i] = basehazT[max1(which(basehazT$time<=U[i])),1]
  }
  
  CumHazC = outer(c(exp(XC %*% coeC)),basehazC1)
  CumHazT = outer(c(exp(XT %*% coeT)),basehazT1)
  SPC = exp(-CumHazC)         #survival prob of C
  SPT = exp(-CumHazT)
  #############################martingale################
  # \int{Y(r)*lambda_c(r,X)/[K_c(r,X)*H(r,X)]}
  q = matrix(NA,n,n)
  q[,1] = -CumHazC[,1]/(SPC[,1]*SPT[,1])
  for(j in 2:n)
    q[j:n,j] = q[j:n,j-1] - (CumHazC[j:n,j] - CumHazC[j:n,j-1])/(SPC[j:n,j]*SPT[j:n,j])
  for(i in 1:(n-1)){
    q[i,(i+1):n] = rep(q[i,i],n-i)
    if(delta[i]==0)
      q[i,(i+1):n] = q[i,(i+1):n] + 1/(SPC[i,i]*SPT[i,i])
  }
  ################## Estimation ##############
  S<-function(u){
    pos = min(sum(u>U)+1,n)
    pos2 = max1(which(U<=u))
    est =  (sum( w1.01.est[pos:n]*w01_est[pos:n] * (s[pos:n]) * (z[pos:n])/((p1_est[pos:n]) * (tp_est[pos:n]) * SPC[pos:n,pos2]))
            + sum(w1.01.est*w01_est * (1-(z)/(tp_est)) * SPT[,pos2])
            - sum(((1-z)*(s-p0_est)/(1-tp_est) + (p0_est)*(z)*(p1_est-s)/((p1_est)*(tp_est))) * SPT[,pos2] * w1.01.est^2/eps1)/(p1.mean-p0.mean)
            + sum( w1.01.est*w01_est * (s) * (z) * SPT[,pos2] * q[,pos]/((p1_est)*(tp_est)))
    )
    return(est/n)
  }
  Sv<-function(u){             
    l = length(u)
    result = numeric(l)
    for(i in 1:l)
      result[i] = S(u[i])
    return(result)
  }
  return(Sv(times))
} 
est.MR.01.0 = function(times,U,delta,z,s,Xtp,Xps,XC,XT,eps0){
  u.order <- order(U)
  U <- U[u.order]
  Xtp <- Xtp[u.order, ]
  Xps <- Xps[u.order, ]
  XT <- XT[u.order, ]
  XC <- XC[u.order, ]
  delta <- delta[u.order]
  z <- z[u.order]
  s <- s[u.order]
  n = length(z)
  fit.tp <- glm(z~Xtp, family=binomial(link = "logit"))
  g0_est <- fit.tp$coef[1]
  g1_est <- fit.tp$coef[-1]
  tp_est <- as.vector(exp(g0_est+ Xtp %*% g1_est)/(1 + exp(g0_est+ Xtp %*% g1_est)))  #estimation of tp
  Xpsint = cbind(1,Xps)
  set00 = subset(1:n,(s==0)&(z==0))
  set01 = subset(1:n,(s==1)&(z==0))
  set10 = subset(1:n,(s==0)&(z==1))
  set11 = subset(1:n,(s==1)&(z==1))
  n00 = length(set00)
  n01 = length(set01)
  n10 = length(set10)
  n11 = length(set11)
  p = length(Xps[1,])
  t1 = t2 = numeric(p+1)
  EM = function(t1,t2){
    theta0 = theta3 = theta4 = numeric(p+1)
    theta1 = t1
    theta2 = t2
    pr1 = pr2 = numeric(length(set11))
    pr3 = pr4 = numeric(length(set00))
    for(j in 1:3){                         ###EM?㷨??theta
      pr1 = 1/(1+exp(Xpsint[set11,]%*%theta1))
      pr2 = 1-pr1
      pr3 = 1/(1+exp(Xpsint[set00,]%*%theta2))
      pr4 = 1-pr3
      V = c(rep(0,length(set11)),rep(1,length(set11)),rep(2,length(set10)),rep(1,length(set01)),
            rep(0,length(set00)),rep(2,length(set00)))
      tmpX<-rbind(Xpsint[set11,],Xpsint[set11,],Xpsint[set10,],Xpsint[set01,],Xpsint[set00,],Xpsint[set00,])
      weight = c(pr1,pr2,rep(1,length(set10) + length(set01)),pr3,pr4)
      fit<-multinom(V~tmpX[,2:(p+1)],data = as.data.frame(tmpX),weights = weight,trace=F)
      theta1<-coefficients(fit)[1,]
      theta2<-coefficients(fit)[2,]
    }
    while (sum(abs(theta1-theta3))/sum(abs(theta1))>10^-4 || sum(abs(theta2-theta4))/sum(abs(theta2))>10^-4){
      theta3 = theta1
      theta4 = theta2
      pr1 = 1/(1+exp(Xpsint[set11,]%*%theta1))
      pr2 = 1-pr1
      pr3 = 1/(1+exp(Xpsint[set00,]%*%theta2))
      pr4 = 1-pr3
      V = c(rep(0,length(set11)),rep(1,length(set11)),rep(2,length(set10)),rep(1,length(set01)),
            rep(0,length(set00)),rep(2,length(set00)))
      tmpX<-rbind(Xpsint[set11,],Xpsint[set11,],Xpsint[set10,],Xpsint[set01,],Xpsint[set00,],Xpsint[set00,])
      weight = c(pr1,pr2,rep(1,length(set10) + length(set01)),pr3,pr4)
      fit<-multinom(V~tmpX[,2:(p+1)],data = as.data.frame(tmpX),weights = weight,trace=F)
      theta1<-coefficients(fit)[1,]
      theta2<-coefficients(fit)[2,]
    }
    return(c(theta1,theta2))
  }
  tmp = EM(t1,t2)
  ps01_est = 1/(1+exp(Xpsint%*%tmp[1:(p+1)])+exp(Xpsint%*%tmp[(p+2):(2*p+2)]))
  ps11_est = exp(Xpsint%*%tmp[1:(p+1)])/(1+exp(Xpsint%*%tmp[1:(p+1)])+exp(Xpsint%*%tmp[(p+2):(2*p+2)]))
  ps00_est = exp(Xpsint%*%tmp[(p+2):(2*p+2)])/(1+exp(Xpsint%*%tmp[1:(p+1)])+exp(Xpsint%*%tmp[(p+2):(2*p+2)]))
  p1_est = 1-ps00_est
  p0_est = ps11_est
  p1.mean = mean(z*(s-p1_est)/tp_est+p1_est)
  p0.mean = mean((1-z)*(s-p0_est)/(1-tp_est)+p0_est)
  w01_est = (p1_est-p0_est)/mean(p1.mean-p0.mean)
  w0.01.est = (eps0 * (1-p0_est))/(eps0 * (p1_est-p0_est) + (1-p1_est))
  ################## PH Model for T and C #########################
  n00 = length(set00)
  U00 = U[set00]
  delta00 = delta[set00]
  cdelta00 = 1-delta00
  XT00 = XT[set00,]
  XC00 = XC[set00,]
  Tdata = data.frame(U00,delta00,XT00)
  Cdata = data.frame(U00,cdelta00,XC00)
  fit.C00<-coxph(Surv(U00, cdelta00) ~ XC00,data=Cdata)
  fit.T00<-coxph(Surv(U00,delta00) ~ XT00,data=Tdata)
  coeC = fit.C00$coefficients
  coeT = fit.T00$coefficients
  
  CumHazT = CumHazC = matrix(0,n,n)            #cumulative hazard of T and C
  basehazC = basehaz(fit.C00,centered = F)
  basehazT = basehaz(fit.T00,centered = F)
  basehazT1 = basehazC1 = numeric(n)           #basehaz for T and C in n time points
  for(i in 1:n){
    basehazC1[i] = basehazC[max1(which(basehazC$time<=U[i])),1]
    basehazT1[i] = basehazT[max1(which(basehazT$time<=U[i])),1]
  }
  
  CumHazC = outer(c(exp(XC %*% coeC)),basehazC1)
  CumHazT = outer(c(exp(XT %*% coeT)),basehazT1)
  SPC = exp(-CumHazC)         #survival prob of C
  SPT = exp(-CumHazT)
  #############################martingale################
  # \int{Y(r)*lambda_c(r,X)/[K_c(r,X)*H(r,X)]}
  q = matrix(NA,n,n)
  q[,1] = -CumHazC[,1]/(SPC[,1]*SPT[,1])
  for(j in 2:n)
    q[j:n,j] = q[j:n,j-1] - (CumHazC[j:n,j] - CumHazC[j:n,j-1])/(SPC[j:n,j]*SPT[j:n,j])
  for(i in 1:(n-1)){
    q[i,(i+1):n] = rep(q[i,i],n-i)
    if(delta[i]==0)
      q[i,(i+1):n] = q[i,(i+1):n] + 1/(SPC[i,i]*SPT[i,i])
  }
  ################## Estimation ##############
  S<-function(u){
    pos = min(sum(u>U)+1,n)
    pos2 = max1(which(U<=u))
    est =  (sum( w0.01.est[pos:n]*w01_est[pos:n] * (1-s[pos:n]) * (1-z[pos:n])/((1-p0_est[pos:n]) * (1-tp_est[pos:n]) * SPC[pos:n,pos2]))
            +  sum(w0.01.est*w01_est * (1-(1-z)/(1-tp_est)) * SPT[,pos2])
            +  sum((z*(s-p1_est)/tp_est + (1-p1_est)*(1-z)*(p0_est-s)/((1-p0_est)*(1-tp_est))) * SPT[,pos2] * w0.01.est^2)/((p1.mean-p0.mean)*eps0)
            +  sum(w0.01.est*w01_est * (1-s) * (1-z) * SPT[,pos2] * q[,pos]/((1-p0_est)*(1-tp_est)))
    )
    return(est/n)
  }
  Sv<-function(u){             
    l = length(u)
    result = numeric(l)
    for(i in 1:l)
      result[i] = S(u[i])
    return(result)
  }
  # p <- ggplot() + xlim(a,b) + ylim(0,1) + geom_function(fun = Sv,color = "red")
  # print(p) 
  return(Sv(times))
} 

est.MR.SA.mono <- function(times, propensity.model, principal.model, censor.model, failure.model, data, 
                           zeta,
                           estimand = c("NACE", "CACE", "AACE", "DACE")) {
  # propensity: Z ~ Xtp
  # principal: S ~ Xps
  # censor: delta ~ XC
  # failure: U ~ XT
  tm <- Sys.time()
  lower_bound_vec <- function(vals, seq) {
    # Find the position of the first element that is not less than 
    # a given value, or length(inc_seq) + 1 if no such element is found.
    # Time complexity: O(m + n) - often faster than O(n log m) by 
    # doing elementwise binary search.
    result <- integer(length(vals))
    p <- 1
    q <- 1
    while (p <= length(vals)) {
      while (q <= length(seq) && seq[q] < vals[p])
        q <- q + 1
      result[p] <- q
      p <- p + 1
    }
    return (result)
  }
  
  upper_bound_vec <- function(vals, seq) {
    # Find the position of the first element that is greater than 
    # a given value, or length(inc_seq) + 1 if no such element is found.
    # Time complexity: O(m + n) - often faster than O(n log m) by 
    # applying the univariate version of upper_bound
    result <- integer(length(vals))
    p <- 1
    q <- 1
    while (p <= length(vals)) {
      while (q <= length(seq) && seq[q] <= vals[p])
        q <- q + 1
      result[p] <- q
      p <- p + 1
    }
    return (result)
  }
  
  mono <- function(seq){
    cur <- Inf
    for (i in 1:length(seq)) {
      seq[i] <- min(cur, seq[i])
      cur <- min(cur, seq[i])
    }
    return (seq)
  }
  
  Z_name <- propensity.model[[2]]
  S_name <- principal.model[[2]]
  Delta_name <- censor.model[[2]]
  U_name <- failure.model[[2]]
  
  data <- data %>% arrange({{U_name}})
  n <- nrow(data)
  z <- data %>% pull({{Z_name}})
  s <- data %>% pull({{S_name}})
  U <- data %>% pull({{U_name}})
  delta <- data %>% pull({{Delta_name}})
  XC <- model.matrix(censor.model, data)[, -1]
  XT <- model.matrix(failure.model, data)[, -1]
  
  fit.tp <- glm(propensity.model, family=binomial(link = "logit"), data = data)
  tp_est <- fit.tp$fitted.values  #estimation of tp
  fit.p1 <- glm(principal.model, family = binomial(link = "logit"), data = data[z == 1, ])
  p1_est <- predict(fit.p1, data, type = "response")
  fit.p0 <- glm(principal.model, family = binomial(link = "logit"), data = data[z == 0, ])
  p0_est <- predict(fit.p0, data, type = "response")
  
  ps00_est = 1-p0_est-(p1_est-p0_est)/(1-zeta)
  ps11_est = p1_est-(p1_est-p0_est)/(1-zeta)
  ps01_est = (p1_est-p0_est)/(1-zeta)
  ps10_est = zeta*(p1_est-p0_est)/(1-zeta)
  p1.mean.dr = mean((z*(s-p1_est)/tp_est+p1_est))
  p0.mean.dr = mean(((1-z)*(s-p0_est)/(1-tp_est)+p0_est))
  # print(min(1-(p1.mean.dr-p0.mean.dr)/pmin(p1.mean.dr,1-p0.mean.dr)))
  # if(zeta>min(1-(p1.mean.dr-p0.mean.dr)/pmin(p1.mean.dr,1-p0.mean.dr))){
  #   print("zeta out of bound")
  # }
  
  p1.mean.zeta = mean((z*(s-p1_est)/tp_est+p1_est)/(1-zeta))
  p0.mean.zeta = mean(((1-z)*(s-p0_est)/(1-tp_est)+p0_est)/(1-zeta))
  
  p1.mean = mean(z*(s-p1_est)/tp_est+p1_est)
  p0.mean = mean((1-z)*(s-p0_est)/(1-tp_est)+p0_est)
  e01.zeta = p1.mean.zeta-p0.mean.zeta
  e00.zeta = 1-p0.mean-p1.mean.zeta+p0.mean.zeta
  e11.zeta = p1.mean-p1.mean.zeta+p0.mean.zeta
  e10.zeta = p1.mean.zeta-p0.mean.zeta-p1.mean+p0.mean
  w01_est = ps01_est/(e01.zeta)
  w00_est = ps00_est/(e00.zeta)
  w11_est = ps11_est/(e11.zeta)
  w10_est = ps10_est/(e10.zeta)
  
  xi1 = xi0 = 0
  eta1 = eta0 = 1
  eps1 <- exp(xi1 * (times/max(times))^eta1)
  eps0 <- exp(xi0 * (times/max(times))^eta0)
  w1c <- sapply(eps1, function(x) (x * ps01_est + x * ps11_est) / (x * ps01_est + ps11_est))
  w0c <- sapply(eps0, function(x) (x * ps01_est + x * ps00_est) / (x * ps01_est + ps00_est))
  w0n <- sapply(eps0, function(x) (ps01_est + ps00_est) / (x * ps01_est + ps00_est))
  w1a <- sapply(eps1, function(x) (ps01_est + ps11_est) / (x * ps01_est + ps11_est))
  
  pos1 <- pmin(lower_bound_vec(times, U), n)
  pos2 <- pmax(upper_bound_vec(times, U) - 1, 1)
  time_id <- 1:length(times)
  
  subsets <- list(
    "00" = s == 0 & z == 0,
    "01" = s == 1 & z == 0,
    "10" = s == 0 & z == 1,
    "11" = s == 1 & z == 1
  )
  
  call_censor <- call(
    "~",
    as.call(list(
      quote(Surv), U_name, as.call(list(quote(`-`), 1, Delta_name))
    )),
    censor.model[[3]]
  )
  
  call_failure <- call(
    "~",
    as.call(list(
      quote(Surv), U_name, Delta_name
    )),
    failure.model[[3]]
  )
  
  new_outer <- function(a, b) {
    return (structure(list(
      a = a,
      b = b
    ), class = 'outer'))
  }
  
  `[.outer` <- function(x, i, j) {
    if (missing(i)) {
      return (x$a * x$b[j])
    }
    if (missing(j)) {
      return (x$a[i] * x$b)
    }
    return (x$a[i] * x$b[j])
  }
  
  calc_cumhaz <- function(ss) {
    
    formula_censor <- eval(call_censor)
    formula_failure <- eval(call_failure)
    
    fit <- list(
      T = coxph(formula_failure, data = data[ss, ]),
      C = coxph(formula_censor, data = data[ss, ])
    )
    tmp_func <- function(fit, X, delta) {
      #basehaz <- basehaz(fit, centered = F)
      #basehaz1 <- basehaz$hazard[pmax(upper_bound_vec(U, basehaz$time) - 1L, 1L)]
      #cumhaz <- outer(c(exp(X %*% fit$coefficients)), basehaz1)
      b_est <- fit$coefficients
      s1 <- cumsum(ss[n:1]*exp(as.matrix(X[n:1,]) %*% b_est))
      dls <- ss * delta / s1[n:1]
      dls[(1:length(dls))[is.nan(dls)]] <- 0
      cumhaz <- new_outer(c(exp(X %*% fit$coefficients)), cumsum(dls))
      #cumhaz <- outer(c(exp(X %*% fit$coefficients)), cumsum(dls))
      return (cumhaz)
    }
    return (list(
      T = tmp_func(fit$T, XT, delta),
      C = tmp_func(fit$C, XC, 1 - delta)
    ))
  }
  
  calc_q <- function(cumhaz) {
    q <- diag(NA, n)
    q[,1] = -cumhaz$C[,1] * exp(cumhaz$C[,1] + cumhaz$T[,1])
    
    for(j in 2:n) {
      q[j:n,j] = q[j:n,j-1] - (cumhaz$C[j:n,j] - cumhaz$C[j:n,j-1]) * exp(cumhaz$C[j:n,j] + cumhaz$T[j:n,j])
    }
    for(i in 1:(n-1)){
      q[i,(i+1):n] = q[i,i] + (1 - delta[i]) * exp(cumhaz$C[i, i] + cumhaz$T[i, i])
    }
    return (q)
  }
  
  subset_ind <- c(
    "NACE" %in% estimand || "CACE" %in% estimand, # 00
    "AACE" %in% estimand || "DACE" %in% estimand, # 01
    "NACE" %in% estimand || "DACE" %in% estimand, # 10
    "CACE" %in% estimand || "AACE" %in% estimand # 11
  )
  
  subsets_to_use <- subsets[subset_ind]
  
  # calculate q only for needed
  
  cumhaz <- lapply(subsets_to_use, calc_cumhaz)
  q <- lapply(cumhaz, calc_q)
  gc()
  
  result <- list()
  
  gen_result <- function(func1, func0) {
    est1 <- mapply(function(x, y, t) func1(x, y, t) / n, pos1, pos2, time_id)
    est0 <- mapply(function(x, y, t) func0(x, y, t) / n, pos1, pos2, time_id)
    z1_mono = mono(est1)
    z0_mono = mono(est0)
    
    return (list(
      z1 = est1,
      z0 = est0,
      diff = est1 - est0,
      z1_mono = z1_mono,
      z0_mono = z0_mono,
      diff_mono = z1_mono - z0_mono
    ))
  }
  
  NACE_1 <- function(x, y, t) {
    s1 <- w00_est * (1 - s) * z
    s2 <- (1 - p1_est) * tp_est
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["10"]]$C[x:n, y]))
    ) + sum(
      w00_est * (1 - z / tp_est) * exp(-cumhaz[["10"]]$T[, y])
    ) + sum(
      zeta * ((1 - z) * (s - p0_est) / (1 - tp_est) + (1 - p0_est) * z * (p1_est - s) / s2) * exp(-cumhaz[["10"]]$T[, y])
    ) / (e00.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["10"]]$T[, y]) * q[["10"]][, x] / s2
    )
  }
  
  NACE_0 <- function(x, y, t) {
    s1 <- w00_est * (1 - s) * (1 - z)
    s2 <- (1 - p0_est) * (1 - tp_est)
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["00"]]$C[x:n, y]))
    ) + sum(
      w00_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["00"]]$T[, y])
    ) - sum(
      (z * (s - p1_est) / tp_est + (1 - p1_est) * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["00"]]$T[, y])
    ) / (e00.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["00"]]$T[, y]) * q[["00"]][, x] / s2
    )
  }
  
  CACE_1 <- function(x, y, t) {
    s1 <- w01_est * s * z
    s2 <- p1_est * tp_est
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["11"]]$C[x:n, y]))
    ) + sum(
      w01_est * (1 - z / tp_est) * exp(-cumhaz[["11"]]$T[, y])
    ) - sum(
      ((1 - z) * (s - p0_est) / (1 - tp_est) + p0_est * z * (p1_est - s) / s2) * exp(-cumhaz[["11"]]$T[, y])
    ) / (e01.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["11"]]$T[, y]) * q[["11"]][, x] / s2
    )
  }
  
  CACE_0 <-  function(x, y, t) {
    s1 <- w01_est * (1 - s) * (1 - z)
    s2 <- (1 - p0_est) * (1 - tp_est)
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["00"]]$C[x:n, y]))
    ) + sum(
      w01_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["00"]]$T[, y])
    ) + sum(
      (z * (s - p1_est) / tp_est + (1 - p1_est) * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["00"]]$T[, y])
    ) / (e01.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["00"]]$T[, y]) * q[["00"]][, x] / s2
    )
  }
  
  AACE_1 <-  function(x, y, t) {
    s1 <- w11_est * s * z
    s2 <- p1_est * tp_est
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["11"]]$C[x:n, y]))
    ) + sum(
      w11_est * (1 - z / tp_est) * exp(-cumhaz[["11"]]$T[, y])
    ) + sum(
      ((1 - z) * (s - p0_est) / (1 - tp_est) + p0_est * z * (p1_est - s) / s2) * exp(-cumhaz[["11"]]$T[, y])
    ) / (e11.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["11"]]$T[, y]) * q[["11"]][, x] / s2
    )
  }
  
  AACE_0 <- function(x, y, t) {
    s1 <- w11_est * s * (1 - z)
    s2 <- p0_est * (1 - tp_est)
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["01"]]$C[x:n, y]))
    ) + sum(
      w11_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["01"]]$T[, y])
    ) + sum(
      zeta * (z * (p1_est - s) / tp_est + p1_est * (1 - z) * (s - p0_est) / s2) * exp(-cumhaz[["01"]]$T[, y])
    ) / (e11.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["01"]]$T[, y]) * q[["01"]][, x] / s2
    )
  }
  
  DACE_1 <-  function(x, y, t) {
    s1 <- w10_est * (1 - s) * z
    s2 <- (1 - p1_est) * tp_est
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["10"]]$C[x:n, y]))
    ) + sum(
      w10_est * (1 - z / tp_est) * exp(-cumhaz[["10"]]$T[, y])
    ) + sum(
      zeta * ((1 - z) * (p0_est - s) / (1 - tp_est) + (1 - p0_est) * z * (s - p1_est) / s2) * exp(-cumhaz[["10"]]$T[, y])
    ) / (e10.zeta * (1-zeta)) + sum(
      s1 * exp(-cumhaz[["10"]]$T[, y]) * q[["10"]][, x] / s2
    )
  }
  
  DACE_0 <- function(x, y, t) {
    s1 <- w10_est * s * (1 - z)
    s2 <- p0_est * (1 - tp_est)
    sum(
      s1[x:n] / (s2[x:n] * exp(-cumhaz[["01"]]$C[x:n, y]))
    ) + sum(
      w10_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["01"]]$T[, y])
    ) + sum(
      zeta * (z * (s - p1_est) / tp_est + p1_est * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["01"]]$T[, y])
    ) / (e10.zeta * (1 - zeta)) + sum(
      s1 * exp(-cumhaz[["01"]]$T[, y]) * q[["01"]][, x] / s2
    )
  }
  
  if ("NACE" %in% estimand) 
    result$NACE <- gen_result(NACE_1, NACE_0)
  if ("CACE" %in% estimand)
    result$CACE <- gen_result(CACE_1, CACE_0)
  if ("AACE" %in% estimand)
    result$AACE <- gen_result(AACE_1, AACE_0)
  if ("DACE" %in% estimand)
    result$DACE <- gen_result(DACE_1, DACE_0)
  result$bound <- min(1-(p1.mean.dr-p0.mean.dr)/pmin(p1.mean.dr,1-p0.mean.dr))
  class(result) <- "est.MR.SA"
  return (result)
}

ACE.SA.mono <- function(times, propensity.model, principal.model, censor.model, failure.model, data, 
                        zeta,
                        estimand = c("NACE", "CACE", "AACE", "DACE"), bootstrap = 50) {
  gc()
  mean_result <- est.MR.SA(times, propensity.model, principal.model, censor.model, failure.model, data, 
                           zeta, estimand)
  if (bootstrap == 0)
    return (mean_result)
  bs_result <- c()
  for (name in estimand) {
    z1.B <- z0.B <- diff.B <- z1_mono.B <- z0_mono.B <- diff_mono.B <- matrix(nrow = length(times), ncol = bootstrap)
    bs_result[[name]] <- list(
      z1.B = z1.B,
      z0.B = z0.B,
      diff.B = diff.B,
      z1_mono.B = z1_mono.B,
      z0_mono.B = z0_mono.B,
      diff_mono.B = diff_mono.B
    )
  }
  gc()
  pb <- progress::progress_bar$new(total = bootstrap, format = "(:spin) [:bar] :percent. :elapsed | :eta")
  cnt <- 0
  for (i in 1:bootstrap) {
    gc()
    data.B <- slice_sample(data, n = nrow(data), replace = T)
    result <- est.MR.SA(times, propensity.model, principal.model, censor.model, failure.model, data.B, 
                        zeta, estimand)
    if (result$bound < zeta) {
      cnt <- cnt + 1
      cat("Resample: ", cnt, " with bound ", result$bound, "\n")
      i <- i - 1
      next
    }
    pb$tick()
    for (name in estimand) {
      bs_result[[name]]$z1.B[, i] <- result[[name]]$z1
      bs_result[[name]]$z0.B[, i] <- result[[name]]$z0
      bs_result[[name]]$diff.B[, i] <- result[[name]]$diff
      bs_result[[name]]$z1_mono.B[, i] <- result[[name]]$z1_mono
      bs_result[[name]]$z0_mono.B[, i] <- result[[name]]$z0_mono
      bs_result[[name]]$diff_mono.B[, i] <- result[[name]]$diff_mono
    }
    gc()
  }
  r_val <- list(
    estimate = mean_result,
    bootstrap = bs_result
  )
  class(r_val) <- "est.MR.SA.bootstrap"
  return (r_val)
}

##########################################################

est.MR.SA <- function(times, propensity.model, principal.model, censor.model, failure.model, data, 
                      xi0, xi1, eta0, eta1,
                      estimand = c("NACE", "CACE", "AACE")) {
  # propensity: Z ~ Xtp
  # principal: S ~ Xps
  # censor: delta ~ XC
  # failure: U ~ XT
  tm <- Sys.time()
  lower_bound_vec <- function(vals, seq) {
    # Find the position of the first element that is not less than 
    # a given value, or length(inc_seq) + 1 if no such element is found.
    # Time complexity: O(m + n) - often faster than O(n log m) by 
    # doing elementwise binary search.
    result <- integer(length(vals))
    p <- 1
    q <- 1
    while (p <= length(vals)) {
      while (q <= length(seq) && seq[q] < vals[p])
        q <- q + 1
      result[p] <- q
      p <- p + 1
    }
    return (result)
  }
  
  upper_bound_vec <- function(vals, seq) {
    # Find the position of the first element that is greater than 
    # a given value, or length(inc_seq) + 1 if no such element is found.
    # Time complexity: O(m + n) - often faster than O(n log m) by 
    # applying the univariate version of upper_bound
    result <- integer(length(vals))
    p <- 1
    q <- 1
    while (p <= length(vals)) {
      while (q <= length(seq) && seq[q] <= vals[p])
        q <- q + 1
      result[p] <- q
      p <- p + 1
    }
    return (result)
  }
  
  mono <- function(seq){
    cur <- Inf
    for (i in 1:length(seq)) {
      seq[i] <- min(cur, seq[i])
      cur <- min(cur, seq[i])
    }
    return (seq)
  }
  
  Z_name <- propensity.model[[2]]
  S_name <- principal.model[[2]]
  Delta_name <- censor.model[[2]]
  U_name <- failure.model[[2]]
  
  data <- data %>% arrange({{U_name}})
  n <- nrow(data)
  z <- data %>% pull({{Z_name}})
  s <- data %>% pull({{S_name}})
  U <- data %>% pull({{U_name}})
  delta <- data %>% pull({{Delta_name}})
  XC <- model.matrix(censor.model, data)[, -1]
  XT <- model.matrix(failure.model, data)[, -1]
  
  fit.tp <- glm(propensity.model, family=binomial(link = "logit"), data = data)
  tp_est <- fit.tp$fitted.values  #estimation of tp
  fit.p1 <- glm(principal.model, family = binomial(link = "logit"), data = data[z == 1, ])
  p1_est <- predict(fit.p1, data, type = "response")
  fit.p0 <- glm(principal.model, family = binomial(link = "logit"), data = data[z == 0, ])
  p0_est <- predict(fit.p0, data, type = "response")
  
  ps01_est <<- pmax(p1_est - p0_est,0)
  ps00_est <<- 1-p1_est
  ps11_est <<- p0_est
  p1.mean <- mean(z*(s-p1_est)/tp_est+p1_est)
  p0.mean <- mean((1-z)*(s-p0_est)/(1-tp_est)+p0_est)
  w00_est <- (1 - p1_est) / (1 - p1.mean)
  w01_est <- (p1_est - p0_est)/(p1.mean - p0.mean)
  w11_est <- p0_est / p0.mean
  
  eps1 <- exp(xi1 * (times/max(times))^eta1)
  eps0 <- exp(xi0 * (times/max(times))^eta0)
  w1c <- sapply(eps1, function(x) (x * ps01_est + x * ps11_est) / (x * ps01_est + ps11_est))
  w0c <- sapply(eps0, function(x) (x * ps01_est + x * ps00_est) / (x * ps01_est + ps00_est))
  w0n <- sapply(eps0, function(x) (ps01_est + ps00_est) / (x * ps01_est + ps00_est))
  w1a <- sapply(eps1, function(x) (ps01_est + ps11_est) / (x * ps01_est + ps11_est))
  
  pos1 <- pmin(lower_bound_vec(times, U), n)
  pos2 <- pmax(upper_bound_vec(times, U) - 1, 1)
  time_id <- 1:length(times)
  
  subsets <- list(
    "00" = s == 0 & z == 0,
    "01" = s == 1 & z == 0,
    "10" = s == 0 & z == 1,
    "11" = s == 1 & z == 1
  )
  
  call_censor <- call(
    "~",
    as.call(list(
      quote(Surv), U_name, as.call(list(quote(`-`), 1, Delta_name))
    )),
    censor.model[[3]]
  )
  
  call_failure <- call(
    "~",
    as.call(list(
      quote(Surv), U_name, Delta_name
    )),
    failure.model[[3]]
  )
  
  calc_cumhaz <- function(ss) {
    
    formula_censor <- eval(call_censor)
    formula_failure <- eval(call_failure)
    
    fit <- list(
      T = coxph(formula_failure, data = data[ss, ]),
      C = coxph(formula_censor, data = data[ss, ])
    )
    
    new_outer <- function(a, b) {
      return (structure(list(
        a = a,
        b = b
      ), class = 'outer'))
    }
    
    `[.outer` <- function(x, i, j) {
      if (missing(i)) {
        return (x$a * x$b[j])
      }
      if (missing(j)) {
        return (x$a[i] * x$b)
      }
      return (x$a[i] * x$b[j])
    }
    
    tmp_func <- function(fit, X, delta) {
      #basehaz <- basehaz(fit, centered = F)
      #basehaz1 <- basehaz$hazard[pmax(upper_bound_vec(U, basehaz$time) - 1L, 1L)]
      #cumhaz <- outer(c(exp(X %*% fit$coefficients)), basehaz1)
      b_est <- fit$coefficients
      s1 <- cumsum(ss[n:1]*exp(as.matrix(X[n:1,]) %*% b_est))
      dls <- ss * delta / s1[n:1]
      dls[(1:length(dls))[is.nan(dls)]] <- 0
      cumhaz <- new_outer(c(exp(X %*% fit$coefficients)), cumsum(dls))
      return (cumhaz)
    }
    return (list(
      T = tmp_func(fit$T, XT, delta),
      C = tmp_func(fit$C, XC, 1 - delta)
    ))
  }
  
  calc_q <- function(cumhaz) {
    q <- diag(NA, n)
    q[,1] = -cumhaz$C[,1] * exp(cumhaz$C[,1] + cumhaz$T[,1])
    
    for(j in 2:n) {
      q[j:n,j] = q[j:n,j-1] - (cumhaz$C[j:n,j] - cumhaz$C[j:n,j-1]) * exp(cumhaz$C[j:n,j] + cumhaz$T[j:n,j])
    }
    for(i in 1:(n-1)){
      q[i,(i+1):n] = q[i,i] + (1 - delta[i]) * exp(cumhaz$C[i, i] + cumhaz$T[i, i])
    }
    return (q)
  }
  
  subset_ind <- c(
    "NACE" %in% estimand || "CACE" %in% estimand, # 00
    "AACE" %in% estimand, # 01
    "NACE" %in% estimand, # 10
    "CACE" %in% estimand || "AACE" %in% estimand # 11
  )
  
  subsets_to_use <- subsets[subset_ind]
  
  # calculate q only for needed
  
  cumhaz <- lapply(subsets_to_use, calc_cumhaz)
  q <- lapply(cumhaz, calc_q)
  gc()
  
  result <- list()
  
  gen_result <- function(func1, func0) {
    est1 <- mapply(function(x, y, t) func1(x, y, t) / n, pos1, pos2, time_id)
    est0 <- mapply(function(x, y, t) func0(x, y, t) / n, pos1, pos2, time_id)
    z1_mono = mono(est1)
    z0_mono = mono(est0)
    
    return (list(
      z1 = est1,
      z0 = est0,
      diff = est1 - est0,
      z1_mono = z1_mono,
      z0_mono = z0_mono,
      diff_mono = z1_mono - z0_mono
    ))
  }
  
  NACE_1 <- function(x, y, t) {
    sum(
      (1 - s[x:n]) * (z[x:n]) / ((1 - p1.mean) * (tp_est[x:n]) * exp(-cumhaz[["10"]]$C[x:n, y]))
    ) + sum(
      (1 - p1_est) / (1 - p1.mean) * (1 - z / tp_est) * exp(-cumhaz[["10"]]$T[, y])
    ) + sum(
      (1 - s) * z * exp(-cumhaz[["10"]]$T[, y]) * q[["10"]][, x] / ((1 - p1.mean) * tp_est)
    )
  }
  
  NACE_0 <- function(x, y, t) {
    s1 <- w00_est * (1 - s) * (1 - z)
    s2 <- (1 - p0_est) * (1 - tp_est)
    sum(
      w0n[x:n ,t] * s1[x:n] / (s2[x:n] * exp(-cumhaz[["00"]]$C[x:n, y]))
    ) + sum(
      w0n[, t] * w00_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["00"]]$T[, y])
    ) - sum(
      w0n[, t]^2 / eps0[t] * (z * (s - p1_est) / tp_est + (1 - p1_est) * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["00"]]$T[, y])
    ) / (1 - p1.mean) + sum(
      w0n[, t] * s1 * exp(-cumhaz[["00"]]$T[, y]) * q[["00"]][, x] / s2
    )
  }
  
  CACE_1 <- function(x, y, t) {
    s1 <- w01_est * s * z
    s2 <- p1_est * tp_est
    sum(
      w1c[x:n, t] * s1[x:n] / (s2[x:n] * exp(-cumhaz[["11"]]$C[x:n, y]))
    ) + sum(
      w1c[, t] * w01_est * (1 - z / tp_est) * exp(-cumhaz[["11"]]$T[, y])
    ) - sum(
      w1c[, t]^2 / eps1[t] *((1 - z) * (s - p0_est) / (1 - tp_est) + p0_est * z * (p1_est - s) / s2) * exp(-cumhaz[["11"]]$T[, y])
    ) / (p1.mean - p0.mean) + sum(
      w1c[, t] * s1 * exp(-cumhaz[["11"]]$T[, y]) * q[["11"]][, x] / s2
    )
  }
  
  CACE_0 <-  function(x, y, t) {
    s1 <- w01_est * (1 - s) * (1 - z)
    s2 <- (1 - p0_est) * (1 - tp_est)
    sum(
      w0c[x:n, t] * s1[x:n] / (s2[x:n] * exp(-cumhaz[["00"]]$C[x:n, y]))
    ) + sum(
      w0c[, t] * w01_est * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["00"]]$T[, y])
    ) + sum(
      w0c[, t]^2 / eps0[t] * (z * (s - p1_est) / tp_est + (1 - p1_est) * (1 - z) * (p0_est - s) / s2) * exp(-cumhaz[["00"]]$T[, y])
    ) / (p1.mean - p0.mean) + sum(
      w0c[, t] * s1 * exp(-cumhaz[["00"]]$T[, y]) * q[["00"]][, x] / s2
    )
  }
  
  AACE_1 <-  function(x, y, t) {
    s1 <- w11_est * s * z
    s2 <- p1_est * tp_est
    sum(
      w1a[x:n, t] * s1[x:n] / (s2[x:n] * exp(-cumhaz[["11"]]$C[x:n, y]))
    ) + sum(
      w1a[, t] * w11_est * (1 - z / tp_est) * exp(-cumhaz[["11"]]$T[, y])
    ) + sum(
      w1a[, t]^2 / eps1[t] * ((1 - z) * (s - p0_est) / (1 - tp_est) + p0_est * z * (p1_est - s) / s2) * exp(-cumhaz[["11"]]$T[, y])
    ) / (p0.mean) + sum(
      w1a[, t] * s1 * exp(-cumhaz[["11"]]$T[, y]) * q[["11"]][, x] / s2
    )
  }
  
  AACE_0 <- function(x, y, t) {
    sum(
      s[x:n] * (1 - z[x:n]) / (p0.mean * (1 - tp_est[x:n]) * exp(-cumhaz[["01"]]$C[x:n, y]))
    ) + sum(
      (p0_est) / (p0.mean) * (1 - (1 - z) / (1 - tp_est)) * exp(-cumhaz[["01"]]$T[, y])
    ) + sum(
      s * (1 - z) * exp(-cumhaz[["01"]]$T[, y]) * q[["01"]][, x] / (p0.mean * (1 - tp_est))
    )
  }
  
  if ("NACE" %in% estimand) 
    result$NACE <- gen_result(NACE_1, NACE_0)
  if ("CACE" %in% estimand)
    result$CACE <- gen_result(CACE_1, CACE_0)
  if ("AACE" %in% estimand)
    result$AACE <- gen_result(AACE_1, AACE_0)
  class(result) <- "est.MR.SA"
  return (result)
}

ACE.SA <- function(times, propensity.model, principal.model, censor.model, failure.model, data, 
                   xi0, xi1, eta0, eta1,
                   estimand = c("NACE", "CACE", "AACE"), bootstrap = 50) {
  mean_result <- est.MR.SA(times, propensity.model, principal.model, censor.model, failure.model, data, 
                           xi0, xi1, eta0, eta1, estimand)
  gc()
  if (bootstrap == 0)
    return (mean_result)
  bs_result <- c()
  for (name in estimand) {
    z1.B <- z0.B <- diff.B <- z1_mono.B <- z0_mono.B <- diff_mono.B <- matrix(nrow = length(times), ncol = bootstrap)
    bs_result[[name]] <- list(
      z1.B = z1.B,
      z0.B = z0.B,
      diff.B = diff.B,
      z1_mono.B = z1_mono.B,
      z0_mono.B = z0_mono.B,
      diff_mono.B = diff_mono.B
    )
  }
  pb <- progress::progress_bar$new(total = bootstrap, format = "(:spin) [:bar] :percent. :elapsed | :eta")
  for (i in 1:bootstrap) {
    gc()
    pb$tick()
    data.B <- slice_sample(data, n = nrow(data), replace = T)
    result <- est.MR.SA(times, propensity.model, principal.model, censor.model, failure.model, data.B, 
                        xi0, xi1, eta0, eta1, estimand)
    for (name in estimand) {
      bs_result[[name]]$z1.B[, i] <- result[[name]]$z1
      bs_result[[name]]$z0.B[, i] <- result[[name]]$z0
      bs_result[[name]]$diff.B[, i] <- result[[name]]$diff
      bs_result[[name]]$z1_mono.B[, i] <- result[[name]]$z1_mono
      bs_result[[name]]$z0_mono.B[, i] <- result[[name]]$z0_mono
      bs_result[[name]]$diff_mono.B[, i] <- result[[name]]$diff_mono
    }
    gc()
  }
  r_val <- list(
    estimate = mean_result,
    bootstrap = bs_result
  )
  class(r_val) <- "est.MR.SA.bootstrap"
  return (r_val)
}

plot.est.MR.SA.bootstrap <- function(res, time = c(0, 1)) {
  df_list <- list()
  for (name in c("NACE", "AACE", "CACE")) {
    df_list <- c(df_list, list(data.frame(
      time = time,
      estimand = name,
      group = "Z = 1",
      type = "Outcome",
      est = res$estimate[[name]][["z1_mono"]],
      lwr = apply(res$bootstrap[[name]][["z1_mono.B"]], 1, quantile,
                  0.025, na.rm = T),
      upr = apply(res$bootstrap[[name]][["z1_mono.B"]], 1, quantile, 
                  0.975, na.rm = T)
    )))
    df_list <- c(df_list, list(data.frame(
      time = time,
      estimand = name,
      group = "Z = 0",
      type = "Outcome",
      est = res$estimate[[name]][["z0_mono"]],
      lwr = apply(res$bootstrap[[name]][["z0_mono.B"]], 1, quantile, 
                  0.025, na.rm = T),
      upr = apply(res$bootstrap[[name]][["z0_mono.B"]], 1, quantile,
                  0.975, na.rm = T)
    )))
    df_list <- c(df_list, list(data.frame(
      time = time,
      estimand = name,
      group = "Difference",
      type = "Difference",
      est = res$estimate[[name]][["diff"]],
      lwr = apply(res$bootstrap[[name]][["diff.B"]], 1, quantile, 
                  0.025, na.rm = T),
      upr = apply(res$bootstrap[[name]][["diff.B"]], 1, quantile, 
                  0.975, na.rm = T)
    )))
  }
  df <- do.call(bind_rows, df_list)
  df$type <- factor(df$type, levels = c("Difference", "Outcome"))
  df$estimand <- ifelse(
    df$estimand == 'NACE', '00',
    ifelse(
      df$estimand == 'CACE', '01', '11'
    )
  )
  df$estimand <- factor(df$estimand, levels = c("00", "01", "11"))
  ggplot(df) + geom_line(aes(time, est, color = group)) +
    geom_ribbon(aes(time, ymin = lwr, ymax = upr, fill = group), alpha = 0.3) +
    facet_grid(type~estimand, scale = "free_y") + 
    ylab("Survival probability")
}
