# Load necessary libraries
library(dplyr)
source("Simulation util.R")

# Set seed for reproducibility
set.seed(123)


generate_group_data <- function(n, sex_probs, ecog_probs, nephrectomy_prob, radiation_prob, 
                                age_median, age_range, metastases_probs, risk_probs) {
  
  # Sex (Male, Female)
  sex <- sample(c("Male", "Female"), size = n, replace = TRUE, prob = sex_probs)
  
  # Age (normally distributed around the median with a range)
  age <- round(rnorm(n, mean = age_median, sd = (age_range[2] - age_range[1]) / 4))
  
  # ECOG performance status (0, 1, 2)
  ecog <- sample(0:2, size = n, replace = TRUE, prob = ecog_probs)
  
  # Previous nephrectomy (Yes, No)
  previous_nephrectomy <- rbinom(n, 1, nephrectomy_prob)
  
  # Previous radiation therapy (Yes, No)
  previous_radiation <- rbinom(n, 1, radiation_prob)
  
  # Common sites of metastases
  metastases <- data.frame(
    Lung = rbinom(n, 1, metastases_probs[1]),
    Lymph_node = rbinom(n, 1, metastases_probs[2]),
    Bone = rbinom(n, 1, metastases_probs[3]),
    Liver = rbinom(n, 1, metastases_probs[4])
  )
  
  # Number of adverse risk factors (0, 1-2, >=3)
  risk_factors <- sample(0:2, size = n, replace = TRUE, prob = risk_probs)
  
  # Create a data frame
  data <- data.frame(Sex = sex, Age = age, ECOG_Performance_Status = ecog, 
                     Previous_Nephrectomy = previous_nephrectomy, 
                     Previous_Radiation = previous_radiation, 
                     Lung_Metastasis = metastases$Lung, 
                     Lymph_Node_Metastasis = metastases$Lymph_node, 
                     Bone_Metastasis = metastases$Bone, 
                     Liver_Metastasis = metastases$Liver, 
                     Risk_Factors = risk_factors)
  
  return(data)
}

#========Generate Data=========
n_bevacizumab <- 369
n_interferon <- 363

# Generate data for each treatment group
data_bevacizumab <- generate_group_data(
  n = 369, 
  sex_probs = c(0.73, 0.27), 
  ecog_probs = c(0.62, 0.36, 0.02), 
  nephrectomy_prob = c(0.85, 0.15), 
  radiation_prob = 0.09, 
  age_median = 61, age_range = c(56, 70), 
  metastases_probs = c(0.68, 0.35, 0.28, 0.20), 
  risk_probs = c(0.26, 0.64, 0.10))

data_interferon <- generate_group_data(
  n = 363,  
  sex_probs = c(0.66, 0.34),  
  ecog_probs = c(0.62, 0.37, 0.01),  
  nephrectomy_prob = c(0.85, 0.15),  
  radiation_prob = c(0.10, 0.90),  
  age_median = 62,  
  age_range = c(55, 70),  
  metastases_probs = c(0.70, 0.36, 0.30, 0.20),  
  risk_probs = c(0.26, 0.64, 0.10)  
)

# Combine both datasets and add a treatment group variable
combined_data <- bind_rows(
  mutate(data_bevacizumab, Group = "bevacizumab"),
  mutate(data_interferon, Group = "interferon")
)



# select_data <- data.frame(
#   X1 = combined_data$Sex == "Male",
#   X2 = combined_data$ECOG_Performance_Status == "0",
#   X3 = combined_data$Previous_Nephrectomy == "1",
#   X4 = combined_data$Previous_Radiation == "1",
#   X5 = combined_data$Lung_Metastasis == "1",
#   X6 = combined_data$Lymph_Node_Metastasis == "1",
#   X7 = combined_data$Bone_Metastasis == "1",
#   X8 = combined_data$Liver_Metastasis == "1",
#   X9 = combined_data$Risk_Factors == "1",
#   X10 = combined_data$Risk_Factors == "2",
#   X11 = (combined_data$Age - mean(combined_data$Age)) / sd(combined_data$Age)
# )

select_data <- data.frame(
  X1 = combined_data$Sex == "Male",
  X2 = combined_data$Previous_Nephrectomy == 1,
  X3 = combined_data$Risk_Factors == 2,
  X4 = (combined_data$Age - mean(combined_data$Age)) / sd(combined_data$Age)
)

covariates_data <- as.data.frame(lapply(select_data, as.numeric))

Z <- as.numeric(combined_data$Group == "bevacizumab")


# covariates_data <- as.data.frame(lapply(select_data, as.numeric))
# covariates_data_mis <- as.data.frame(lapply(select_data_mis, as.numeric))
#============Functions used in generating outcome ========================
gen_binary_data <- function(Z, data, strata_coef_matrix, outcome_coef_matrix,  true_model_id){
  data_matrix <- as.matrix(cbind(1, data))
  S <- sample_S(strata_coef_matrix, data_matrix, c(0, 1, 3))
  #table(S)
  # Z <- rbinom(N, 1, 0.5)
  D <- ifelse(bitwAnd(S, 2^(1-Z)), 1, 0)
  # print(c(sum(D[Z==0]), sum(D[Z==1])))  # our aim is to get (206, 168)
  Y_results <- sample_Y_binary(outcome_coef_matrix, data_matrix)
  Y_all <- Y_results[[1]]
  Y_all_prob <- Y_results[[2]]
  Y <- Y_all[, true_model_id[1]] * (S == 0 & Z == 0) + 
    Y_all[, true_model_id[2]] * (S == 0 & Z == 1) + 
    Y_all[, true_model_id[3]] * (S == 1 & Z == 0) +
    Y_all[, true_model_id[4]] * (S == 1 & Z == 1) +
    Y_all[, true_model_id[5]] * (S == 3 & Z == 0) + 
    Y_all[, true_model_id[6]] * (S == 3 & Z == 1) 
  sim_data <- data.frame(
    S, Z, D, Y, as.matrix(data)
  )
  true_strata_mean <- c(sum(Y_all[, 1] * (S == 0))/ sum(S==0),  sum(Y_all[, 2] * (S == 0))/ sum(S==0), 
                        sum(Y_all[, 3] * (S == 1))/ sum(S==1), sum(Y_all[, 4] * (S == 1))/ sum(S==1),
                        sum(Y_all[, 5] * (S == 3))/ sum(S==3), sum(Y_all[, 6] * (S == 3))/ sum(S==3))
  
  #print(c(sum(Y[Z==0]), sum(Y[Z==1])))
  return(list(
    data = sim_data,
    strata_coef_matrix = strata_coef_matrix,
    outcome_coef_matrix = outcome_coef_matrix,
    true_strata_mean = true_strata_mean
  ))
}

# phi is hazard for each strata
# gen_sur_data <- function(
    #     Z, data, strata_coef_matrix, outcome_coef_matrix, phi, 
#     true_model_id, beta_C
# ) {
#   set.seed(0)
#   data_matrix <- as.matrix(cbind(1, data))
#   S <- sample_S(strata_coef_matrix, data_matrix, c(0, 1, 3))
#   D <- ifelse(bitwAnd(S, 2^(1-Z)), 1, 0)
#   T_all <- sample_Y_cox(outcome_coef_matrix, phi, data_matrix)
#   T <- T_all[, true_model_id[1]] * (S == 0 & Z == 0) + 
#     T_all[, true_model_id[2]] * (S == 0 & Z == 1) + 
#     T_all[, true_model_id[3]] * (S == 1 & Z == 0) +
#     T_all[, true_model_id[4]] * (S == 1 & Z == 1) +
#     T_all[, true_model_id[5]] * (S == 3 & Z == 0) + 
#     T_all[, true_model_id[6]] * (S == 3 & Z == 1) +
#   C <- rexp(nrow(data_matrix), exp(data_matrix %*% beta_C))
#   Y <- pmin(T, C)
#   delta <- ifelse(T < C, 1, 0)
#   sim_data <- data.frame(
#     S, Z, D, T, C, Y, delta, as.matrix(data)
#   )
#   return (list(
#     data = sim_data,
#     strata_coef_matrix = strata_coef_matrix,
#     outcome_coef_matrix = outcome_coef_matrix,
#     phi = phi
#   ))
# }


gen_sur_data_eight <- function(
    Z, data, strata_coef_matrix, outcome_coef_matrix, phi, 
    true_model_id, beta_C
) {
  set.seed(0)
  data_matrix <- as.matrix(cbind(1, data))
  S <- sample_S(strata_coef_matrix, data_matrix, c(0, 1, 2, 3))
  D <- ifelse(bitwAnd(S, 2^(1-Z)), 1, 0)
  T_all <- sample_Y_cox(outcome_coef_matrix, phi, data_matrix)
  T <- T_all[, true_model_id[1]] * (S == 0 & Z == 0) + 
    T_all[, true_model_id[2]] * (S == 0 & Z == 1) + 
    T_all[, true_model_id[3]] * (S == 1 & Z == 0) +
    T_all[, true_model_id[4]] * (S == 1 & Z == 1) +
    T_all[, true_model_id[5]] * (S == 2 & Z == 0) + 
    T_all[, true_model_id[6]] * (S == 2 & Z == 1) +
    T_all[, true_model_id[7]] * (S == 3 & Z == 0) +
    T_all[, true_model_id[8]] * (S == 3 & Z == 1) 
  C <- rexp(nrow(data_matrix), exp(data_matrix %*% beta_C))
  Y <- pmin(T, C)
  delta <- ifelse(T < C, 1, 0)
  sim_data <- data.frame(
    S, Z, D, T, C, Y, delta, data
  )
  return (list(
    data = sim_data,
    strata_coef_matrix = strata_coef_matrix,
    outcome_coef_matrix = outcome_coef_matrix,
    phi = phi
  ))
}

gen_sur_PI_data <- function(
    Z, data, ice_coef, PI_outcome_coef_matrix, phi, 
    beta_C, true_model_id = c(1,2,3,4)
) {
  set.seed(0)
  data_matrix <- as.matrix(cbind(1, data))
  # S <- sample_S(strata_coef_matrix, data_matrix, c(0, 1, 3))
  # D <- ifelse(bitwAnd(S, 2^(1-Z)), 1, 0)
  D <- sample_D(ice_coef, data_matrix, Z) # ICE
  T_all <- sample_Y_cox(PI_outcome_coef_matrix, phi, data_matrix)
  T <- T_all[, true_model_id[1]] * (D == 0 & Z == 0) + 
    T_all[, true_model_id[2]] * (D == 1 & Z == 0) + 
    T_all[, true_model_id[3]] * (D == 0 & Z == 1) +
    T_all[, true_model_id[4]] * (D == 1 & Z == 1) 
  C <- rexp(nrow(data_matrix), exp(data_matrix %*% beta_C))
  Y <- pmin(T, C)
  delta <- ifelse(T < C, 1, 0)
  sim_data <- data.frame(
    Z, D, T, C, Y, delta, as.matrix(data)
  )   # No longer return S as we don't know.
  return (list(
    data = sim_data,
    PI_outcome_coef_matrix = PI_outcome_coef_matrix,
    phi = phi
  ))
}

##==============For whole dataset===========

# four_strata_coef_matrix = t(matrix(
#   c(
#     0, 0, 0, 0, 0,               # U = (1, 1) 0, -0.1, -0.1, -0.1, 0
#     -1.4, 0.1, 0.2, 0.3, 0.5,    # U = (0, 1) -1.4, 0, 0.1, 0.2, 0.5
#     -2.3, 0, 0, 0, 0.5,          # U = (1, 0) -2.3, -0.1, -0.1, -0.1, 0.5
#     0, 0.1, 0.1, 0.1, 0          # U = (0, 0) 0, 0, 0, 0, 0
#   ), ncol = 5, byrow = TRUE
# ))

four_strata_coef_matrix = t(matrix(
  c(
     0, -0.1, -0.1, -0.1, 0,
     -1.4, 0, 0.1, 0.2, 0.5,
    -2.3, -0.1, -0.1, -0.1, 0.5,
    0, 0, 0, 0, 0
  ), ncol = 5, byrow = TRUE
))


outcome_coef_eight = t(matrix(
  c(
    -5, -0.3, 0.3,  1, 1,      # (1, 1) Z = 0  (ER)   (PI2)  (G1)
    -5.2, -0.3, 0.3,  1, 1,    # (1, 1) Z = 1  (ER)          (G1)
    -4.8, -0.3, 0.3, 1, 1,     # (0, 1) Z = 0         (PI2)  (G1)
    -5.5, -0.5, 0.6, 1.5, 1.5,     # (0, 1) Z = 1         (PI1)  (G2)
    -4, -0.5, 0.6, 0.5, 1.0,   # (1, 0) Z = 0  (mono) 
    -4.5, -0.5, 0.6, 0.5, 0.5, # (1, 0) Z = 1  (mono) 
    -4.5, -0.6, 0.5, 0.9, 0.5,   # (0, 0) Z = 0 !!!            (G3)
    -5.5, -0.5, 0.6, 1.4, 1.5      # (0, 0) Z = 1 !!!     (PI1)  (G2)
  ), ncol = 5, byrow = TRUE
))

beta_C_eight <- c(-5, rep(-1, 4))  # eight means no ER, no monotonicity
phi_eight = c(2, 2, 2, 1.6, 2.5, 2.5, 2.2, 1.6)
#phi_ER = c(10,10,10,10)
true_model_id_eight <- c(1, 2,3,4,5,6,7,8)

sim_noER_nomo_result <- gen_sur_data_eight(Z, covariates_data, four_strata_coef_matrix, outcome_coef_eight, phi_eight, true_model_id_eight, beta_C_eight)
data_noER_nomo <- sim_noER_nomo_result$data

print_info <- function(df){
  print(c(nrow(df), sum(df$Z==0), sum(df$Z==1)))
  print(table(df$S))
  print(table(df$S)/nrow(df))
  print("Non-complicane -- D:")
  print(c(sum(df$D==1), sum(df$D[df$Z==0]==1), sum(df$D[df$Z==1]==1)))
  print("Event -- delta:")
  print(c(sum(df$delta==1), sum(df$delta[df$Z==0]==1), sum(df$delta[df$Z==1]==1)))
  print(c(sum(df$delta==1)/nrow(df) , sum(df$delta[df$Z==0]==1)/sum(df$Z==0), sum(df$delta[df$Z==1]==1)/sum(df$Z==1)))
  print("Censoring Rate:")
  print(c(sum(df$delta==0), sum(df$delta[df$Z==0]==0), sum(df$delta[df$Z==1]==0)))
  print(c(sum(df$delta==0)/nrow(df), sum(df$delta[df$Z==0]==0)/sum(df$Z==0), sum(df$delta[df$Z==1]==0)/sum(df$Z==1)))
  print("T:")
  print(summary(df$T))
}

print_info(data_noER_nomo) 

data <- data_noER_nomo %>%
    rename(
      arm        = Z,
      ages       = X4,
      SEX        = X1,
      nephrectomy = X2,
      riskfac    = X3,
      survtime   = Y,
      death      = delta
    )

saveRDS(data, file = "data.rds")
