library(PStrata)
library(dplyr)
library(ggplot2)

data <- readRDS("data.rds")
# Specify model formulas -----------------------------------------------------
S.formula <- arm + D ~ ages + SEX + nephrectomy + riskfac
Y.formula <- survtime + death ~ ages + SEX + nephrectomy + riskfac

# Fit principal stratification model  ------------------------
fit <- PStrata(
  S.formula         = S.formula,
  Y.formula         = Y.formula,
  Y.family          = survival(method = "Cox"),
  data              = data,
  strata            = c(a = "11", n = "00", c = "01"), # Monotonicity, no ER
  # strata            = c(a = "11*", n = "00", c = "01"), #  Monotonicity,  ER 
  # strata            = c(a = "11", n = "00", c = "01", d = "10"), # No monotonicity, no ER
  # strata            = c(a = "11*", n = "00", c = "01", d = "10"), # No monotonicity,  ER 
  prior_intercept   = prior_normal(0, 1),
  prior_coefficient = prior_normal(0, 1),
  warmup            = 2000,
  iter              = 4000,
  chains            = 3
)

# Summarize fitted model -----------------------------------------------------
summary(fit)   # stratum-specific summaries

# Posterior predictive survival probability curves --------------------------------------
fit_outcome <- PSOutcome(fit)
plot(fit_outcome) + 
  xlab("time") + 
  ylab("survival probability")

# Extract outcome array ------------------------------------------------------
outcome_array <- fit_outcome$outcome_array



#===========================================================
#=============================PI===========================
# Install and load package ---------------------------------------------------
devtools::install_github("chaochengstat/mrPStrata")
library(mrPStrata)

# Implement the multiply robust estimators------------------------------------
res = mrPStrata(times=seq(0,7,0.1),
                data = data,
                Xpi_names = c("ages","SEX","nephrectomy","riskfac"),
                Xe_names = c("ages","SEX","nephrectomy","riskfac"),
                Xc_names = c("ages","SEX","nephrectomy","riskfac"),
                Xt_names = c("ages","SEX","nephrectomy","riskfac"),
                Z_name = "arm",
                S_name = "D", 
                U_name ="survtime",
                delta_name = "death",  
                B=50)

# Extract results by stratum -------------------------------------------------
print(res$Compliers) 
print(res$Always_Takers)
print(res$Never_Takers)

# Visualize results ----------------------------------------------------------
plot.psce(res)

#Sensitivity analysis for the monotonicity------------------------------------
res_mo = mrPStrata_MO_SA(times=seq(0,7,0.1),
                         data = data,
                         Xpi_names = c("ages","SEX","nephrectomy","riskfac"),
                         Xe_names = c("ages","SEX","nephrectomy","riskfac"),
                         Xc_names = c("ages","SEX","nephrectomy","riskfac"),
                         Xt_names = c("ages","SEX","nephrectomy","riskfac"),
                         Z_name = "arm",
                         S_name = "D", 
                         U_name ="survtime",
                         delta_name = "death",  
                         zeta=0.2,
                         B=50)

print(res_mo$Compliers) 
print(res_mo$Always_Takers)
print(res_mo$Never_Takers)

# Visualize results ----------------------------------------------------------
plot.psce(res_mo)

#Sensitivity analysis for the principal ignorability--------------------------
res_PI = mrPStrata_PI_SA(times=seq(0,7,0.1),
                         data = data,
                         Xpi_names = c("ages","SEX","nephrectomy","riskfac"),
                         Xe_names = c("ages","SEX","nephrectomy","riskfac"),
                         Xc_names = c("ages","SEX","nephrectomy","riskfac"),
                         Xt_names = c("ages","SEX","nephrectomy","riskfac"),
                         Z_name = "arm",
                         S_name = "D", 
                         U_name ="survtime",
                         delta_name = "censortime",  
                         xi0 = log(0.9),
                         xi1 = log(0.9),
                         eta0=1,
                         eta1=1,
                         B=50)


