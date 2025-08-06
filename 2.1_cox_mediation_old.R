# Mediation Analysis

# Fit Cox models first for the total effect of sdoh_score on MORTSTAT?

# Move directly to mediation analysis (e.g. using mediation package or lavaan)?
  
#  Do both, starting with total effect?



# 1 Setup: Packages and Directories-------

# Set working directory to analysis folder
setwd("/Users/dengshuyue/Desktop/SDOH/analysis")

# Define directory structure
dir <- list()
dir$root    <- getwd()
dir$data    <- file.path(dir$root, "data")
dir$output  <- file.path(dir$root, "output")
dir$code    <- file.path(dir$root, "code")


# Load Required Packages------
# List of required packages
want <- c("dplyr", "survey", "foreign", "Hmisc", "data.table", "tidyr", 
          "tibble", "readr", "flextable", "officer", "usethis", "gert", "survival")

# Install any missing packages
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)

# Load all required packages
lapply(want, function(pkg) require(pkg, character.only = TRUE))
rm(want, need)

# 1.0 Load merged dataset -----
df <- fread(file.path(dir$data, "SODH_diet_mort4.csv"))

# Check SDOH score and mortality variables
names(df)

df$FI <- ifelse(df$FS == 0, 1, 0)  # 1 = food insecure, 0 = food secure

# 1.1 Create Survey Design Object-----
# Create survey design
nhanes_design <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wt10,
  nest = TRUE,
  data = df
)

names(df)



# 1.2 cov list ------
# Covariates for Cox regression and mediation analysis
covariates_base <- c(
  "RIDAGEYR",       # Age
  "SEX",            # Sex
  "RACE",           # Race/ethnicity
  "household_size" # Household size
 # "EDU",
 # "SDDSRVYR"        # NHANES cycle
)

covariates_behavior <- c(
  "SMK_AVG",        # Avg cigarettes per day
  "SMK",            # Former smoker 
  "ALCG2",          # Alcohol use category
  "met_hr"          # Physical activity (quintile)
)

covariates_clinical <- c(
  "bmic",           # BMI category
  "DIABE",          # Diabetes
  "HYPERTEN",       # Hypertension
  "chol_rx",        # High cholesterol
  "CVD",            # Cardiovascular disease
  "cancer"          # Cancer history
)

# add behavior cov
cov_base_behavior <- c(covariates_base, covariates_behavior)

# Combine for full adjustment
covariates_full <- c(covariates_base, covariates_behavior, covariates_clinical)


# 2. Without mediation, SDOH predict moretality-----

# 2.1 Base model: SDOH score + demographic covariates ------
cox_sdoh_base <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                   paste(covariates_base, collapse = " + "))),
  design = nhanes_design
)
summary(cox_sdoh_base)

# 2.2 Extended model: SDOH score + demographic + lifestyle covariates-------
# Combine covariates for clarity
cov_base_behavior <- c(covariates_base, covariates_behavior)

cox_sdoh_base_behavior <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                   paste(cov_base_behavior, collapse = " + "))),
  design = nhanes_design
)
summary(cox_sdoh_base_behavior)

# 2.3 Fully adjusted model: SDOH score + demographic + lifestyle + clinical covariates ----
cox_sdoh_full <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)
summary(cox_sdoh_full)


# 3. Analysis based on Hypo -----
# 3.1 H1 model----- 
# FI → mortality, independent of SDOH, diet quality, and other covariates

fs_sdoh_diet_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score + ahei_total +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)

# fs_mt_depr <- svycoxph(
#  as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score + ahei_total + probable_depression +", 
#                   paste(covariates_full, collapse = " + "))),
#  design = nhanes_design
# )

summary(fs_sdoh_diet_mt)

# report in missing 
summary(df$ahei_total)
summary(df$sdoh_score)
summary(df$probable_depression)

# gut check not reported
df %>%
  group_by(SDDSRVYR) %>%
  summarise(
    total = n(),
    missing = sum(is.na(probable_depression)),
    missing_pct = round(mean(is.na(probable_depression)) * 100, 1)
  )


# 3.2 H2 model same model as H1 ------
# 2. Diet quality →  mortality, independent of SDOH, FS, and other covariates.

FI_sdoh_diet_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score + ahei_total +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)
summary(FI_sdoh_diet_mt)


# 3.3 H3: SDOH → FI or AHEI → Mortality (Separate Mediation) -----
# hypo 3. Food insecurity, diet quality separately mediate SDOH  → mortality.

###### model separately with only full cov list -----
# Step 1: Total effect model (SDOH → Mortality)
sdoh_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)


sdoh_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)

summary(sdoh_mt)

# Step 2: SDOH + FI → Mortality
sdoh_FI_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FI +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)

summary(sdoh_FI_mt)

# Step 3: SDOH + AHEI → Mortality
sdoh_ahei_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + ahei_total +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)

summary(sdoh_ahei_mt)

# Step 4: Extract SDOH coefficients
beta_total    <- coef(sdoh_mt)["sdoh_score"]
beta_FI       <- coef(sdoh_FI_mt)["sdoh_score"]
beta_ahei     <- coef(sdoh_ahei_mt)["sdoh_score"]

# Step 5: Calculate proportion mediated
prop_FI   <- (beta_total - beta_FI) / beta_total
prop_ahei <- (beta_total - beta_ahei) / beta_total

# Step 6: Display results
cat("Proportion of SDOH effect on mortality mediated by FI:   ", round(prop_FI * 100, 2), "%\n")
cat("Proportion of SDOH effect on mortality mediated by AHEI: ", round(prop_ahei * 100, 2), "%\n")

# Optional: Save for summary table or export
mediation_results <- list(
  beta_total = beta_total,
  beta_FI = beta_FI,
  prop_FI = prop_FI,
  beta_ahei = beta_ahei,
  prop_ahei = prop_ahei
)


###### (not use) gut check try to replicate what lu did ------
# df$FI <- ifelse(df$FS == 0, 1, 0)  # 1 = food insecure, 0 = food secure

# covariates_base_lu <- c("RIDAGEYR", "SEX", "RACE", "EDU")

# formula_FI_mort <- as.formula(
#  paste("Surv(py, MORTSTAT) ~ FI +", paste(covariates_base_lu, collapse = " + "))
#)


# cox_fs_mort <- svycoxph(formula_fs_mort, design = nhanes_design)

# Step 4: View results
# summary(cox_fs_mort)

##### create a function for cov set ------

# Define updated function to run mediation + Wald test
run_mediation_analysis <- function(covariate_list, design_obj) {
  # Step 1: Total effect model
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Step 2: SDOH + FI
  model_FI <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FI +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Step 3: SDOH + AHEI
  model_ahei <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + ahei_total +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Extract beta and SEs
  beta_total <- coef(model_total)["sdoh_score"]
  se_total <- sqrt(vcov(model_total)["sdoh_score", "sdoh_score"])
  
  beta_FI <- coef(model_FI)["sdoh_score"]
  se_FI <- sqrt(vcov(model_FI)["sdoh_score", "sdoh_score"])
  
  beta_ahei <- coef(model_ahei)["sdoh_score"]
  se_ahei <- sqrt(vcov(model_ahei)["sdoh_score", "sdoh_score"])
  
  # Proportion mediated
  prop_FI <- (beta_total - beta_FI) / beta_total
  prop_ahei <- (beta_total - beta_ahei) / beta_total
  
  # Wald test for FI mediation
  z_FI <- (beta_total - beta_FI) / sqrt(se_total^2 + se_FI^2)
  pval_FI <- 2 * (1 - pnorm(abs(z_FI)))
  
  # Wald test for AHEI mediation
  z_ahei <- (beta_total - beta_ahei) / sqrt(se_total^2 + se_ahei^2)
  pval_ahei <- 2 * (1 - pnorm(abs(z_ahei)))
  
  # Return results
  return(list(
    beta_total = beta_total,
    beta_FI = beta_FI,
    se_FI = se_FI,
    prop_FI = prop_FI,
    pval_FI = pval_FI,
    beta_ahei = beta_ahei,
    se_ahei = se_ahei,
    prop_ahei = prop_ahei,
    pval_ahei = pval_ahei
  ))
}


# Run for each covariate set
result_base <- run_mediation_analysis(covariates_base, nhanes_design)
result_behavior <- run_mediation_analysis(cov_base_behavior, nhanes_design)
result_full <- run_mediation_analysis(covariates_full, nhanes_design)



# Combine into a table with HRs and p-values
mediation_table <- data.frame(
  Covariate_Set = c("Base", "Base + Behavior", "Full"),
  
  HR_Total    = round(exp(c(result_base$beta_total,
                            result_behavior$beta_total,
                            result_full$beta_total)), 3),
  
  HR_FI       = round(exp(c(result_base$beta_FI,
                            result_behavior$beta_FI,
                            result_full$beta_FI)), 3),
  
  Prop_FI     = round(c(result_base$prop_FI,
                        result_behavior$prop_FI,
                        result_full$prop_FI) * 100, 2),
  
  P_FI        = signif(c(result_base$pval_FI,
                         result_behavior$pval_FI,
                         result_full$pval_FI), 3),
  
  HR_AHEI     = round(exp(c(result_base$beta_ahei,
                            result_behavior$beta_ahei,
                            result_full$beta_ahei)), 3),
  
  Prop_AHEI   = round(c(result_base$prop_ahei,
                        result_behavior$prop_ahei,
                        result_full$prop_ahei) * 100, 2),
  
  P_AHEI      = signif(c(result_base$pval_ahei,
                         result_behavior$pval_ahei,
                         result_full$pval_ahei), 3)
)

# Print final table
print(mediation_table)


# 4. Food insecurity and diet quality together mediate SDOH  → mortality ------ 
# (i.e., tested in the same model).


# Function to run joint mediation analysis (FI + AHEI together)
run_joint_mediation <- function(covariate_list, design_obj) {
  # Step 1: Total effect model
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Step 2: Joint mediation model (SDOH + FI + AHEI)
  model_joint <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FI + ahei_total +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Extract beta and SEs
  beta_total <- coef(model_total)["sdoh_score"]
  beta_joint <- coef(model_joint)["sdoh_score"]
  
  se_total <- sqrt(vcov(model_total)["sdoh_score", "sdoh_score"])
  se_joint <- sqrt(vcov(model_joint)["sdoh_score", "sdoh_score"])
  
  # Proportion mediated
  prop_joint <- (beta_total - beta_joint) / beta_total
  
  # Wald test
  z_joint <- (beta_total - beta_joint) / sqrt(se_total^2 + se_joint^2)
  pval_joint <- 2 * (1 - pnorm(abs(z_joint)))
  
  # Return results
  return(list(
    beta_total = beta_total,
    beta_joint = beta_joint,
    HR_total = exp(beta_total),
    HR_joint = exp(beta_joint),
    prop_joint = prop_joint,
    pval_joint = pval_joint
  ))
}


# Run the joint mediation for each covariate set
result_joint_base <- run_joint_mediation(covariates_base, nhanes_design)
result_joint_behavior <- run_joint_mediation(cov_base_behavior, nhanes_design)
result_joint_full <- run_joint_mediation(covariates_full, nhanes_design)

# Create summary table
joint_mediation_table <- data.frame(
  Covariate_Set = c("Base", "Base + Behavior", "Full"),
  HR_Total = round(c(result_joint_base$HR_total,
                     result_joint_behavior$HR_total,
                     result_joint_full$HR_total), 3),
  HR_Joint = round(c(result_joint_base$HR_joint,
                     result_joint_behavior$HR_joint,
                     result_joint_full$HR_joint), 3),
  Prop_Joint = round(c(result_joint_base$prop_joint,
                       result_joint_behavior$prop_joint,
                       result_joint_full$prop_joint) * 100, 2),
  P_Joint = signif(c(result_joint_base$pval_joint,
                     result_joint_behavior$pval_joint,
                     result_joint_full$pval_joint), 3)
)

print(joint_mediation_table)


# 5. Diet quality mediates  food insecurity → mortality,  adjusting for SDOH + covariates.-------
run_FI_mediation_by_ahei <- function(covariate_list, design_obj) {
  # Total effect: FI on mortality (adjusted for SDOH + covariates)
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Adjusted model: FI + AHEI + SDOH + covariates
  model_mediator <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + ahei_total + sdoh_score +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Extract beta and SEs for FI
  beta_total <- coef(model_total)["FI"]
  beta_mediator <- coef(model_mediator)["FI"]
  
  se_total <- sqrt(vcov(model_total)["FI", "FI"])
  se_mediator <- sqrt(vcov(model_mediator)["FI", "FI"])
  
  # Proportion mediated
  prop_mediated <- (beta_total - beta_mediator) / beta_total
  
  # Wald test
  z <- (beta_total - beta_mediator) / sqrt(se_total^2 + se_mediator^2)
  pval <- 2 * (1 - pnorm(abs(z)))
  
  # Return results
  return(list(
    beta_total = beta_total,
    beta_mediator = beta_mediator,
    HR_total = exp(beta_total),
    HR_mediator = exp(beta_mediator),
    prop_mediated = prop_mediated,
    pval = pval
  ))
}

result_FI_base <- run_FI_mediation_by_ahei(covariates_base, nhanes_design)
result_FI_behavior <- run_FI_mediation_by_ahei(cov_base_behavior, nhanes_design)
result_FI_full <- run_FI_mediation_by_ahei(covariates_full, nhanes_design)

# Create results table
FI_ahei_mediation_table <- data.frame(
  Covariate_Set = c("Base", "Base + Behavior", "Full"),
  HR_FI_Total = round(c(result_FI_base$HR_total,
                        result_FI_behavior$HR_total,
                        result_FI_full$HR_total), 3),
  HR_FI_Adjusted_AHEI = round(c(result_FI_base$HR_mediator,
                                result_FI_behavior$HR_mediator,
                                result_FI_full$HR_mediator), 3),
  Prop_Mediated = round(c(result_FI_base$prop_mediated,
                          result_FI_behavior$prop_mediated,
                          result_FI_full$prop_mediated) * 100, 2),
  P_Value = signif(c(result_FI_base$pval,
                     result_FI_behavior$pval,
                     result_FI_full$pval), 3)
)

print(FI_ahei_mediation_table)



# H6. Exploratory: depressive mediates food insecurity → mortality, adjusting for SDOH + covariates.-----

run_FI_mediation_by_depression <- function(covariate_list, design_obj) {
  # Total effect model: FI → Mortality (adjusted for SDOH + covariates)
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Mediator-adjusted model: FI → Mortality, adjusted for depression
  model_mediator <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + probable_depression + sdoh_score +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Extract coefficients and SEs for FI
  beta_total <- coef(model_total)["FI"]
  beta_mediator <- coef(model_mediator)["FI"]
  
  se_total <- sqrt(vcov(model_total)["FI", "FI"])
  se_mediator <- sqrt(vcov(model_mediator)["FI", "FI"])
  
  # Proportion mediated
  prop_mediated <- (beta_total - beta_mediator) / beta_total
  
  # Wald test
  z <- (beta_total - beta_mediator) / sqrt(se_total^2 + se_mediator^2)
  pval <- 2 * (1 - pnorm(abs(z)))
  
  return(list(
    beta_total = beta_total,
    beta_mediator = beta_mediator,
    HR_total = exp(beta_total),
    HR_mediator = exp(beta_mediator),
    prop_mediated = prop_mediated,
    pval = pval
  ))
}



result_dep_base <- run_FI_mediation_by_depression(covariates_base, nhanes_design)
result_dep_behavior <- run_FI_mediation_by_depression(cov_base_behavior, nhanes_design)
result_dep_full <- run_FI_mediation_by_depression(covariates_full, nhanes_design)

# Create summary table
depression_mediation_table <- data.frame(
  Covariate_Set = c("Base", "Base + Behavior", "Full"),
  HR_FI_Total = round(c(result_dep_base$HR_total,
                        result_dep_behavior$HR_total,
                        result_dep_full$HR_total), 3),
  HR_FI_Adjusted_Dep = round(c(result_dep_base$HR_mediator,
                               result_dep_behavior$HR_mediator,
                               result_dep_full$HR_mediator), 3),
  Prop_Mediated = round(c(result_dep_base$prop_mediated,
                          result_dep_behavior$prop_mediated,
                          result_dep_full$prop_mediated) * 100, 2),
  P_Value = signif(c(result_dep_base$pval,
                     result_dep_behavior$pval,
                     result_dep_full$pval), 3)
)

print(depression_mediation_table)



# 7 summary -----
# Build long-form table across hypotheses and covariate sets
combined_mediation_long <- data.frame(
  Hypothesis = rep(c(
    "H3a: FI mediates SDOH → Mortality",
    "H3b: AHEI mediates SDOH → Mortality",
    "H4: FI + AHEI jointly mediate SDOH → Mortality",
    "H5: AHEI mediates FI → Mortality",
    "H6: Depression mediates FI → Mortality"
  ), each = 3),
  
  Covariate_Set = rep(c("Base", "Base + Behavior", "Full"), times = 5),
  
  HR_Total = round(c(
    exp(result_base$beta_total),       # H3a
    exp(result_behavior$beta_total),
    exp(result_full$beta_total),
    
    exp(result_base$beta_total),       # H3b
    exp(result_behavior$beta_total),
    exp(result_full$beta_total),
    
    result_joint_base$HR_total,        # H4
    result_joint_behavior$HR_total,
    result_joint_full$HR_total,
    
    result_FI_base$HR_total,           # H5
    result_FI_behavior$HR_total,
    result_FI_full$HR_total,
    
    result_dep_base$HR_total,          # H6
    result_dep_behavior$HR_total,
    result_dep_full$HR_total
  ), 3),
  
  HR_Adjusted = round(c(
    exp(result_base$beta_FI),          # H3a
    exp(result_behavior$beta_FI),
    exp(result_full$beta_FI),
    
    exp(result_base$beta_ahei),        # H3b
    exp(result_behavior$beta_ahei),
    exp(result_full$beta_ahei),
    
    result_joint_base$HR_joint,        # H4
    result_joint_behavior$HR_joint,
    result_joint_full$HR_joint,
    
    result_FI_base$HR_mediator,        # H5
    result_FI_behavior$HR_mediator,
    result_FI_full$HR_mediator,
    
    result_dep_base$HR_mediator,       # H6
    result_dep_behavior$HR_mediator,
    result_dep_full$HR_mediator
  ), 3),
  
  Prop_Mediated = round(c(
    result_base$prop_FI,               # H3a
    result_behavior$prop_FI,
    result_full$prop_FI,
    
    result_base$prop_ahei,             # H3b
    result_behavior$prop_ahei,
    result_full$prop_ahei,
    
    result_joint_base$prop_joint,      # H4
    result_joint_behavior$prop_joint,
    result_joint_full$prop_joint,
    
    result_FI_base$prop_mediated,      # H5
    result_FI_behavior$prop_mediated,
    result_FI_full$prop_mediated,
    
    result_dep_base$prop_mediated,     # H6
    result_dep_behavior$prop_mediated,
    result_dep_full$prop_mediated
  ) * 100, 2),
  
  P_Value = signif(c(
    result_base$pval_FI,               # H3a
    result_behavior$pval_FI,
    result_full$pval_FI,
    
    result_base$pval_ahei,             # H3b
    result_behavior$pval_ahei,
    result_full$pval_ahei,
    
    result_joint_base$pval_joint,      # H4
    result_joint_behavior$pval_joint,
    result_joint_full$pval_joint,
    
    result_FI_base$pval,               # H5
    result_FI_behavior$pval,
    result_FI_full$pval,
    
    result_dep_base$pval,              # H6
    result_dep_behavior$pval,
    result_dep_full$pval
  ), 3)
)

# Show full table
print(combined_mediation_long)




