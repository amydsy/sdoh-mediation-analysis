# Mediation Analysis

# other analysis might explore later if have time
# Move directly to mediation analysis (e.g. using mediation package or lavaan)?

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
want <- c("dplyr", "survey", "foreign", "Hmisc", "data.table", "tidyr", "stringr",
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

# ---- FUNCTION 1: Separate Mediation (FI and AHEI) ----

run_separate_mediation <- function(covariate_list, design_obj) {
  sdoh_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  sdoh_FI_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FI +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  sdoh_AHEI_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + ahei_total +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  beta_total <- coef(sdoh_mort)["sdoh_score"]
  se_total <- sqrt(vcov(sdoh_mort)["sdoh_score", "sdoh_score"])
  
  beta_FI <- coef(sdoh_FI_mort)["sdoh_score"]
  se_FI <- sqrt(vcov(sdoh_FI_mort)["sdoh_score", "sdoh_score"])
  
  beta_AHEI <- coef(sdoh_AHEI_mort)["sdoh_score"]
  se_AHEI <- sqrt(vcov(sdoh_AHEI_mort)["sdoh_score", "sdoh_score"])
  
  prop_FI <- (beta_total - beta_FI) / beta_total
  prop_AHEI <- (beta_total - beta_AHEI) / beta_total
  
  z_FI <- (beta_total - beta_FI) / sqrt(se_total^2 + se_FI^2)
  z_AHEI <- (beta_total - beta_AHEI) / sqrt(se_total^2 + se_AHEI^2)
  
  pval_FI <- 2 * (1 - pnorm(abs(z_FI)))
  pval_AHEI <- 2 * (1 - pnorm(abs(z_AHEI)))
  
  return(list(
    beta_total = beta_total,
    se_total = se_total,
    beta_FI = beta_FI,
    se_FI = se_FI,
    prop_FI = prop_FI,
    pval_FI = pval_FI,
    beta_AHEI = beta_AHEI,
    se_AHEI = se_AHEI,
    prop_AHEI = prop_AHEI,
    pval_AHEI = pval_AHEI
  ))
}


# 4. Food insecurity and diet quality together mediate SDOH  → mortality ------ 
# (i.e., tested in the same model).

# ---- FUNCTION 2: Joint Mediation (FI + AHEI) ----
run_joint_mediation <- function(covariate_list, design_obj) {
  # Step 1: Total effect model
  sdoh_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Step 2: Joint mediator model
  sdoh_FI_AHEI_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FI + ahei_total +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Extract betas and SEs
  beta_total <- coef(sdoh_mort)["sdoh_score"]
  beta_joint <- coef(sdoh_FI_AHEI_mort)["sdoh_score"]
  
  se_total <- sqrt(vcov(sdoh_mort)["sdoh_score", "sdoh_score"])
  se_joint <- sqrt(vcov(sdoh_FI_AHEI_mort)["sdoh_score", "sdoh_score"])
  
  # Compute mediated proportion and p-value
  prop_joint <- (beta_total - beta_joint) / beta_total
  z <- (beta_total - beta_joint) / sqrt(se_total^2 + se_joint^2)
  pval <- 2 * (1 - pnorm(abs(z)))
  
  # Optional: assign models to global environment (for inspection)
  # assign("sdoh_mort", sdoh_mort, envir = .GlobalEnv)
  # assign("sdoh_FI_AHEI_mort", sdoh_FI_AHEI_mort, envir = .GlobalEnv)
  
  return(list(
    beta_total = beta_total,
    se_total = se_total,        # ✅ added
    beta_joint = beta_joint,
    se_joint = se_joint,        # ✅ added
    HR_total = exp(beta_total),
    HR_joint = exp(beta_joint),
    prop_joint = prop_joint,
    pval_joint = pval
  ))
}



# 5. Diet quality mediates  food insecurity → mortality,  adjusting for SDOH + covariates.-------
# H5

# ---- FUNCTION 3: AHEI mediates FI → Mortality ----

run_FI_mediation_by_ahei <- function(covariate_list, design_obj) {
  # Step 1: Total effect of FI on mortality
  FI_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Step 2: Mediation by AHEI
  FI_AHEI_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + ahei_total + sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Extract betas and SEs for FI
  beta_total <- coef(FI_mort)["FI"]
  beta_mediator <- coef(FI_AHEI_mort)["FI"]
  
  se_total <- sqrt(vcov(FI_mort)["FI", "FI"])
  se_mediator <- sqrt(vcov(FI_AHEI_mort)["FI", "FI"])
  
  # Calculate proportion mediated and p-value
  prop <- (beta_total - beta_mediator) / beta_total
  z <- (beta_total - beta_mediator) / sqrt(se_total^2 + se_mediator^2)
  pval <- 2 * (1 - pnorm(abs(z)))
  
  # Optional: expose model objects to global environment for debugging
  # assign("FI_mort", FI_mort, envir = .GlobalEnv)
  # assign("FI_AHEI_mort", FI_AHEI_mort, envir = .GlobalEnv)
  
  return(list(
    beta_total = beta_total,
    se_total = se_total,             # ✅ added
    beta_mediator = beta_mediator,
    se_mediator = se_mediator,       # ✅ added
    HR_total = exp(beta_total),
    HR_mediator = exp(beta_mediator),
    prop_mediated = prop,
    pval = pval
  ))
}


# 6. Exploratory: depressive mediates food insecurity → mortality, adjusting for SDOH + covariates.-----
# H6
# ---- FUNCTION 4: Depression mediates FI → Mortality (complete-case version) -----

run_FI_mediation_by_depression <- function(covariate_list, design_obj) {
  
  # Step 1: Identify variables used in both models
  vars_needed <- c("FI", "sdoh_score", "probable_depression", covariate_list, "MORTSTAT", "py")
  
  # Step 2: Subset to complete cases
  complete_cases <- design_obj$variables[complete.cases(design_obj$variables[, ..vars_needed]), ]
  
  # Step 3: Subset survey design to matched sample
  design_cc <- subset(design_obj, SEQN %in% complete_cases$SEQN)
  
  # Step 4: Fit models on matched sample
  FI_sdoh_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_cc
  )
  
  FI_depression_mort <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + probable_depression + sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_cc
  )
  
  # Step 5: Extract estimates and SEs
  beta_total <- coef(FI_sdoh_mort)["FI"]
  beta_mediator <- coef(FI_depression_mort)["FI"]
  
  se_total <- sqrt(vcov(FI_sdoh_mort)["FI", "FI"])
  se_mediator <- sqrt(vcov(FI_depression_mort)["FI", "FI"])
  
  prop <- (beta_total - beta_mediator) / beta_total
  z <- (beta_total - beta_mediator) / sqrt(se_total^2 + se_mediator^2)
  pval <- 2 * (1 - pnorm(abs(z)))
  
  # Optional: Assign models to global env for verification
  # assign("FI_sdoh_mort", FI_sdoh_mort, envir = .GlobalEnv)
  # assign("FI_depression_mort", FI_depression_mort, envir = .GlobalEnv)
  
  return(list(
    beta_total = beta_total,
    se_total = se_total,             # ✅ added
    beta_mediator = beta_mediator,
    se_mediator = se_mediator,       # ✅ added
    HR_total = exp(beta_total),
    HR_mediator = exp(beta_mediator),
    prop_mediated = prop,
    pval = pval
  ))
}


# 7.1 Run Models for All Covariate Sets -----
# Run all hypotheses with 3 covariate sets

# H3a & H3b: SDOH → Mortality, mediated by FI and AHEI separately
result_sdoh_base <- run_separate_mediation(covariates_base, nhanes_design)
result_sdoh_behavior <- run_separate_mediation(cov_base_behavior, nhanes_design)
result_sdoh_full <- run_separate_mediation(covariates_full, nhanes_design)

# H4: Joint mediation by FI + AHEI
result_joint_base <- run_joint_mediation(covariates_base, nhanes_design)
result_joint_behavior <- run_joint_mediation(cov_base_behavior, nhanes_design)
result_joint_full <- run_joint_mediation(covariates_full, nhanes_design)

# H5: AHEI mediates FI → Mortality
result_FI_base <- run_FI_mediation_by_ahei(covariates_base, nhanes_design)
result_FI_behavior <- run_FI_mediation_by_ahei(cov_base_behavior, nhanes_design)
result_FI_full <- run_FI_mediation_by_ahei(covariates_full, nhanes_design)

# H6: Depression mediates FI → Mortality
result_dep_base <- run_FI_mediation_by_depression(covariates_base, nhanes_design)
result_dep_behavior <- run_FI_mediation_by_depression(cov_base_behavior, nhanes_design)
result_dep_full <- run_FI_mediation_by_depression(covariates_full, nhanes_design)


# 7.2 Build Summary Table for All Hypotheses -----

# Summary of all hypotheses in long form
combined_mediation_long <- data.frame(
  Hypothesis = rep(c(
    "H3a: FI mediates SDOH → Mortality",
    "H3b: AHEI mediates SDOH → Mortality",
    "H4: FI + AHEI jointly mediate SDOH → Mortality",
    "H5: AHEI mediates FI → Mortality",
    "H6: Depression mediates FI → Mortality"
  ), each = 3),
  
  Covariate_Set = rep(c("Base", "Base + Behavior", "Full"), times = 5)
)

# Helper function to compute CI
get_ci <- function(beta, se) {
  lower <- exp(beta - 1.96 * se)
  upper <- exp(beta + 1.96 * se)
  return(c(round(lower, 2), round(upper, 2)))
}

# HR Total
betas_total <- c(
  result_sdoh_base$beta_total,
  result_sdoh_behavior$beta_total,
  result_sdoh_full$beta_total,
  
  result_sdoh_base$beta_total,
  result_sdoh_behavior$beta_total,
  result_sdoh_full$beta_total,
  
  log(result_joint_base$HR_total),
  log(result_joint_behavior$HR_total),
  log(result_joint_full$HR_total),
  
  log(result_FI_base$HR_total),
  log(result_FI_behavior$HR_total),
  log(result_FI_full$HR_total),
  
  log(result_dep_base$HR_total),
  log(result_dep_behavior$HR_total),
  log(result_dep_full$HR_total)
)

ses_total <- c(
  result_sdoh_base$se_total,
  result_sdoh_behavior$se_total,
  result_sdoh_full$se_total,
  
  result_sdoh_base$se_total,
  result_sdoh_behavior$se_total,
  result_sdoh_full$se_total,
  
  result_joint_base$se_total,
  result_joint_behavior$se_total,
  result_joint_full$se_total,
  
  result_FI_base$se_total,
  result_FI_behavior$se_total,
  result_FI_full$se_total,
  
  result_dep_base$se_total,
  result_dep_behavior$se_total,
  result_dep_full$se_total
)

cis_total <- t(mapply(get_ci, betas_total, ses_total))
combined_mediation_long$HR_Total <- round(exp(betas_total), 2)
combined_mediation_long$HR_Total_LCL <- cis_total[,1]
combined_mediation_long$HR_Total_UCL <- cis_total[,2]

# HR Adjusted
betas_adj <- c(
  result_sdoh_base$beta_FI,
  result_sdoh_behavior$beta_FI,
  result_sdoh_full$beta_FI,
  
  result_sdoh_base$beta_AHEI,
  result_sdoh_behavior$beta_AHEI,
  result_sdoh_full$beta_AHEI,
  
  log(result_joint_base$HR_joint),
  log(result_joint_behavior$HR_joint),
  log(result_joint_full$HR_joint),
  
  log(result_FI_base$HR_mediator),
  log(result_FI_behavior$HR_mediator),
  log(result_FI_full$HR_mediator),
  
  log(result_dep_base$HR_mediator),
  log(result_dep_behavior$HR_mediator),
  log(result_dep_full$HR_mediator)
)

ses_adj <- c(
  result_sdoh_base$se_FI,
  result_sdoh_behavior$se_FI,
  result_sdoh_full$se_FI,
  
  result_sdoh_base$se_AHEI,
  result_sdoh_behavior$se_AHEI,
  result_sdoh_full$se_AHEI,
  
  result_joint_base$se_joint,
  result_joint_behavior$se_joint,
  result_joint_full$se_joint,
  
  result_FI_base$se_mediator,
  result_FI_behavior$se_mediator,
  result_FI_full$se_mediator,
  
  result_dep_base$se_mediator,
  result_dep_behavior$se_mediator,
  result_dep_full$se_mediator
)

cis_adj <- t(mapply(get_ci, betas_adj, ses_adj))
combined_mediation_long$HR_Adjusted <- round(exp(betas_adj), 2)
combined_mediation_long$HR_Adj_LCL <- cis_adj[,1]
combined_mediation_long$HR_Adj_UCL <- cis_adj[,2]

# Add Prop Mediated and P-value (same as before)
combined_mediation_long$Prop_Mediated <- round(c(
  result_sdoh_base$prop_FI,
  result_sdoh_behavior$prop_FI,
  result_sdoh_full$prop_FI,
  
  result_sdoh_base$prop_AHEI,
  result_sdoh_behavior$prop_AHEI,
  result_sdoh_full$prop_AHEI,
  
  result_joint_base$prop_joint,
  result_joint_behavior$prop_joint,
  result_joint_full$prop_joint,
  
  result_FI_base$prop_mediated,
  result_FI_behavior$prop_mediated,
  result_FI_full$prop_mediated,
  
  result_dep_base$prop_mediated,
  result_dep_behavior$prop_mediated,
  result_dep_full$prop_mediated
) * 100, 2)

combined_mediation_long$P_Value <- signif(c(
  result_sdoh_base$pval_FI,
  result_sdoh_behavior$pval_FI,
  result_sdoh_full$pval_FI,
  
  result_sdoh_base$pval_AHEI,
  result_sdoh_behavior$pval_AHEI,
  result_sdoh_full$pval_AHEI,
  
  result_joint_base$pval_joint,
  result_joint_behavior$pval_joint,
  result_joint_full$pval_joint,
  
  result_FI_base$pval,
  result_FI_behavior$pval,
  result_FI_full$pval,
  
  result_dep_base$pval,
  result_dep_behavior$pval,
  result_dep_full$pval
), 2)



print(combined_mediation_long)

##### Create CI helper function ----
get_hr_ci <- function(beta, se) {
  HR <- exp(beta)
  LCL <- exp(beta - 1.96 * se)
  UCL <- exp(beta + 1.96 * se)
  return(round(c(HR, LCL, UCL), 2))
}

##### Extract log(HR) and SE for FI and AHEI from model ----
beta_FI <- coef(fs_sdoh_diet_mt)["FI"]
se_FI <- sqrt(vcov(fs_sdoh_diet_mt)["FI", "FI"])
ci_FI <- get_hr_ci(beta_FI, se_FI)

beta_AHEI <- coef(fs_sdoh_diet_mt)["ahei_total"]
se_AHEI <- sqrt(vcov(fs_sdoh_diet_mt)["ahei_total", "ahei_total"])
ci_AHEI <- get_hr_ci(beta_AHEI, se_AHEI)

##### Create new rows for Hypothesis 1 and 2 with CI columns ----
hypo1_row <- data.frame(
  Hypothesis = "H1: FI → Mortality (+ SDOH, AHEI, covariates)",
  Covariate_Set = "Full",
  HR_Total = ci_FI[1],
  HR_Total_LCL = ci_FI[2],
  HR_Total_UCL = ci_FI[2],
  HR_Adjusted = NA,
  HR_Adj_LCL = NA,
  HR_Adj_UCL = NA,
  Prop_Mediated = NA,
  P_Value = signif(summary(fs_sdoh_diet_mt)$coefficients["FI", "Pr(>|z|)"], 2)
)

hypo2_row <- data.frame(
  Hypothesis = "H2: AHEI → Mortality (+ SDOH, FI, covariates)",
  Covariate_Set = "Full",
  HR_Total = ci_AHEI[1],
  HR_Total_LCL = ci_AHEI[2],
  HR_Total_UCL = ci_AHEI[2],
  HR_Adjusted = NA,
  HR_Adj_LCL = NA,
  HR_Adj_UCL = NA,
  Prop_Mediated = NA,
  P_Value = signif(summary(fs_sdoh_diet_mt)$coefficients["ahei_total", "Pr(>|z|)"], 2)
)

##### Combine with updated mediation table ----
combined_mediation_long <- rbind(
  hypo1_row,
  hypo2_row,
  combined_mediation_long
)

print(combined_mediation_long)

##### Create formatted HR_Total column with CI ----
combined_mediation_long <- combined_mediation_long %>%
  mutate(
    HR_Total = ifelse(
      !is.na(HR_Total),
      sprintf("%.2f (%.2f, %.2f)", HR_Total, HR_Total_LCL, HR_Total_UCL),
      NA
    ),
    HR_Adjusted = ifelse(
      !is.na(HR_Adjusted),
      sprintf("%.2f (%.2f, %.2f)", HR_Adjusted, HR_Adj_LCL, HR_Adj_UCL),
      NA
    )
  )


combined_mediation_long <- combined_mediation_long %>%
  select(-HR_Total_LCL, -HR_Total_UCL, -HR_Adj_LCL, -HR_Adj_UCL)

print(combined_mediation_long)

# export to doc 

# Define output path
output_path_csv <- "/Users/dengshuyue/Desktop/SDOH/analysis/output/mediation_summary_table.csv"
output_path_docx <- "/Users/dengshuyue/Desktop/SDOH/analysis/output/mediation_summary_table.docx"

# Save CSV version
combined_export <- combined_mediation_long %>%
  mutate(Hypothesis = str_replace_all(Hypothesis, "→", "->"))

write_csv(combined_export, output_path_csv)

# Create flextable
mediation_flex <- flextable(combined_mediation_long)

# Optional: autofit columns for better display
# mediation_flex <- autofit(mediation_flex)

# Create Word document with the table
doc <- read_docx() %>%
  body_add_par("Table. Mediation Analysis Summary Across Hypotheses", style = "heading 1") %>%
  body_add_flextable(mediation_flex)

# Save Word document
print(doc, target = output_path_docx)


# 8 FINAL CHECK  -----

verify_mediation <- function(model_total, model_adjusted, exposure_var = "sdoh_score") {
  # 1. Sample size consistency
  n_total <- nrow(model_total$y)
  n_adjusted <- nrow(model_adjusted$y)
  same_sample <- n_total == n_adjusted
  
  # 2. Coefficients and SEs
  beta_total <- coef(model_total)[exposure_var]
  se_total <- sqrt(vcov(model_total)[exposure_var, exposure_var])
  
  beta_adj <- coef(model_adjusted)[exposure_var]
  se_adj <- sqrt(vcov(model_adjusted)[exposure_var, exposure_var])
  
  # 3. Derived values
  HR_total <- exp(beta_total)
  HR_adjusted <- exp(beta_adj)
  prop_mediated <- (beta_total - beta_adj) / beta_total
  
  # 4. Delta method for indirect effect
  z <- (beta_total - beta_adj) / sqrt(se_total^2 + se_adj^2)
  p_val <- 2 * (1 - pnorm(abs(z)))
  
  # 5. Output
  cat("=====================================\n")
  cat("Verification for:", exposure_var, "\n")
  cat("Sample size equal?:", same_sample, "\n")
  cat("N Total:", n_total, " | N Adjusted:", n_adjusted, "\n")
  cat("HR Total:", round(HR_total, 3), "\n")
  cat("HR Adjusted:", round(HR_adjusted, 3), "\n")
  cat("Proportion Mediated (%):", round(prop_mediated * 100, 2), "\n")
  cat("P-Value (delta method):", signif(p_val, 3), "\n")
  cat("=====================================\n\n")
}

verify_mediation(sdoh_mort, sdoh_FI_mort, exposure_var = "sdoh_score")

verify_mediation(sdoh_mort, sdoh_AHEI_mort, exposure_var = "sdoh_score")

verify_mediation(FI_mort, FI_AHEI_mort, exposure_var = "FI")

verify_mediation(FI_sdoh_mort, FI_depression_mort, exposure_var = "FI")





