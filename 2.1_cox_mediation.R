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
df <- fread(file.path(dir$data, "SODH_diet_mort3.csv"))

# Check SDOH score and mortality variables
names(df)

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
  "FSDHH",          # Household size
  "SDDSRVYR"        # NHANES cycle
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
# FS → mortality, independent of SDOH, diet quality, and other covariates

fs_sdoh_diet_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ FS + sdoh_score + ahei_total +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)

# fs_mt_depr <- svycoxph(
#  as.formula(paste("Surv(py, MORTSTAT) ~ FS + sdoh_score + ahei_total + probable_depression +", 
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

fs_sdoh_diet_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ FS + sdoh_score + ahei_total +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)
summary(fs_sdoh_diet_mt)


# 3.3 H3: SDOH → FS or AHEI → Mortality (Separate Mediation) -----
# hypo 3. Food insecurity, diet quality separately mediate SDOH  → mortality.

###### model separately with only full cov list -----
# Step 1: Total effect model (SDOH → Mortality)
sdoh_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)

summary(sdoh_mt)

# Step 2: SDOH + FS → Mortality
sdoh_fs_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FS +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)

summary(sdoh_fs_mt)

# Step 3: SDOH + AHEI → Mortality
sdoh_ahei_mt <- svycoxph(
  as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + ahei_total +", 
                   paste(covariates_full, collapse = " + "))),
  design = nhanes_design
)

summary(sdoh_ahei_mt)

# Step 4: Extract SDOH coefficients
beta_total    <- coef(sdoh_mt)["sdoh_score"]
beta_fs       <- coef(sdoh_fs_mt)["sdoh_score"]
beta_ahei     <- coef(sdoh_ahei_mt)["sdoh_score"]

# Step 5: Calculate proportion mediated
prop_fs   <- (beta_total - beta_fs) / beta_total
prop_ahei <- (beta_total - beta_ahei) / beta_total

# Step 6: Display results
cat("Proportion of SDOH effect on mortality mediated by FS:   ", round(prop_fs * 100, 2), "%\n")
cat("Proportion of SDOH effect on mortality mediated by AHEI: ", round(prop_ahei * 100, 2), "%\n")

# Optional: Save for summary table or export
mediation_results <- list(
  beta_total = beta_total,
  beta_fs = beta_fs,
  prop_fs = prop_fs,
  beta_ahei = beta_ahei,
  prop_ahei = prop_ahei
)


##### create a function for cov set ------
# Define function to run mediation analysis
run_mediation_analysis <- function(covariate_list, design_obj) {
  # Step 1: Total effect model
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Step 2: SDOH + FS
  model_fs <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FS +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Step 3: SDOH + AHEI
  model_ahei <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + ahei_total +", 
                     paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  # Extract coefficients
  beta_total <- coef(model_total)["sdoh_score"]
  beta_fs <- coef(model_fs)["sdoh_score"]
  beta_ahei <- coef(model_ahei)["sdoh_score"]
  
  # Proportion mediated
  prop_fs <- (beta_total - beta_fs) / beta_total
  prop_ahei <- (beta_total - beta_ahei) / beta_total
  
  # Return results
  return(list(
    beta_total = beta_total,
    beta_fs = beta_fs,
    prop_fs = prop_fs,
    beta_ahei = beta_ahei,
    prop_ahei = prop_ahei
  ))
}

# Run for each covariate set
result_base <- run_mediation_analysis(covariates_base, nhanes_design)
result_behavior <- run_mediation_analysis(cov_base_behavior, nhanes_design)
result_full <- run_mediation_analysis(covariates_full, nhanes_design)

# Display results
cat("\n--- Proportion Mediated by Covariate Set ---\n")
cat("Base Covariates:\n")
cat("FS:   ", round(result_base$prop_fs * 100, 2), "%\n")
cat("AHEI: ", round(result_base$prop_ahei * 100, 2), "%\n\n")

cat("Base + Behavior Covariates:\n")
cat("FS:   ", round(result_behavior$prop_fs * 100, 2), "%\n")
cat("AHEI: ", round(result_behavior$prop_ahei * 100, 2), "%\n\n")

cat("Full Covariates:\n")
cat("FS:   ", round(result_full$prop_fs * 100, 2), "%\n")
cat("AHEI: ", round(result_full$prop_ahei * 100, 2), "%\n")

# Create data frame with results
mediation_table <- data.frame(
  Covariate_Set = c("Base", "Base + Behavior", "Full"),
  Beta_Total    = c(result_base$beta_total,
                    result_behavior$beta_total,
                    result_full$beta_total),
  Beta_FS       = c(result_base$beta_fs,
                    result_behavior$beta_fs,
                    result_full$beta_fs),
  Prop_FS       = round(c(result_base$prop_fs,
                          result_behavior$prop_fs,
                          result_full$prop_fs) * 100, 2),
  Beta_AHEI     = c(result_base$beta_ahei,
                    result_behavior$beta_ahei,
                    result_full$beta_ahei),
  Prop_AHEI     = round(c(result_base$prop_ahei,
                          result_behavior$prop_ahei,
                          result_full$prop_ahei) * 100, 2)
)

# Print the table
print(mediation_table)

# 4. Food insecurity and diet quality together mediate SDOH  → mortality 
# (i.e., tested in the same model).

# 5. Diet quality mediates  food insecurity → mortality,  adjusting for SDOH + covariates.

# 6. Exploratory: depressive mediates food insecurity → mortality, adjusting for SDOH + covariates.









