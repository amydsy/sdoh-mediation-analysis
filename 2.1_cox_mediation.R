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

# ---- FUNCTION 1: Separate Mediation (FI and AHEI) ----
run_separate_mediation <- function(covariate_list, design_obj) {
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  model_FI <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FI +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  model_AHEI <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + ahei_total +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  # Extract
  beta_total <- coef(model_total)["sdoh_score"]
  se_total <- sqrt(vcov(model_total)["sdoh_score", "sdoh_score"])
  
  beta_FI <- coef(model_FI)["sdoh_score"]
  se_FI <- sqrt(vcov(model_FI)["sdoh_score", "sdoh_score"])
  
  beta_AHEI <- coef(model_AHEI)["sdoh_score"]
  se_AHEI <- sqrt(vcov(model_AHEI)["sdoh_score", "sdoh_score"])
  
  # Proportion mediated
  prop_FI <- (beta_total - beta_FI) / beta_total
  prop_AHEI <- (beta_total - beta_AHEI) / beta_total
  
  # Wald tests
  z_FI <- (beta_total - beta_FI) / sqrt(se_total^2 + se_FI^2)
  z_AHEI <- (beta_total - beta_AHEI) / sqrt(se_total^2 + se_AHEI^2)
  
  pval_FI <- 2 * (1 - pnorm(abs(z_FI)))
  pval_AHEI <- 2 * (1 - pnorm(abs(z_AHEI)))
  
  return(list(
    beta_total = beta_total,
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
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  model_joint <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ sdoh_score + FI + ahei_total +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  beta_total <- coef(model_total)["sdoh_score"]
  beta_joint <- coef(model_joint)["sdoh_score"]
  
  se_total <- sqrt(vcov(model_total)["sdoh_score", "sdoh_score"])
  se_joint <- sqrt(vcov(model_joint)["sdoh_score", "sdoh_score"])
  
  prop_joint <- (beta_total - beta_joint) / beta_total
  z <- (beta_total - beta_joint) / sqrt(se_total^2 + se_joint^2)
  pval <- 2 * (1 - pnorm(abs(z)))
  
  return(list(
    beta_total = beta_total,
    beta_joint = beta_joint,
    HR_total = exp(beta_total),
    HR_joint = exp(beta_joint),
    prop_joint = prop_joint,
    pval_joint = pval
  ))
}


# H5. Diet quality mediates  food insecurity → mortality,  adjusting for SDOH + covariates.-------

# ---- FUNCTION 3: AHEI mediates FI → Mortality ----
run_FI_mediation_by_ahei <- function(covariate_list, design_obj) {
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  model_adj <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + ahei_total + sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  beta_total <- coef(model_total)["FI"]
  beta_mediator <- coef(model_adj)["FI"]
  
  se_total <- sqrt(vcov(model_total)["FI", "FI"])
  se_mediator <- sqrt(vcov(model_adj)["FI", "FI"])
  
  prop <- (beta_total - beta_mediator) / beta_total
  z <- (beta_total - beta_mediator) / sqrt(se_total^2 + se_mediator^2)
  pval <- 2 * (1 - pnorm(abs(z)))
  
  return(list(
    beta_total = beta_total,
    beta_mediator = beta_mediator,
    HR_total = exp(beta_total),
    HR_mediator = exp(beta_mediator),
    prop_mediated = prop,
    pval = pval
  ))
}


# H6. Exploratory: depressive mediates food insecurity → mortality, adjusting for SDOH + covariates.-----

# ---- FUNCTION 4: Depression mediates FI → Mortality ----
run_FI_mediation_by_depression <- function(covariate_list, design_obj) {
  model_total <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  model_adj <- svycoxph(
    as.formula(paste("Surv(py, MORTSTAT) ~ FI + probable_depression + sdoh_score +", paste(covariate_list, collapse = " + "))),
    design = design_obj
  )
  
  beta_total <- coef(model_total)["FI"]
  beta_mediator <- coef(model_adj)["FI"]
  
  se_total <- sqrt(vcov(model_total)["FI", "FI"])
  se_mediator <- sqrt(vcov(model_adj)["FI", "FI"])
  
  prop <- (beta_total - beta_mediator) / beta_total
  z <- (beta_total - beta_mediator) / sqrt(se_total^2 + se_mediator^2)
  pval <- 2 * (1 - pnorm(abs(z)))
  
  return(list(
    beta_total = beta_total,
    beta_mediator = beta_mediator,
    HR_total = exp(beta_total),
    HR_mediator = exp(beta_mediator),
    prop_mediated = prop,
    pval = pval
  ))
}


# 7.1 Run Models for All Cov Sets -----
# Run all hypotheses with 3 covariate sets
result_base <- run_separate_mediation(covariates_base, nhanes_design)
result_behavior <- run_separate_mediation(cov_base_behavior, nhanes_design)
result_full <- run_separate_mediation(covariates_full, nhanes_design)

result_joint_base <- run_joint_mediation(covariates_base, nhanes_design)
result_joint_behavior <- run_joint_mediation(cov_base_behavior, nhanes_design)
result_joint_full <- run_joint_mediation(covariates_full, nhanes_design)

result_FI_base <- run_FI_mediation_by_ahei(covariates_base, nhanes_design)
result_FI_behavior <- run_FI_mediation_by_ahei(cov_base_behavior, nhanes_design)
result_FI_full <- run_FI_mediation_by_ahei(covariates_full, nhanes_design)

result_dep_base <- run_FI_mediation_by_depression(covariates_base, nhanes_design)
result_dep_behavior <- run_FI_mediation_by_depression(cov_base_behavior, nhanes_design)
result_dep_full <- run_FI_mediation_by_depression(covariates_full, nhanes_design)

# 7.2 Build Summary Table
# Summary of all hypotheses in long form
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
    
    exp(result_base$beta_total),       # H3b (same base)
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
    
    exp(result_base$beta_AHEI),        # H3b
    exp(result_behavior$beta_AHEI),
    exp(result_full$beta_AHEI),
    
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
    result_base$prop_FI,
    result_behavior$prop_FI,
    result_full$prop_FI,
    
    result_base$prop_AHEI,
    result_behavior$prop_AHEI,
    result_full$prop_AHEI,
    
    result_joint_base$prop_joint,
    result_joint_behavior$prop_joint,
    result_joint_full$prop_joint,
    
    result_FI_base$prop_mediated,
    result_FI_behavior$prop_mediated,
    result_FI_full$prop_mediated,
    
    result_dep_base$prop_mediated,
    result_dep_behavior$prop_mediated,
    result_dep_full$prop_mediated
  ) * 100, 2),
  
  P_Value = signif(c(
    result_base$pval_FI,
    result_behavior$pval_FI,
    result_full$pval_FI,
    
    result_base$pval_AHEI,
    result_behavior$pval_AHEI,
    result_full$pval_AHEI,
    
    result_joint_base$pval_joint,
    result_joint_behavior$pval_joint,
    result_joint_full$pval_joint,
    
    result_FI_base$pval,
    result_FI_behavior$pval,
    result_FI_full$pval,
    
    result_dep_base$pval,
    result_dep_behavior$pval,
    result_dep_full$pval
  ), 3)
)

# View results
print(combined_mediation_long)


# export to doc 

# Define output path
output_path_csv <- "/Users/dengshuyue/Desktop/SDOH/analysis/output/mediation_summary_table.csv"
output_path_docx <- "/Users/dengshuyue/Desktop/SDOH/analysis/output/mediation_summary_table.docx"

# Save CSV version
write_csv(combined_mediation_long, output_path_csv)

# Create flextable
mediation_flex <- flextable(combined_mediation_long)

# Optional: autofit columns for better display
mediation_flex <- autofit(mediation_flex)

# Create Word document with the table
doc <- read_docx() %>%
  body_add_par("Table. Mediation Analysis Summary Across Hypotheses", style = "heading 1") %>%
  body_add_flextable(mediation_flex)

# Save Word document
print(doc, target = output_path_docx)


