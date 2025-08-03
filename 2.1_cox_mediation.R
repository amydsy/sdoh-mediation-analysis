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
  "FSDHH", # Household size
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


# 2. Fit Survey-Weighted Cox Model -----

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







