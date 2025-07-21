
# Setup: Packages and Directories-------


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
want <- c("dplyr", "survey", "foreign", "Hmisc", "data.table")

# Install any missing packages
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)

# Load all required packages
lapply(want, function(pkg) require(pkg, character.only = TRUE))
rm(want, need)


# Generate Demographic Summary Table Using Survey Design-------
# Load NHANES dataset
df <- fread(file.path(dir$data, "SODH_diet_mort.csv"))

# Remove rows with missing survey design elements
df <- df %>% filter(!is.na(wt10) & !is.na(sdmvstra) & !is.na(sdmvpsu))

# Set up the survey design
nhanes_design <- svydesign(
  id = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wt10,
  data = df,
  nest = TRUE
)

# === Categorical variables ===
cat_vars <- c("SEX", "RACE", "EDU", "pir", "SNAP", "SMK", "ALCG2", "bmic")

cat_results <- lapply(cat_vars, function(v) {
  tab <- svytable(as.formula(paste0("~", v)), design = nhanes_design)
  pct <- prop.table(tab) * 100
  means <- svymean(as.formula(paste0("~", v)), design = nhanes_design, na.rm = TRUE)
  se <- SE(means) * 100
  
  data.frame(
    Variable = v,
    Category = names(tab),
    Count = as.vector(tab),
    Mean_or_Percent = round(as.vector(pct), 1),
    SE = round(as.vector(se)[names(tab)], 2),
    # Mean = NA_real_,  # for compatibility
    Type = "Categorical"
  )
}) %>% bind_rows()

# === Binary variables ===
binary_vars <- c("DIABETES", "CVD", "dm_rx", "chol_rx", "angina", "cancer", "lung_disease", "MORTSTAT")

binary_results <- lapply(binary_vars, function(v) {
  # Unweighted count for category 1 (TRUE or 1)
  count <- sum(df[[v]] == 1, na.rm = TRUE)
  
  # Survey-weighted mean and SE
  mean_obj <- svymean(as.formula(paste0("~", v)), nhanes_design, na.rm = TRUE)
  
  data.frame(
    Variable = v,
    Category = "1",
    Count = count,
    Mean_or_Percent = round(100 * coef(mean_obj)[[1]], 1),
    SE = round(100 * SE(mean_obj)[[1]], 2),
    Type = "Binary"
  )
}) %>% bind_rows()


# === Continuous variables ===
cont_vars <- c("RIDAGEYR", "met_hr", "bmi", "hba1c", "sbp", "dbp", "hdl", "ldl", "tg", "HEI2015_TOTAL_SCORE")

cont_results <- lapply(cont_vars, function(v) {
  # Unweighted non-missing count
  count <- sum(!is.na(df[[v]]))
  
  # Survey-weighted mean and SE
  stat <- svymean(as.formula(paste0("~", v)), nhanes_design, na.rm = TRUE)
  
  data.frame(
    Variable = v,
    Count = count,
    Mean_or_Percent = round(coef(stat)[[1]], 2),
    SE = round(SE(stat)[[1]], 2),
    Type = "Continuous"
  )
}) %>% bind_rows()


# === Combine all results ===

demo_summary <- bind_rows(cat_results, binary_results, cont_results)

# === Save and view ===

head(demo_summary, 40)

write_csv(demo_summary, "/Users/dengshuyue/Desktop/SDOH/analysis/output/demo_summary_r.csv")

