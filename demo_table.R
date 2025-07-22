
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
want <- c("dplyr", "survey", "foreign", "Hmisc", "data.table", "tidyr", 
          "tibble", "readr")

# Install any missing packages
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)

# Load all required packages
lapply(want, function(pkg) require(pkg, character.only = TRUE))
rm(want, need)


# Generate Demographic Summary Table Using Survey Design-------
# Load NHANES dataset
df <- fread(file.path(dir$data, "SODH_diet_mort.csv"))

df <- df %>% filter(!is.na(wt10) & !is.na(sdmvstra) & !is.na(sdmvpsu))

nhanes_design <- svydesign(
  id = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wt10,
  data = df,
  nest = TRUE
)

# === Variable name mapping (original variable names only) ===
variable_labels <- c(
  SEX = "sex", RACE = "Race", EDU = "edu", pir = "Family income to poverty ratio",
  SNAP = "SNAP", SMK = "smk", ALCG2 = "Drinking", bmic = "BMI",
  DIABE = "Diabetes", CVD = "CVD", dm_rx = "DM_rx", chol_rx = "Cholestory",
  angina = "ang", cancer = "Cancer", lung_disease = "lung-disease", MORTSTAT = "Death",
  RIDAGEYR = "Age, years", met_hr = "Physical activity", bmi = "BMI", hba1c = "hba1c",
  sbp = "sbp", dbp = "dbp", hdl = "hdl", ldl = "ldl", tg = "tg",
  HEI2015_TOTAL_SCORE = "hei2015_"
)

# Helper function to map variable names once only
map_variable_labels_once <- function(var_vector) {
  sapply(var_vector, function(v) {
    if (v %in% names(variable_labels)) {
      variable_labels[[v]]
    } else {
      v
    }
  }, USE.NAMES = FALSE)
}

# === Categorical variables ===
cat_vars <- c("SEX", "RACE", "EDU", "pir", "SNAP", "SMK", "ALCG2", "bmic")

cat_results <- lapply(cat_vars, function(v) {
  unweighted_tab <- table(df[[v]])
  tab <- svytable(as.formula(paste0("~", v)), design = nhanes_design)
  pct <- prop.table(tab) * 100
  means <- svymean(as.formula(paste0("~", v)), design = nhanes_design, na.rm = TRUE)
  se <- SE(means) * 100
  
  data.frame(
    Variable = v,
    Category = names(tab),
    Count = as.vector(unweighted_tab[names(tab)]),
    Mean_or_Percent = round(as.vector(pct), 1),
    SE = round(as.vector(se)[names(tab)], 2),
    Type = "Categorical"
  )
}) %>% bind_rows()

# === Binary variables ===
binary_vars <- c("DIABE", "CVD", "dm_rx", "chol_rx", "angina", "cancer", "lung_disease", "MORTSTAT")

binary_results <- lapply(binary_vars, function(v) {
  count <- sum(nhanes_design$variables[[v]] == 1, na.rm = TRUE)
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
  count <- sum(!is.na(df[[v]]))
  stat <- svymean(as.formula(paste0("~", v)), nhanes_design, na.rm = TRUE)
  
  data.frame(
    Variable = v,
    Count = count,
    Mean_or_Percent = round(coef(stat)[[1]], 2),
    SE = round(SE(stat)[[1]], 2),
    Type = "Continuous"
  )
}) %>% bind_rows()

# === Rename variables once ===
cat_results$Variable <- map_variable_labels_once(cat_results$Variable)
binary_results$Variable <- map_variable_labels_once(binary_results$Variable)
cont_results$Variable <- map_variable_labels_once(cont_results$Variable)

# === Category label mapping ===
category_labels <- list(
  sex = c("1" = "Male", "2" = "Female"),
  Race = c("1" = "Non-Hispanic White", "2" = "Non-Hispanic Black", "3" = "Hispanic", "4" = "Other"),
  edu = c("1" = "Less than high school", "2" = "High school or equivalent", "3" = "Some college", "4" = "College or above"),
  `Family income to poverty ratio` = c("1" = "< 1.3", "2" = "1.3~2.99", "3" = ">=3"),
  SNAP = c("1" = "Not participant", "2" = "Participant", "3" = "income eligible non-participant"),
  smk = c("1" = "Nonsmokers", "2" = "Former smokers", "3" = "<15 cigarettes/day", "4" = "15-24.9 cigarettes/day", "5" = "≥ 25 cigarettes/day"),
  Drinking = c("1" = "Nondrinkers", "2" = "Moderate drinker", "3" = "Heavy drinker", "4" = "Missing"),
  BMI = c("1" = "BMI <18.5 kg/m2", "2" = "18-24.9 kg/m2", "3" = "25-29.9 kg/m2", "4" = "BMI ≥30 kg/m2")
)

cat_results <- cat_results %>%
  rowwise() %>%
  mutate(
    Category = if (Variable %in% names(category_labels)) {
      label_map <- category_labels[[Variable]]
      ifelse(Category %in% names(label_map), label_map[Category], Category)
    } else {
      Category
    }
  ) %>%
  ungroup()

# === Final merge ===
demo_summary <- bind_rows(cat_results, binary_results, cont_results)

# Define which variables should show Mean (SD) instead of Count (%)
mean_sd_vars <- c(
  "Age, years", "Physical activity", "hba1c", 
  "sbp", "dbp", "hdl", "ldl", "tg", "hei2015_"
)

# Create unified display column
demo_summary <- demo_summary %>%
  mutate(
    `Primary population` = case_when(
      Variable %in% mean_sd_vars ~ paste0(Mean_or_Percent, " (", SE, ")"),
      TRUE ~ paste0(Count, " (", Mean_or_Percent, "%)")
    )
  ) %>%
  select(Variable, Category, `Primary population`)


#### group those varibales -----

# Define groupings
variable_groups <- list(
  "Sociodemographics" = c("sex", "Race", "edu", "Family income to poverty ratio", "SNAP"),
  "Health Behaviors" = c("smk", "Drinking", "Physical activity"),
  "Clinical Characteristics" = c("BMI", "hba1c", "Diabetes", "DM_rx", "CVD", "ang", "Cancer", "Cholestory", "lung-disease", "Death"),
  "Dietary & Physiologic Measures" = c("Age, years", "sbp", "dbp", "hdl", "ldl", "tg", "hei2015_")
)

group_lookup <- enframe(variable_groups, name = "Group", value = "Variable") %>%
  unnest(cols = c(Variable))

# Attach the group label and sort:
demo_summary <- demo_summary %>%
  left_join(group_lookup, by = "Variable") %>%
  arrange(factor(Group, levels = names(variable_groups)), 
          factor(Variable, levels = unlist(variable_groups)), 
          Category)

# Display cleaner table with group headers
demo_summary <- demo_summary %>%
  mutate(
    Variable = ifelse(duplicated(Variable), "", Variable),
    Group = ifelse(duplicated(Group), "", Group)
  ) %>%
  select(Group, Variable, Category, `Primary population`)


# === View first rows ===
head(demo_summary, 50)

write_csv(demo_summary, "/Users/dengshuyue/Desktop/SDOH/analysis/output/demo_summary_r.csv")

