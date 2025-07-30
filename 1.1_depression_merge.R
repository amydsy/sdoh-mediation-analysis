# merge PHQ-9 depression data 

# Set working directory to analysis folder
setwd("/Users/dengshuyue/Desktop/SDOH/analysis")

# Define directory structure
dir <- list()
dir$root    <- getwd()
dir$data    <- file.path(dir$root, "data")
dir$output  <- file.path(dir$root, "output")
dir$code    <- file.path(dir$root, "code")

# 1 Load Required Packages------
# List of required packages
want <- c("dplyr", "haven", "foreign", "survey")

# Install any missing packages
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)

# Load all required packages
lapply(want, function(pkg) require(pkg, character.only = TRUE))
rm(want, need)


# 2.1 Define PHQ-9 variable names
phq_vars <- c(
  "DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050",
  "DPQ060", "DPQ070", "DPQ080", "DPQ090"
)

# Define cycles and file suffixes
cycles <- c("DPQ_E", "DPQ_F", "DPQ_G", "DPQ_H", "DPQ_I", "DPQ_J")  # 2005–2018
years <- c("0506", "0708", "0910", "1112", "1314", "1516", "1718")

# Function to load and clean PHQ-9 data from one file
read_phq_file <- function(filename, cycle) {
  path <- file.path(dir$data, paste0(filename, ".XPT"))
  if (!file.exists(path)) return(NULL)
  
  df <- read_xpt(path) %>%
    select(SEQN, one_of(phq_vars)) %>%
    mutate(cycle = cycle)
  
  return(df)
}

# 2.2 Load and bind all PHQ-9 files ------ 
phq_list <- Map(read_phq_file, cycles, years)
phq_data <- bind_rows(phq_list)

# Assign score: valid values are 0–3. Set others (7=Refused, 9=Don't know) to NA
phq_data <- phq_data %>%
  mutate(across(all_of(phq_vars), ~ ifelse(.x %in% 0:3, .x, NA_integer_)))

# Compute total PHQ-9 score and depression indicator (e.g., score ≥ 10)
phq_data <- phq_data %>%
  rowwise() %>%
  mutate(
    phq9_total = sum(c_across(all_of(phq_vars)), na.rm = FALSE),
    phq9_complete = all(!is.na(c_across(all_of(phq_vars)))),
    probable_depression = ifelse(phq9_complete & phq9_total >= 10, 1,
                                 ifelse(phq9_complete, 0, NA_integer_))
  ) %>%
  ungroup()

# View summary
summary(phq_data$phq9_total)
table(phq_data$probable_depression, useNA = "always")

write_csv(phq_data, file.path(dir$data, "depression_combined.csv"))

# 3.0 Merge with SODH_diet_mort.csv and esimtate prevalence of depression for gut check -------

# 3.1. Read in your depression data and main dataset
phq_data <- read_csv(file.path(dir$data, "depression_combined.csv"))
main_data <- read_csv(file.path(dir$data, "SODH_diet_mort2.csv"))

# 3.2 Filter to the most recent cycle before merging ------
phq_data <- phq_data %>%
  mutate(cycle_num = as.numeric(cycle)) %>%
  arrange(SEQN, desc(cycle_num)) %>%
  group_by(SEQN) %>%
  slice(1) %>%
  ungroup()

# 3.3 Merge on SEQN (respondent ID) -----
merged_data <- main_data %>%
  left_join(phq_data %>% select(SEQN, probable_depression), by = "SEQN")

write_csv(merged_data, file.path(dir$data, "SODH_diet_mort_depr.csv"))

# 3.4 Create survey design object
nhanes_design <- merged_data %>%
  filter(!is.na(sdmvpsu) & !is.na(sdmvstra) & !is.na(wt)) %>%
  svydesign(
    ids = ~sdmvpsu,
    strata = ~sdmvstra,
    weights = ~wt,
    nest = TRUE,
    data = .
  )

# Estimate weighted prevalence
svymean(~probable_depression, nhanes_design, na.rm = TRUE)

#                      mean    SE
# probable_depression 0.083453 0.003

# Note: The weighted prevalence of probable depression (PHQ-9 ≥ 10) is ~8.3%,
# which aligns with national estimates (~7.6%–8.1%) from CDC (NCHS Data Briefs No. 172 & 303).










