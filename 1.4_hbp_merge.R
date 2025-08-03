# merge hypertension related var

# 1. Setup -------
library(haven)
library(dplyr)
library(readr)

setwd("/Users/dengshuyue/Desktop/SDOH/analysis")

# Define directory structure
dir <- list()
dir$data <- file.path(getwd(), "data")
dir$data_bpq <- file.path(getwd(), "data/bpq")  # NEW path for BPQ files

# Define NHANES BPQ files and corresponding cycles
bpq_files <- paste0("BPQ_", c("C", "D", "E", "F", "G", "H", "I", "J"))  # 03–04 to 17–18
years <- c("0304", "0506", "0708", "0910", "1112", "1314", "1516", "1718")

# 2. Function to read each BPQ file ------
read_bpq_file <- function(file, cycle) {
  path <- file.path(dir$data_bpq, paste0(file, ".XPT"))  # use new path here
  if (!file.exists(path)) return(NULL)
  
  read_xpt(path) %>%
    select(SEQN, BPQ020, BPQ050A) %>%
    mutate(cycle = cycle)
}

# 3. Load and combine all files ------
bpq_list <- Map(read_bpq_file, bpq_files, years)
bpq_data <- bind_rows(bpq_list)

# 4. Deduplicate to keep latest per SEQN ------
bpq_data <- bpq_data %>%
  mutate(cycle_num = as.numeric(cycle)) %>%
  arrange(SEQN, desc(cycle_num)) %>%
  group_by(SEQN) %>%
  slice(1) %>%
  ungroup()

# 5. Save combined BPQ data ------
write_csv(bpq_data, file.path(dir$data, "bpq_combined.csv"))

# 6. Merge with main dataset ------
main_data <- read_csv(file.path(dir$data, "SODH_diet_mort2.csv"))

merged_data <- main_data %>%
  left_join(bpq_data, by = "SEQN")

# 7. Create final hypertension indicator ------
merged_data <- merged_data %>%
  mutate(
    HYPERTEN = case_when(
      !is.na(BPQ020) & BPQ020 == 1 ~ 1,                        # self-reported diagnosis
      !is.na(BPQ050A) & BPQ050A == 1 ~ 1,                      # currently taking meds
      !is.na(sbp) & sbp >= 130 ~ 1,                            # SBP ≥ 130 mm Hg
      !is.na(dbp) & dbp >= 85 ~ 1,                             # DBP ≥ 85 mm Hg
      TRUE ~ 0
    )
  )

summary(merged_data$HYPERTEN)

# gut check 
# Filter non-missing design variables
nhanes_design <- merged_data %>%
  filter(!is.na(sdmvpsu), !is.na(sdmvstra), !is.na(wt)) %>%
  svydesign(
    ids = ~sdmvpsu,
    strata = ~sdmvstra,
    weights = ~wt,
    nest = TRUE,
    data = .
  )

# Estimate weighted prevalence of hypertension
svymean(~HYPERTEN, nhanes_design, na.rm = TRUE)

# HYPERTEN 0.45356 0.0054

# 8. Save final merged dataset ------
write_csv(merged_data, file.path(dir$data, "SODH_diet_mort3.csv"))





