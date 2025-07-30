# Health access 

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
          "tibble", "readr", "flextable", "officer", "usethis", "gert")

# Install any missing packages
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)

# Load all required packages
lapply(want, function(pkg) require(pkg, character.only = TRUE))
rm(want, need)

# 1. List your HUQ files from the correct folder -------
huq_files <- list.files(path = "data/health_access", pattern = "HUQ_.*\\.xpt$", full.names = TRUE)

# Read and combine HUQ files with relevant variables
huq_combined <- map_dfr(huq_files, read_xpt) %>%
  select(SEQN, HUQ030, HUQ040)

huq_combined <- huq_combined %>%
  mutate(
    sdoh_access = case_when(
      HUQ030 == 2 ~ 1,   # No usual place → unfavorable
      HUQ030 %in% c(1, 3) ~ 0,  # Yes or more than one place → favorable
      TRUE ~ NA_real_
    )
  )

# Save to your project data folder
fwrite(huq_combined, file.path("data", "huq_combined.csv"))

# end 



