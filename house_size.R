# household size

# Household Size Extraction Script -------

# 1 Setup: Packages and Directories -------

# Set working directory to analysis folder
setwd("/Users/dengshuyue/Desktop/SDOH/analysis")

# Define directory structure
dir <- list()
dir$root    <- getwd()
dir$data    <- file.path(dir$root, "data")
dir$output  <- file.path(dir$root, "output")
dir$code    <- file.path(dir$root, "code")

# Load Required Packages ------
want <- c("dplyr", "haven", "data.table", "purrr")  # use haven instead of foreign for modern .xpt
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)
lapply(want, function(pkg) require(pkg, character.only = TRUE))
rm(want, need)

# 2. List your demographic files (household size lives here) -------
demo_files <- list.files(path = "data/household_size", pattern = "\\.xpt$", full.names = TRUE)

# 3. Read and combine only needed variables: SEQN and DMDHHSIZ -------
hhsize_combined <- map_dfr(demo_files, ~ read_xpt(.x) %>%
                             select(SEQN, DMDHHSIZ))

# 4. Rename and save -------
hhsize_combined <- hhsize_combined %>%
  rename(household_size = DMDHHSIZ)

# Optional: quick check
summary(hhsize_combined$household_size)

# 5. Save to your project data folder -------
fwrite(hhsize_combined, file.path("data", "household_size_combined.csv"))





