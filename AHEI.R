# AHEI estimate 

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
want <- c("dplyr", "haven", "foreign", "survey", "purrr")

# Install any missing packages
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)

# Load all required packages
lapply(want, function(pkg) require(pkg, character.only = TRUE))
rm(want, need)

# 1.0 Define NHANES cycles and suffixes ------
nhanes_cycles <- list(
  "2003-2004" = "C", "2005-2006" = "D", "2007-2008" = "E", "2009-2010" = "F",
  "2011-2012" = "G", "2013-2014" = "H", "2015-2016" = "I", "2017-2018" = "J"
)


# 1.1 Function to read and label DR1TOT or DR1IFF files -------
read_nhanes_file <- function(prefix, suffix) {
  file_path <- file.path(dir$data, paste0(prefix, "_", suffix, ".XPT"))
  if (file.exists(file_path)) {
    df <- read_xpt(file_path)
    df$cycle <- suffix  # Add a label to identify cycle
    return(df)
  } else {
    warning("Missing: ", file_path)
    return(NULL)
  }
}

# 1.2 Read all cycles ----- 
dr1tot_list <- lapply(nhanes_cycles, function(suffix) read_nhanes_file("DR1TOT", suffix))
dr1iff_list <- lapply(nhanes_cycles, function(suffix) read_nhanes_file("DR1IFF", suffix))

# Combine all cycles
dr1tot_all <- bind_rows(dr1tot_list)
dr1iff_all <- bind_rows(dr1iff_list)

nrow(dr1tot_all)    # Should be ~77,000 (1 row per person)
nrow(dr1iff_all)    # Much larger (multiple rows per person)

# Merge on SEQN and cycle
diet_merged <- inner_join(dr1tot_all, dr1iff_all, by = c("SEQN", "cycle"))

# 1.3 Load FPED files (only 2005–2018 cycles) ------

# Define cycle short codes and paths
cycles <- c("0506", "0708", "0910", "1112", "1314", "1516", "1718")

# Create a named vector: names are "2005-2006", etc.; values are full file paths
fped_files <- setNames(
  file.path("/Users/dengshuyue/Desktop/SDOH/analysis/data/fped", paste0("fped_", cycles, ".sas7bdat")),
  paste0("20", substr(cycles, 1, 2), "-20", substr(cycles, 3, 4))
)

# Function to read each FPED file and tag with cycle info
read_fped <- function(cycle, path) {
  df <- read_sas(path)
  df$cycle <- cycle
  df
}

# Read and bind all FPED files
fped_all <- map2_dfr(names(fped_files), fped_files, read_fped)

fped_all$FOODCODE
glimpse(fped_all$cycle)

# Optional: Save vector of included FPED cycles
fped_cycles <- names(fped_files)

# 1.4 Merge FPED with DR1IFF (food-level data) ------

# Create a named vector to map letter codes to year labels
cycle_map <- c(
  "C" = "2005-2006",
  "D" = "2007-2008",
  "E" = "2009-2010",
  "F" = "2011-2012",
  "G" = "2013-2014",
  "H" = "2015-2016",
  "I" = "2017-2018"
)

# Apply mapping to dr1iff_all
dr1iff_all <- dr1iff_all %>%
  mutate(cycle = cycle_map[cycle])

dr1iff_fped <- dr1iff_all %>%
  filter(cycle %in% fped_cycles) %>%
  left_join(fped_all, by = c("DR1IFDCD" = "FOODCODE", "cycle"))

dr1iff_all$DR1IFDCD
glimpse(dr1iff_all$cycle)

# 1.5 OR use existing FPED file hei9918 ------

# hei_9918 <- haven::read_sas(file.path(dir$data, "hei9918.sas7bdat"))

# glimpse(hei_9918)

# rm(hei_9918)

# Note: hei_9918 includes most FPED-derived food group variables needed for AHEI.
# Missing components: red/processed meat, EPA+DHA, PUFA (% energy), alcohol — need to pull from other files.

# 2. calculate AHEI ------

glimpse(dr1iff_fped)

# Check dimensions and sample
dim(dr1iff_fped)
glimpse(dr1iff_fped)

# Check that FPED variables (e.g., g_whole1, add_sugars1) are included
names(dr1iff_fped)


# Start with a copy of your dataset
ahei <- dr1iff_fped

# Total energy intake in kcal
ahei$kcal <- ahei$kcal  # already available

#  🔥 🔥 🔥 🔥 🔥 🔥2.1. AHEI Component Calculations  -------

ahei <- dr1iff_fped %>%
  mutate(
    kcal = DR1IKCAL,  # Total energy intake
    
    # 1. Vegetables (excluding potatoes and juice)
    ahei_veg = pmin((V_DRKGR + V_REDOR_OTHER + V_OTHER + V_LEGUMES) / (5 / 1000 * kcal), 1) * 10,
    
    # 2. Fruits (excluding juice)
    ahei_fruit = pmin((F_CITMLB + F_OTHER) / (4 / 1000 * kcal), 1) * 10,
    
    # 3. Whole Grains (grams/day)
    ahei_wholegrains = pmin(G_WHOLE, 75) / 75 * 10,
    
    # 4. SSB + Fruit Juice (reverse scored)
    ahei_ssbjuice = pmax(0, pmin(1 - ((A_DRINKS + F_JUICE) / (1.5 / 1000 * kcal)), 1)) * 10,
    
    # 5. Nuts and Legumes
    ahei_nutslegumes = pmin((PF_NUTSDS + PF_LEGUMES) / (1 / 1000 * kcal), 1) * 10,
    
    # 6. Red and Processed Meat (reverse scored)
    ahei_redmeat = pmax(0, pmin(1 - ((PF_MEAT + PF_CUREDMEAT) / (1.5 / 1000 * kcal)), 1)) * 10,
    
    # 7. Trans fat — not available in NHANES, skip or set to NA
    ahei_transfat = NA_real_,
    
    # 8. Long-chain omega-3 fats (EPA + DHA)
    epa_dha = DR1IP205 + DR1IP225 + DR1IP226,
    ahei_omega3 = pmin(epa_dha, 250) / 250 * 10,
    
    # 9. PUFA (% of total energy)
    pufa = DR1IP182 + DR1IP183,
    ahei_pufa = pmin(pmax((pufa / kcal * 100 - 2) / (10 - 2), 0), 1) * 10,
    
    # 10. Sodium (reverse scored by decile)
    ahei_sodium = 10 - ntile(DR1ISODI, 10),
    
    # 11. Alcohol (gender-specific scoring)
    alcohol = DR1IALCO,
    ahei_alcohol = ifelse(
      MODCODE == 1,  # Male
      ifelse(alcohol >= 0.5 & alcohol <= 2.0, 10,
             ifelse(alcohol < 0.5, alcohol / 0.5 * 2.5,
                    pmax(0, 2.5 - (alcohol - 2) / 1.5 * 2.5))),
      ifelse(alcohol >= 0.5 & alcohol <= 1.5, 10,
             ifelse(alcohol < 0.5, alcohol / 0.5 * 2.5,
                    pmax(0, 2.5 - (alcohol - 1.5) / 1.0 * 2.5)))
    )
  ) %>%
  mutate(
    # Final AHEI score (excluding transfat)
    ahei_total = rowSums(select(., starts_with("ahei_"))[ , -which(names(select(., starts_with("ahei_"))) == "ahei_transfat")], na.rm = TRUE)
  )

ahei %>% summarise(
  mean = mean(ahei_total, na.rm = TRUE),
  sd = sd(ahei_total, na.rm = TRUE)
)

summary(ahei$ahei_total)

table(cut(ahei$ahei_total, breaks = c(0, 20, 40, 60, 80, 100, 120)))


# Get column names that start with 'ahei_' but exclude 'ahei_total'
ahei_cols <- setdiff(grep("^ahei_", names(ahei), value = TRUE), "ahei_total")

# Count how many components are non-missing for each row
ahei$non_missing_components <- rowSums(!is.na(ahei[, ahei_cols]))

# Show summary
summary(ahei$non_missing_components)















