# AHEI estimate 

# Set working directory to analysis folder
setwd("/Users/dengshuyue/Desktop/SDOH/analysis")

# Define directory structure
dir <- list()
dir$root    <- getwd()
dir$data    <- file.path(dir$root, "data/nhanes_deit")
dir$output  <- file.path(dir$root, "output")
dir$code    <- file.path(dir$root, "code")

# 1 Load Required Packages------
# List of required packages
want <- c("dplyr", "haven", "foreign", "survey", "purrr", "ggplot2")

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

fped_all$cycle

fped_all$FOODCODE
glimpse(fped_all$cycle)

# Optional: Save vector of included FPED cycles
fped_cycles <- names(fped_files)

# 1.4 Merge FPED with DR1IFF (food-level data) ------

# Create a named vector to map letter codes to year labels
cycle_map <- c(
  "D" = "2005-2006",
  "E" = "2007-2008",
  "F" = "2009-2010",
  "G" = "2011-2012",
  "H" = "2013-2014",
  "I" = "2015-2016",
  "J" = "2017-2018"
)

# Apply mapping to dr1iff_all
dr1iff_all <- dr1iff_all %>%
  mutate(cycle = cycle_map[cycle])


dr1iff_fped <- dr1iff_all %>%
  filter(cycle %in% fped_cycles) %>%
  left_join(fped_all, by = c("DR1IFDCD" = "FOODCODE", "cycle"))

dr1iff_all$DR1IFDCD
glimpse(dr1iff_all$cycle)

summary(dr1iff_fped$SEQN)


## gut check if cycle is correctly included ----

dr1iff_fped %>% 
  filter(SEQN == 93702) %>% 
  select(SEQN, cycle)

dr1iff_fped %>% 
  filter(SEQN == 21009) %>% 
  select(SEQN, cycle)


# 🔍 Check missing values in FPED variables
missing_summary <- dr1iff_fped %>%
  summarise(
    pct_missing_f_juice = mean(is.na(F_JUICE)) * 100,
    pct_missing_a_drinks = mean(is.na(A_DRINKS)) * 100
  )

print(missing_summary)

# Show unmatched food codes
unmatched_codes <- dr1iff_fped %>%
  filter(is.na(F_JUICE) | is.na(A_DRINKS)) %>%
  distinct(DR1IFDCD, cycle)

print(unmatched_codes)


# 1.5 OR use existing FPED file hei9918 ------

# hei_9918 <- haven::read_sas(file.path(dir$data, "hei9918.sas7bdat"))

# glimpse(hei_9918)

# rm(hei_9918)

# Note: hei_9918 includes most FPED-derived food group variables needed for AHEI.
# Missing components: red/processed meat, EPA+DHA, PUFA (% energy), alcohol — need to pull from other files.

# 2. calculate AHEI ------
# --- AHEI Scoring: Adults Only (Age ≥ 20) --- #

# Load demographics and filter to adults
demo_adult <- fread(file.path(dir$data, "SODH_diet_mort_depr.csv")) %>%
  select(SEQN, RIDAGEYR) %>%
  filter(RIDAGEYR >= 20)

# Prepare base dataset with valid adults only
ahei <- dr1iff_fped %>%
  semi_join(demo_adult, by = "SEQN")  #%>%         # restrict to adults

length(unique(ahei$SEQN))




# 2.1 VEGETABLES ---------------------------------------------------------
ahei_veg <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    kcal_total = sum(DR1IKCAL, na.rm = TRUE),
    veg_total_cup = sum(V_DRKGR + V_REDOR_TOMATO + V_REDOR_OTHER + V_OTHER + V_LEGUMES + PF_LEGUMES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    veg_per_1000kcal = veg_total_cup / (kcal_total / 1000),
    veg_per_1000kcal = ifelse(veg_per_1000kcal > 10, NA, veg_per_1000kcal),
    ahei_veg = pmin(veg_per_1000kcal / 5, 1) * 10
  )

# 2.2 WHOLE FRUIT --------------------------------------------------------
ahei_fruit <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    kcal_total = sum(DR1IKCAL, na.rm = TRUE),
    fruit_total_cup = sum(F_TOTAL - F_JUICE, na.rm = TRUE),  # Exclude juice
    .groups = "drop"
  ) %>%
  mutate(
    fruit_per_1000kcal = fruit_total_cup / (kcal_total / 1000),
    fruit_per_1000kcal = ifelse(fruit_per_1000kcal > 10, NA, fruit_per_1000kcal),
    ahei_fruit = pmin(fruit_per_1000kcal / 2, 1) * 10
  )

# 2.3 WHOLE GRAINS -------------------------------------------------------
ahei_grain <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    kcal_total = sum(DR1IKCAL, na.rm = TRUE),
    wholegrains_oz_eq = sum(G_WHOLE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    wholegrains_g = wholegrains_oz_eq * 28.35,  # Convert oz-eq to grams
    wholegrains_g = ifelse(wholegrains_g > 150, NA, wholegrains_g),  # Optional trim
    ahei_wholegrains = pmin(wholegrains_g / 75, 1) * 10
  )

# 2.4 SSBs + FRUIT JUICE -------------------------------------------------
ahei_ssb <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    ssb_servings = sum(A_DRINKS, na.rm = TRUE),
    fruit_juice_cups = sum(F_JUICE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ssb_fj_servings = ssb_servings + fruit_juice_cups,
    ahei_ssb = case_when(
      ssb_fj_servings >= 1 ~ 0,
      ssb_fj_servings < 1  ~ (1 - ssb_fj_servings) * 10
    )
  )


# 2.5 NUTS & LEGUMES ---------------------------------------------------
ahei_nutslegumes <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    nuts_legumes_oz = sum(PF_NUTSDS + PF_LEGUMES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ahei_nutslegumes = pmin(nuts_legumes_oz / 1, 1) * 10
  )


# 2.6 Red and processed meat ---------------------------------------------------

ahei_meat <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    red_proc_oz = sum(PF_MEAT + PF_CUREDMEAT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ahei_redprocmeat = pmin(red_proc_oz / 1.5, 1),
    ahei_redprocmeat = (1 - ahei_redprocmeat) * 10
  )

# 2.7 Trans fat (usually omitted or proxied) -----------------------------------

# 2.8 LONG-CHAIN OMEGA-3 (EPA+DHA) ---------------------------------------------

ahei_longn3 <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    long_chain_n3 = sum(DR1IP204 + DR1IP225 + DR1IP205, na.rm = TRUE),  # grams/day
    .groups = "drop"
  ) %>%
  mutate(
    long_chain_n3 = ifelse(long_chain_n3 > 5, NA, long_chain_n3),  # trim extreme
    ahei_longn3 = pmin(long_chain_n3 / 0.25, 1) * 10
  )

# --- 2.9 PUFAs (as % of energy) ---
ahei_pufa <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    kcal_total = sum(DR1IKCAL, na.rm = TRUE),
    pufa_g = sum(DR1IPFAT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pufa_energy_pct = (pufa_g * 9) / kcal_total * 100,
    ahei_pufa = case_when(
      pufa_energy_pct >= 10 ~ 10,
      pufa_energy_pct <= 2 ~ 0,
      TRUE ~ (pufa_energy_pct - 2) / (10 - 2) * 10
    )
  )

# --- 2.10 SODIUM (reverse scored) ---
ahei_sodium <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    sodium_mg = sum(DR1ISODI, na.rm = TRUE),
    kcal_total = sum(DR1IKCAL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sodium_per_1000kcal = sodium_mg / (kcal_total / 1000),
    ahei_sodium = case_when(
      sodium_per_1000kcal <= 1000 ~ 10,
      sodium_per_1000kcal >= 2300 ~ 0,
      TRUE ~ (2300 - sodium_per_1000kcal) / (2300 - 1000) * 10
    )
  )

# --- 2.11 ALCOHOL ---
ahei_alcohol <- ahei %>%
  group_by(SEQN) %>%
  summarise(
    alcohol_g = sum(DR1IALCO, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ahei_alcohol = case_when(
      alcohol_g == 0 ~ 2.5,
      alcohol_g > 0 & alcohol_g <= 13 ~ 10,
      alcohol_g > 13 & alcohol_g <= 26 ~ (26 - alcohol_g) / (26 - 13) * 10,
      alcohol_g > 26 ~ 0,
      TRUE ~ NA_real_
    )
  )



# 2.12 COMBINE ALL COMPONENTS --------------------------------------------
ahei_combined <- list(
  ahei_veg, ahei_fruit, ahei_grain, ahei_ssb, ahei_nutslegumes,
  ahei_meat, ahei_longn3, ahei_pufa, ahei_sodium, ahei_alcohol
) %>%
  reduce(full_join, by = "SEQN") %>%
  mutate(
    ahei_total = rowSums(select(., starts_with("ahei_")), na.rm = TRUE)
  )

# OPTIONAL: Summary
summary(ahei_combined$ahei_total)
summary(ahei_combined$ahei_pufa)

# Save AHEI combined scores to data folder
fwrite(ahei_combined, file = file.path(dir$data, "ahei_combined.csv"))

# Histogram of AHEI total score
ggplot(ahei_combined, aes(x = ahei_total)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black", alpha = 0.8) +
  labs(
    title = "Distribution of AHEI Total Scores",
    x = "AHEI Total Score",
    y = "Count"
  ) +
  theme_minimal()













