# AHEI 99-03


# 1. ---- Paths (edit if needed) ----


# Set working directory to analysis folder (edit if needed)
setwd("/Users/dengshuyue/Desktop/SDOH/analysis")

# Directory map
dir <- list(
  root   = getwd(),
  data   = file.path(getwd(), "data"),
  output = file.path(getwd(), "output"),
  code   = file.path(getwd(), "code"),
  fped   = file.path(getwd(), "data", "fped"),
  nhanes = file.path(getwd(), "data", "nhanes_deit")
)

# Ensure required directories exist
invisible(lapply(dir[c("data","output","code","fped","nhanes")], dir.create, showWarnings = FALSE, recursive = TRUE))

# R session defaults (stable, reproducible)
options(stringsAsFactors = FALSE, scipen = 999, repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(148)  # your lucky number :)

# Quiet package loader/installer
pkgs <- c("dplyr", "haven", "foreign", "survey", "purrr", "ggplot2", "readr", "stringr", "tidyr", "rlang")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, quiet = TRUE)

suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))

# Helpful message
cat("Project root:", dir$root, "\nData dir:", dir$data, "\nFPED/MPED dir:", dir$fped, "\n\n")


# 2 MPED MyPyramid (1999–2004) paths & reads ------------------------------------------

# Paths (use your dir list)
pyr_tot_9902_path   <- file.path(dir$fped, "pyr_tot.sas7bdat")  # MPED v1.0 totals (1999–2002)
pyr_tot_d1_0304_path<- file.path(dir$fped, "MyPyrEquivDB2", "data", "intakes", "pyr_tot_d1.sas7bdat") # MPED 2.0 (2003–2004) Day 1
pyr_tot_d2_0304_path<- file.path(dir$fped, "MyPyrEquivDB2", "data", "intakes", "pyr_tot_d2.sas7bdat") # MPED 2.0 (2003–2004) Day 2

# Safe readers
safe_read_sas <- function(path) {
  if (!file.exists(path)) {
    message("Missing file: ", path)
    return(NULL)
  }
  haven::read_sas(path)
}

# Read MPED totals
pyr_tot_9902   <- safe_read_sas(pyr_tot_9902_path)
pyr_tot_d1_0304<- safe_read_sas(pyr_tot_d1_0304_path)
pyr_tot_d2_0304<- safe_read_sas(pyr_tot_d2_0304_path)

# Tag days for 2003–2004 (only if files exist)
if (!is.null(pyr_tot_d1_0304)) pyr_tot_d1_0304 <- dplyr::mutate(pyr_tot_d1_0304, DAY = 1)
if (!is.null(pyr_tot_d2_0304)) pyr_tot_d2_0304 <- dplyr::mutate(pyr_tot_d2_0304, DAY = 2)

# Harmonize possible day variable name in 1999–2002 (often DAY or DAYNO)
if (!is.null(pyr_tot_9902)) {
  nm <- names(pyr_tot_9902)
  if ("DAYNO" %in% nm && !("DAY" %in% nm)) {
    names(pyr_tot_9902)[nm == "DAYNO"] <- "DAY"
  }
}

# Combine 2003–2004 day files if present
pyr_tot_0304 <- dplyr::bind_rows(
  if (!is.null(pyr_tot_d1_0304)) pyr_tot_d1_0304,
  if (!is.null(pyr_tot_d2_0304)) pyr_tot_d2_0304
)

# Quick status message
cat(
  "Loaded MPED files:\n",
  "- 1999–2002 totals: ", ifelse(is.null(pyr_tot_9902), "NO", "YES"), "\n",
  "- 2003–2004 Day1:   ", ifelse(is.null(pyr_tot_d1_0304), "NO", "YES"), "\n",
  "- 2003–2004 Day2:   ", ifelse(is.null(pyr_tot_d2_0304), "NO", "YES"), "\n\n",
  sep = ""
)


# 3) NHANES diet totals (1999–2004) -> merge with MPED totals (DAY 1 only) ----

# Helper to read XPT with a friendly message
safe_read_xpt <- function(path) {
  if (!file.exists(path)) { message("Missing file: ", path); return(NULL) }
  haven::read_xpt(path)
}

# Expected DR totals stored directly in dir$nhanes (no cycle subfolders)
dr_paths <- tibble::tribble(
  ~label,       ~dr1,                                        ~dr2,
  "1999-2000",  file.path(dir$nhanes, "DRXTOT.XPT"),          NA_character_,
  "2001-2002",  file.path(dir$nhanes, "DRXTOT_B.XPT"),        NA_character_,
  "2003-2004",  file.path(dir$nhanes, "DR1TOT_C.XPT"),        NA_character_  # Day 1 only
)

read_dr <- function(dr_path, day_num) {
  if (is.na(dr_path)) return(NULL)
  df <- safe_read_xpt(dr_path)
  if (is.null(df)) return(NULL)
  dplyr::mutate(df, DAY = day_num)
}

# Read DR Day 1 for all three entries
dr_list <- dr_paths %>%
  mutate(
    DR1 = purrr::map(dr1, ~ read_dr(.x, 1)),
    DR2 = purrr::map(dr2, ~ read_dr(.x, 2))  # all NA under our Day-1-only rule
  )

# Bind all Day-1 rows together
dr_all_d1 <- dplyr::bind_rows(dr_list$DR1)

# Join helper (SEQN+DAY if MPED has DAY)
join_mped <- function(dr_df, mped_df) {
  if (is.null(dr_df) || is.null(mped_df)) return(NULL)
  if ("DAY" %in% names(mped_df)) dplyr::left_join(dr_df, mped_df, by = c("SEQN","DAY"))
  else                           dplyr::left_join(dr_df, mped_df, by = "SEQN")
}

# For 1999–2002 use pyr_tot_9902; for 2003–2004 use Day-1 file only
pyr_tot_0304_d1 <- if (!is.null(pyr_tot_d1_0304)) dplyr::mutate(pyr_tot_d1_0304, DAY = 1) else NULL

# Split DR by period and join to the right MPED totals
dr_9900_0102 <- dr_all_d1 %>% dplyr::filter(file.exists(file.path(dir$nhanes, "DRXTOT.XPT")) | 
                                              file.exists(file.path(dir$nhanes, "DRXTOT_B.XPT")))
# safer explicit split:
dr_9900 <- read_dr(file.path(dir$nhanes, "DRXTOT.XPT"),   1)
dr_0102 <- read_dr(file.path(dir$nhanes, "DRXTOT_B.XPT"), 1)
dr_0304 <- read_dr(file.path(dir$nhanes, "DR1TOT_C.XPT"), 1)

dr_9902_mped <- join_mped(dplyr::bind_rows(dr_9900, dr_0102), pyr_tot_9902)
dr_0304_mped <- join_mped(dr_0304,                              pyr_tot_0304_d1)

# Combine Day-1 rows and collapse to person-level
dr_9904_mped <- dplyr::bind_rows(
  dr_9902_mped %||% dplyr::tibble(),
  dr_0304_mped %||% dplyr::tibble()
)

mped_9904_person <- dr_9904_mped %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarize(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

cat("Persons in MPED 1999–2004 (Day 1 only): ", nrow(mped_9904_person), "\n")


# Peek helpers to find columns for AHEI mapping next
peek <- function(pattern) sort(grep(pattern, names(mped_9904_person), value = TRUE, ignore.case = TRUE))
list(
  veg   = peek("veg|veget"),
  fruit = peek("fruit|juic"),
  grain = peek("grain|whole|whol"),
  nuts  = peek("nut"),
  legum = peek("legum|bean|pea"),
  meat  = peek("meat|proc|cured"),
  ssb   = peek("soda|soft|ssb|sweet|addsug|sugar")
)

# Save intermediate (optional)
readr::write_rds(mped_9904_person, file.path(dir$output, "mped_1999_2004_person_day1.rds"))

# 4) Build AHEI inputs (MPED + nutrients) -------------------------------------

# Helper: return first existing column name from candidates
pick_col <- function(df, candidates) {
  hits <- candidates[candidates %in% names(df)]
  if (length(hits)) hits[1] else NA_character_
}

# Helper: coalesce-first-existing column to numeric (or 0)
take <- function(df, candidates, default = NA_real_) {
  nm <- pick_col(df, candidates)
  if (is.na(nm)) rep(default, nrow(df)) else df[[nm]]
}

# ---- Pull nutrients from the day-1 merged DR totals (already joined in dr_9904_mped) ----
nutrients_9904 <- dr_9904_mped %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarize(
    energy_kcal = mean(dplyr::coalesce(!!sym(pick_col(., c("DRXTKCAL","DR1IKCAL"))), NA_real_), na.rm = TRUE),
    sodium_mg   = mean(dplyr::coalesce(!!sym(pick_col(., c("DRDTSODI","DR1ISODI"))), NA_real_), na.rm = TRUE),
    alcohol_g   = mean(dplyr::coalesce(!!sym(pick_col(., c("DRXTALCO","DR1IALCO"))), NA_real_), na.rm = TRUE),
    pufa_g      = mean(dplyr::coalesce(!!sym(pick_col(., c("DRXTPFAT","DR1IPFAT"))), NA_real_), na.rm = TRUE),
    # Long-chain n-3 (EPA+DHA); if unavailable they’ll be NA (component can be dropped or left NA)
    epa_mg      = mean(dplyr::coalesce(!!sym(pick_col(., c("DRXTP205","DR1IP205","EPA_MG"))), NA_real_), na.rm = TRUE),
    dha_mg      = mean(dplyr::coalesce(!!sym(pick_col(., c("DRXTP226","DR1IP226","DR1IP225","DHA_MG"))), NA_real_), na.rm = TRUE),
    .groups = "drop"
  )

# ---- MPED food groups (cups/oz-eq) pulled from mped_9904_person ---------------
# NOTE: Add or tweak candidate names if your peek() shows different names.
ahei_input_mped <- mped_9904_person %>%
  dplyr::transmute(
    SEQN,
    WTDRD1 = take(., c("WTDRD1","WTDR4YR")),  # keep Day 1 weight if present (for later survey design)
    
    # Vegetables (cups), exclude potatoes and vegetable juice if available
    veg_total_cup = take(., c("V_TOTAL","VEG_TOTAL","VEG_TOTAL_CEQ","V_TOTAL_CUP")),
    potato_cup    = take(., c("POTATOES_CEQ","POTATO_CUP","V_POTATO"), 0),
    veg_juice_cup = take(., c("VEG_JUICE_CEQ","V_JUICE","VEG_JUICE"), 0),
    veg_cup_eq    = pmax(veg_total_cup - potato_cup - veg_juice_cup, 0),
    
    # Whole fruit (cups), exclude fruit juice
    fruit_total_cup = take(., c("F_TOTAL","FRUIT_TOTAL","FRUIT_TOTAL_CEQ","F_TOTAL_CUP")),
    fruit_juice_cup = take(., c("FRUIT_JUICE_CEQ","F_JUICE","FRUIT_JUICE"), 0),
    fruit_cup_eq    = pmax(fruit_total_cup - fruit_juice_cup, 0),
    
    # Whole grains (oz-eq → grams)
    wholegr_oz_eq = take(., c("G_WHOLE","WHOLE_GRAINS_OZ_EQ","G_WHOLE_OZ")),
    wholegr_g     = wholegr_oz_eq * 28.3495,
    
    # Nuts & legumes (servings)
    nuts_oz_eq       = take(., c("PF_NUTSDS","NUTS_OZ_EQ","NUTS_OZ"), 0),
    legumes_cup_eq   = take(., c("PF_LEGUMES","LEGUMES_CUP_EQ","LEGUME_CUP"), 0),
    nuts_legumes_serv = nuts_oz_eq + (legumes_cup_eq / 0.5),
    
    # Red + processed meat (servings; 1 oz-eq = 1 serving)
    red_oz_eq        = take(., c("PF_MEAT","RED_MEAT_OZ_EQ","RED_MEAT_OZ"), 0),
    proc_oz_eq       = take(., c("PF_CUREDMEAT","PROC_MEAT_OZ_EQ","PROC_MEAT_OZ"), 0),
    redproc_serv     = red_oz_eq + proc_oz_eq,
    
    # SSB + 100% fruit juice (AHEI groups juice with SSB). MPED doesn’t give SSB servings; default NA here.
    fruit_juice_serv = fruit_juice_cup,         # cups -> treat as servings 1:1 for the SSB+juice component
    ssb_serv         = NA_real_,                # can be built later from DRXIFF beverage codes
    ssb_juice_serv   = ssb_serv + fruit_juice_serv
  ) %>%
  # Join nutrients
  dplyr::left_join(nutrients_9904, by = "SEQN")


# 🍎🍎🍎🍎🍎🍎🍎🍎🍎🍎🍎🍎
# check how we exclude potato?
# other vegetable to exclude?? 
# 🍎🍎🍎🍎🍎🍎🍎🍎🍎🍎🍎🍎

# 5) AHEI component scoring -----------------------------------------------------

# Vegetables (max at ≥5 cups/1000 kcal)
ahei_veg <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   veg_per_1000kcal = veg_cup_eq / (energy_kcal / 1000),
                   veg_per_1000kcal = dplyr::if_else(veg_per_1000kcal > 10, NA_real_, veg_per_1000kcal),
                   ahei_veg = pmin(veg_per_1000kcal / 5, 1) * 10
  )

# Whole fruit (max at ≥2 cups/1000 kcal)
ahei_fruit <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   fruit_per_1000kcal = fruit_cup_eq / (energy_kcal / 1000),
                   fruit_per_1000kcal = dplyr::if_else(fruit_per_1000kcal > 10, NA_real_, fruit_per_1000kcal),
                   ahei_fruit = pmin(fruit_per_1000kcal / 2, 1) * 10
  )

# Whole grains (max at 75 g/d)
ahei_grain <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   wholegrains_g = dplyr::if_else(wholegr_g > 150, NA_real_, wholegr_g),
                   ahei_wholegrains = pmin(wholegrains_g / 75, 1) * 10
  )

# Nuts & legumes (max at ≥1 serving/d)
ahei_nutslegumes <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   ahei_nutslegumes = pmin(nuts_legumes_serv / 1, 1) * 10
  )

# Red + processed meat (reverse; 0 at ≥1.5 servings/d)
ahei_meat <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   ahei_redprocmeat = (1 - pmin(redproc_serv / 1.5, 1)) * 10
  )

# Long-chain n-3 (EPA+DHA) (max at ≥0.25 g/d)
ahei_longn3 <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   long_chain_n3 = (dplyr::coalesce(epa_mg,0) + dplyr::coalesce(dha_mg,0)) / 1000, # mg→g
                   long_chain_n3 = dplyr::if_else(long_chain_n3 > 5, NA_real_, long_chain_n3),
                   ahei_longn3 = pmin(long_chain_n3 / 0.25, 1) * 10
  )

# PUFA (% energy) (0 at ≤2%; 10 at ≥10%)
ahei_pufa <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   pufa_energy_pct = (pufa_g * 9) / energy_kcal * 100,
                   ahei_pufa = dplyr::case_when(
                     is.na(pufa_energy_pct) ~ NA_real_,
                     pufa_energy_pct >= 10 ~ 10,
                     pufa_energy_pct <= 2  ~ 0,
                     TRUE ~ (pufa_energy_pct - 2) / (10 - 2) * 10
                   )
  )

# Sodium (reverse; 10 at ≤1000 mg/1000 kcal; 0 at ≥2300 mg/1000 kcal)
ahei_sodium <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   sodium_per_1000kcal = sodium_mg / (energy_kcal / 1000),
                   ahei_sodium = dplyr::case_when(
                     is.na(sodium_per_1000kcal) ~ NA_real_,
                     sodium_per_1000kcal <= 1000 ~ 10,
                     sodium_per_1000kcal >= 2300 ~ 0,
                     TRUE ~ (2300 - sodium_per_1000kcal) / (2300 - 1000) * 10
                   )
  )

# Alcohol (same heuristic you used before; optionally sex-specific cutpoints if you join DEMO)
ahei_alcohol <- ahei_input_mped %>%
  dplyr::transmute(SEQN,
                   ahei_alcohol = dplyr::case_when(
                     is.na(alcohol_g)        ~ NA_real_,
                     alcohol_g == 0          ~ 2.5,
                     alcohol_g > 0 & alcohol_g <= 13 ~ 10,
                     alcohol_g > 13 & alcohol_g <= 26 ~ (26 - alcohol_g) / (26 - 13) * 10,
                     alcohol_g > 26          ~ 0
                   )
  )

# SSB + juice (MPED gave juice; SSB is NA for now; fill later if you compute from DRXIFF)
ahei_ssb <- ahei_input_mped %>%
  dplyr::transmute(SEQN, ahei_ssb = NA_real_)

# Combine
ahei_9904 <- list(
  ahei_veg, ahei_fruit, ahei_grain, ahei_nutslegumes,
  ahei_meat, ahei_longn3, ahei_pufa, ahei_sodium, ahei_alcohol, ahei_ssb
) %>% purrr::reduce(dplyr::left_join, by = "SEQN") %>%
  dplyr::mutate(ahei_total = rowSums(dplyr::across(dplyr::starts_with("ahei_")), na.rm = TRUE))


# 6) Save + QC -----------------------------------------------------------------

readr::write_csv(ahei_9904, file.path(dir$output, "ahei_1999_2004_day1.csv"))

summary(ahei_9904$ahei_total)

# Quick histogram
ggplot2::ggplot(ahei_9904, ggplot2::aes(x = ahei_total)) +
  ggplot2::geom_histogram(binwidth = 5) +
  ggplot2::labs(title = "AHEI (1999–2004, Day 1)", x = "AHEI total", y = "Count") +
  ggplot2::theme_minimal()


