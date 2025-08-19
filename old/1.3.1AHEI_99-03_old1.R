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

names(pyr_tot_0304)

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















#### try to add individual level DR1IFF 

# 3) NHANES diet totals (1999–2004) + IFF (Day 1) -> merge with MPED totals ----

# Helper to read XPT with a friendly message
safe_read_xpt <- function(path) {
  if (!file.exists(path)) { message("Missing file: ", path); return(NULL) }
  haven::read_xpt(path)
}

# Tiny helpers (no coalesce)
pick_name <- function(nms, candidates) {
  hits <- candidates[candidates %in% nms]
  if (length(hits)) hits[1] else NA_character_
}
mean_preserve_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# ---- A) DR totals (Day 1) as before ------------------------------------------
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

dr_list <- dr_paths %>%
  dplyr::mutate(
    DR1 = purrr::map(dr1, ~ read_dr(.x, 1)),
    DR2 = purrr::map(dr2, ~ read_dr(.x, 2))
  )

dr_all_d1 <- dplyr::bind_rows(dr_list$DR1)

# ---- B) MPED totals join (Day 1) ---------------------------------------------
join_mped <- function(dr_df, mped_df) {
  if (is.null(dr_df) || is.null(mped_df)) return(NULL)
  if ("DAY" %in% names(mped_df)) dplyr::left_join(dr_df, mped_df, by = c("SEQN","DAY"))
  else                           dplyr::left_join(dr_df, mped_df, by = "SEQN")
}

pyr_tot_0304_d1 <- if (!is.null(pyr_tot_d1_0304)) dplyr::mutate(pyr_tot_d1_0304, DAY = 1) else NULL

dr_9900 <- read_dr(file.path(dir$nhanes, "DRXTOT.XPT"),   1)
dr_0102 <- read_dr(file.path(dir$nhanes, "DRXTOT_B.XPT"), 1)
dr_0304 <- read_dr(file.path(dir$nhanes, "DR1TOT_C.XPT"), 1)

dr_9902_mped <- if (!is.null(pyr_tot_9902)) join_mped(dplyr::bind_rows(dr_9900, dr_0102), pyr_tot_9902) else NULL
dr_0304_mped <- if (!is.null(pyr_tot_0304_d1)) join_mped(dr_0304, pyr_tot_0304_d1) else NULL

dr_9904_mped <- dplyr::bind_rows(
  if (!is.null(dr_9902_mped)) dr_9902_mped else dplyr::tibble(),
  if (!is.null(dr_0304_mped)) dr_0304_mped else dplyr::tibble()
)

mped_9904_person <- dr_9904_mped %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarize(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

cat("Persons in MPED 1999–2004 (Day 1 only): ", nrow(mped_9904_person), "\n")


# ---- C) Individual foods (1999–2004): read & standardize ---------

# Example: read just one file
path <- file.path(iff_dir, "DR1IFF_C.xpt")
df <- read_xpt(path)

# See the column names
names(df)

iff_dir <- dir$nhanes

iff_files <- tibble::tribble(
  ~cycle,        ~file,
  "1999-2000",   file.path(iff_dir, "DRXIFF.xpt"),
  "2001-2002",   file.path(iff_dir, "DRXIFF_B.xpt"),
  "2003-2004",   file.path(iff_dir, "DR1IFF_C.xpt")
)

pick_name <- function(nms, candidates) {
  hits <- candidates[candidates %in% nms]
  if (length(hits)) hits[1] else NA_character_
}

read_iff_min <- function(path, cycle_label) {
  df <- haven::read_xpt(path)
  nms <- names(df)
  
  seqn_nm  <- pick_name(nms, c("SEQN"))
  fcode_nm <- pick_name(nms, c("DR1IFDCD","DRDIFDCD","DRXIFDCD","FOODCODE","IFCODE"))
  grams_nm <- pick_name(nms, c("DR1IGRMS","DRXIGRMS","GRAMS"))
  wweia_nm <- pick_name(nms, c("DR1IWWEIA","DRXIWWEIA","WWEIACAT","WWEIACODE"))
  desc_nm  <- pick_name(nms, c("DR1IFDCDTX","DRXIFDCDTX","DR1I_FDNAME","FDNAM","FD_NAME"))
  
  out <- dplyr::tibble(
    SEQN      = if (!is.na(seqn_nm))   df[[seqn_nm]]  else NA_real_,
    FOODCODE  = if (!is.na(fcode_nm))  df[[fcode_nm]] else NA_real_,
    GRAMS     = if (!is.na(grams_nm))  df[[grams_nm]] else NA_real_,
    WWEIA_CAT = if (!is.na(wweia_nm))  df[[wweia_nm]] else NA_real_,
    DESC      = if (!is.na(desc_nm))   as.character(df[[desc_nm]]) else NA_character_,
    DAY       = 1,
    cycle     = cycle_label
  )
  
  if (is.na(wweia_nm)) message("No WWEIA category in ", basename(path), " (expected for 1999–2004).")
  out
}

iff_all <- purrr::map2_dfr(iff_files$file, iff_files$cycle, read_iff_min)
cat("IFF rows stacked (1999–2004): ", nrow(iff_all), "\n")



iff_all$FOODCODE























# ---- D) Define Wang-style flags via WWEIA categories (FILL THESE VECTORS) ----
# >>>>>>> Replace with your WWEIA category codes for 1999–2004 <<<<<<<
W_sandwiches      <- c()  # e.g., 9110, 9120 ...
W_veg_soups       <- c()  # e.g., 7230 ...
W_baby_foods      <- c()  # e.g., 9200:9299 ...
W_veg_with_sauces <- c()  # pasta/rice w/ vegetable sauces ...
W_tomato_mixtures <- c()

W_tomato_juices   <- c()
W_tomato_sauces   <- c()  # includes catsup
W_olives_pickles  <- c()
W_potato_items    <- c()
W_starchy_veg     <- c()  # corn/peas/plantains etc.

W_soft_drinks     <- c()  # regular soda, sweet tea, fruitades
W_fruit_drinks    <- c()  # fruit drinks (not 100% juice)
W_sweetened_other <- c()  # energy/sports/sweetened coffee drinks etc.
W_juice_100       <- c()  # 100% juice

iff_flags <- iff_all %>%
  dplyr::mutate(
    is_half_sandwich   = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_sandwiches,
    is_half_veg_soup   = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_veg_soups,
    is_half_baby_food  = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_baby_foods,
    is_half_with_sauce = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_veg_with_sauces,
    is_half_tom_mix    = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_tomato_mixtures,
    
    is_exc_tom_juice   = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_tomato_juices,
    is_exc_tom_sauce   = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_tomato_sauces,
    is_exc_oliv_pick   = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_olives_pickles,
    is_exc_potato      = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_potato_items,
    is_exc_starchy     = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_starchy_veg,
    
    is_ssb             = !is.na(WWEIA_CAT) & WWEIA_CAT %in% c(W_soft_drinks, W_fruit_drinks, W_sweetened_other),
    is_100pct_juice    = !is.na(WWEIA_CAT) & WWEIA_CAT %in% W_juice_100
  )

# ---- E) SSB + 100% juice servings from IFF (8 oz vs 4 oz) --------------------
iff_bev_person <- iff_flags %>%
  dplyr::transmute(
    SEQN, DAY,
    oz = ifelse(is.na(GRAMS), NA_real_, GRAMS / 28.3495),
    ssb_serv_8oz   = ifelse(!is.na(oz) & is_ssb,          oz / 8, NA_real_),
    juice_serv_4oz = ifelse(!is.na(oz) & is_100pct_juice, oz / 4, NA_real_)
  ) %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarise(
    ssb_serv_8oz   = mean_preserve_na(ssb_serv_8oz),
    juice_serv_4oz = mean_preserve_na(juice_serv_4oz),
    ssb_juice_serv = if (all(is.na(c(ssb_serv_8oz, juice_serv_4oz)))) NA_real_
    else (ifelse(is.na(ssb_serv_8oz), 0, ssb_serv_8oz) + ifelse(is.na(juice_serv_4oz), 0, juice_serv_4oz)),
    .groups = "drop"
  )

# ---- F) Adjusted vegetable cups via per-food MPED/FPED (if available) -------
# Try to read per-food equivalents (MPED 2.0 2003–04 example path)
pyr_food_d1_path <- file.path(dir$fped, "MyPyrEquivDB2", "data", "intakes", "pyr_food_d1.sas7bdat")
pyr_food_d1 <- if (file.exists(pyr_food_d1_path)) safe_read_sas(pyr_food_d1_path) else NULL

veg_adj_person <- {
  if (!is.null(pyr_food_d1)) {
    nms <- names(pyr_food_d1)
    fcode_pf <- pick_name(nms, c("FOODCODE","DR1IFDCD","FDCD"))
    veg_pf   <- pick_name(nms, c("V_TOTAL","VEG_TOT","VEG_CUP_EQ","VEG_CUP"))
    day_pf   <- pick_name(nms, c("DAY"))
    seqn_pf  <- pick_name(nms, c("SEQN"))
    
    if (is.na(fcode_pf) || is.na(veg_pf) || is.na(seqn_pf)) {
      message("Per-food MPED file present but expected columns not found; veg_cup_eq_adj set to NA.")
      dplyr::tibble(SEQN = unique(iff_flags$SEQN), veg_cup_eq_adj = NA_real_)
    } else {
      pf_join <- pyr_food_d1 %>%
        dplyr::transmute(
          SEQN     = .data[[seqn_pf]],
          DAY      = if (!is.na(day_pf)) .data[[day_pf]] else 1,
          FOODCODE = .data[[fcode_pf]],
          VEG_CUPS_FOOD = .data[[veg_pf]]
        ) %>%
        dplyr::filter(DAY == 1)
      
      veg_weighted <- iff_flags %>%
        dplyr::mutate(
          veg_weight = dplyr::case_when(
            is_exc_tom_juice | is_exc_tom_sauce | is_exc_oliv_pick | is_exc_potato | is_exc_starchy ~ 0,
            is_half_sandwich | is_half_veg_soup | is_half_baby_food | is_half_with_sauce | is_half_tom_mix ~ 0.5,
            TRUE ~ 1.0
          )
        )
      
      veg_weighted %>%
        dplyr::inner_join(pf_join, by = c("SEQN","DAY","FOODCODE")) %>%
        dplyr::mutate(
          veg_cup_adj_food = ifelse(is.na(VEG_CUPS_FOOD), NA_real_, VEG_CUPS_FOOD * veg_weight)
        ) %>%
        dplyr::group_by(SEQN) %>%
        dplyr::summarise(veg_cup_eq_adj = mean_preserve_na(veg_cup_adj_food), .groups = "drop")
    }
  } else {
    message("No per-food MPED/FPED equivalents found; veg_cup_eq_adj set to NA.")
    dplyr::tibble(SEQN = unique(iff_flags$SEQN), veg_cup_eq_adj = NA_real_)
  }
}

# ---- G) Save/peek & write intermediate (optional) ----------------------------
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

# Save intermediates (optional)
readr::write_rds(mped_9904_person, file.path(dir$output, "mped_1999_2004_person_day1.rds"))
readr::write_rds(iff_bev_person,   file.path(dir$output, "iff_beverages_day1.rds"))
readr::write_rds(veg_adj_person,   file.path(dir$output, "veg_adjusted_day1.rds"))









# =========================
# 4) Build AHEI inputs (MPED + nutrients), Day 1 only — NO coalesce -------
# =========================

# helper: first existing column name (returns NA_character_ if none)
pick_col <- function(df, candidates) {
  hits <- candidates[candidates %in% names(df)]
  if (length(hits)) hits[1] else NA_character_
}

# helper: mean while preserving NA if all missing
mean_preserve_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# ---- Nutrients from DR totals (already joined in dr_9904_mped) ----
# Standardize the source columns ONCE (no coalesce), then summarize per SEQN.
en_col  <- pick_col(dr_9904_mped, c("DRXTKCAL","DR1TKCAL","DR1IKCAL"))
na_col  <- pick_col(dr_9904_mped, c("DRDTSODI","DR1TSODI","DR1ISODI"))
alc_col <- pick_col(dr_9904_mped, c("DRXTALCO","DR1TALCO","DR1IALCO"))
puf_col <- pick_col(dr_9904_mped, c("DRXTPFAT","DR1TPFAT","DR1IPFAT"))
epa_col <- pick_col(dr_9904_mped, c("DRXTP205","DR1TP205","DR1IP205"))
dha_col <- pick_col(dr_9904_mped, c("DRXTP226","DR1TP226","DR1IP226"))
dpa_col <- pick_col(dr_9904_mped, c("DRXTP225","DR1TP225","DR1IP225"))

dr_std <- dr_9904_mped %>%
  dplyr::transmute(
    SEQN,
    energy_kcal = if (!is.na(en_col))  .data[[en_col]]  else NA_real_,
    sodium_mg   = if (!is.na(na_col))  .data[[na_col]]  else NA_real_,
    alcohol_g   = if (!is.na(alc_col)) .data[[alc_col]] else NA_real_,
    pufa_g      = if (!is.na(puf_col)) .data[[puf_col]] else NA_real_,
    epa_g       = if (!is.na(epa_col)) .data[[epa_col]] else NA_real_,
    dha_g       = if (!is.na(dha_col)) .data[[dha_col]] else NA_real_,
    dpa_g       = if (!is.na(dpa_col)) .data[[dpa_col]] else NA_real_
  )

nutrients_9904 <- dr_std %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarise(
    energy_kcal = mean_preserve_na(energy_kcal),
    sodium_mg   = mean_preserve_na(sodium_mg),
    alcohol_g   = mean_preserve_na(alcohol_g),
    pufa_g      = mean_preserve_na(pufa_g),
    epa_g       = mean_preserve_na(epa_g),
    dha_g       = mean_preserve_na(dha_g),
    dpa_g       = mean_preserve_na(dpa_g),
    .groups = "drop"
  )

# --- helpers for food groups (avoid turning NA into 0)

# Vectorized: if total is NA -> NA; if subtract is NA -> treat as 0
sub_or_keep <- function(total, subtract) {
  ifelse(is.na(total), NA_real_, total - ifelse(is.na(subtract), 0, subtract))
}

# =========================
# MPED food groups -> working inputs (NO coalesce)
# =========================
# Keep valid Day-1 weighted recalls (no coalesce to WTDR4YR), join nutrients, and require kcal present.

names(mped_9904_person)

ahei_input_mped <- mped_9904_person %>%
  dplyr::filter(!is.na(WTDRD1) & WTDRD1 > 0) %>%
  dplyr::left_join(nutrients_9904, by = "SEQN") %>%
  dplyr::filter(!is.na(energy_kcal) & energy_kcal > 0) %>%
  dplyr::transmute(
    SEQN,
    WTDRD1,
    RIAGENDR,  # <-- carry sex forward (1 = Male, 2 = Female)
    
    # --- Vegetables: non-starchy (V_TOTAL - V_POTATO); if V_TOTAL is NA -> NA
    V_TOTAL   = if ("V_TOTAL"   %in% names(mped_9904_person))   V_TOTAL   else NA_real_,
    V_POTATO  = if ("V_POTATO"  %in% names(mped_9904_person))   V_POTATO  else NA_real_,
    veg_cup_eq = sub_or_keep(V_TOTAL, V_POTATO),
    
    # --- Fruits: use F_TOTAL as-is here (if you later get a juice var, subtract it before scoring)
    F_TOTAL    = if ("F_TOTAL"   %in% names(mped_9904_person))   F_TOTAL   else NA_real_,
    fruit_cup_eq = F_TOTAL,
    
    # --- Whole grains: oz-eq -> grams (if G_WHL missing -> NA)
    wholegr_oz_eq = if ("G_WHL" %in% names(mped_9904_person)) G_WHL else NA_real_,
    wholegr_g     = ifelse(is.na(wholegr_oz_eq), NA_real_, wholegr_oz_eq * 28.3495),
    
    # --- Nuts & legumes (keep NA if missing)
    nuts_oz_eq        = if ("M_NUTSD" %in% names(mped_9904_person)) M_NUTSD else NA_real_,   # oz-eq
    legumes_cup_eq    = if ("LEGUMES" %in% names(mped_9904_person)) LEGUMES else NA_real_,   # cups
    nuts_legumes_serv = ifelse(is.na(nuts_oz_eq) | is.na(legumes_cup_eq),
                               ifelse(is.na(nuts_oz_eq) & is.na(legumes_cup_eq), NA_real_,
                                      (ifelse(is.na(nuts_oz_eq), 0, nuts_oz_eq) + ifelse(is.na(legumes_cup_eq), 0, legumes_cup_eq) / 0.5)),
                               nuts_oz_eq + (legumes_cup_eq / 0.5)),
    
    # --- Red + processed meat (oz-eq -> servings; if both missing -> NA)
    red_oz_eq   = if ("M_MEAT"  %in% names(mped_9904_person)) M_MEAT  else NA_real_,
    proc_oz_eq  = if ("M_FRANK" %in% names(mped_9904_person)) M_FRANK else NA_real_,
    redproc_oz_eq = ifelse(is.na(red_oz_eq) & is.na(proc_oz_eq), NA_real_,
                           ifelse(is.na(red_oz_eq), 0, red_oz_eq) + ifelse(is.na(proc_oz_eq), 0, proc_oz_eq)),
    redproc_serv  = ifelse(is.na(redproc_oz_eq), NA_real_, redproc_oz_eq / 3.527),
    
    # --- SSB + 100% fruit juice: placeholder (kept NA so participant drops if component missing)
    ssb_serv         = NA_real_,
    fruit_juice_serv = NA_real_,     # set later from DRXIFF if available
    ssb_juice_serv   = ifelse(is.na(ssb_serv) & is.na(fruit_juice_serv), NA_real_,
                              ifelse(is.na(ssb_serv), 0, ssb_serv) + ifelse(is.na(fruit_juice_serv), 0, fruit_juice_serv)),
    
    # --- Nutrients (carried through)
    energy_kcal, sodium_mg, alcohol_g, pufa_g, epa_g, dha_g
  )

# Sanity messages
cat("AHEI input rows after valid Day-1 & non-missing kcal: ", nrow(ahei_input_mped), "\n")

# mean(ahei_input_mped$veg_cup_eq, na.rm = TRUE) # optional peek


# =========================
# 5) AHEI component scoring — keep NA as NA (no NA->0)-------
# =========================

# ---- helpers ----
lin_pos <- function(x, min0, max10) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (x - min0)/(max10 - min0))) * 10)
lin_rev <- function(x, min10, max0) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (max0 - x)/(max0 - min10))) * 10)

# If you have DEMO (RIAGENDR: 1=Male, 2=Female) merged in ahei_input_mped, this will be used.
# If RIAGENDR is missing, whole-grain score becomes NA to enforce complete-case exclusion.
sex_vec <- if ("RIAGENDR" %in% names(ahei_input_mped)) ahei_input_mped$RIAGENDR else rep(NA_real_, nrow(ahei_input_mped))

# ---- 5.1) Vegetables (servings/day; 0 -> 0 pts; 5 -> 10 pts) ----

# Vegetables (past density spec): max at ≥5 cups/1000 kcal

# 1 serving = 0.5 cup; ensure veg_cup_eq is non-starchy veg
ahei_veg <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    veg_serv = ifelse(is.na(veg_cup_eq), NA_real_, veg_cup_eq / 0.5),
    ahei_veg = lin_pos(veg_serv, 0, 5)
  )

summary(ahei_veg)

# ---- 5.2) Fruit (servings/day; 0 -> 0 pts; 4 -> 10 pts) ----
# 1 serving = 0.5 cup; exclude juice upstream if you have it
ahei_fruit <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    fruit_serv = ifelse(is.na(fruit_cup_eq), NA_real_, fruit_cup_eq / 0.5),
    ahei_fruit = lin_pos(fruit_serv, 0, 4)
  )

summary(ahei_fruit)

# ---- 5.3) Whole grains (g/day; women 75g -> 10, men 90g -> 10) ----
max_wholegr <- ifelse(sex_vec == 2, 75,
                      ifelse(sex_vec == 1, 90, NA_real_))

ahei_grain <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    max_wholegr = max_wholegr,            # for transparency
    ahei_wholegrains = ifelse(is.na(wholegr_g) | is.na(max_wholegr),
                              NA_real_,
                              pmin(wholegr_g / max_wholegr, 1) * 10)
  )

# ---- 5.4) SSB + 100% fruit juice (servings/day; 0 -> 10 pts; ≥1 -> 0 pts) ----
# Expect a single combined servings/day variable (ssb_juice_serv).
# One SSB serving = 8 oz; one 100% juice serving = 4 oz (convert upstream).
ahei_ssb <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    ahei_ssb = lin_rev(ssb_juice_serv, 0, 1)
  )

# ---- 5.5) Nuts & legumes (servings/day; 0 -> 0 pts; 1 -> 10 pts) ----
ahei_nutslegumes <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    ahei_nutslegumes = lin_pos(nuts_legumes_serv, 0, 1)
  )

# ---- 5.6) Red + processed meat (servings/day; 0 -> 10 pts; ≥1.5 -> 0 pts) ----
ahei_meat <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    ahei_redprocmeat = lin_rev(redproc_serv, 0, 1.5)
  )

# ---- 5.7) Long-chain n-3 (EPA + DHA), mg/day; 0 -> 0 pts; 250 mg -> 10 pts ----
ahei_longn3 <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    lc_n3_mg = ifelse(is.na(epa_g) & is.na(dha_g),
                      NA_real_,
                      (ifelse(is.na(epa_g), 0, epa_g) + ifelse(is.na(dha_g), 0, dha_g)) * 1000),
    ahei_longn3 = lin_pos(lc_n3_mg, 0, 250)
  )

# ---- 5.8) PUFA % energy; ≤2%E -> 0 pts; ≥10%E -> 10 pts ----
ahei_pufa <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    pufa_energy_pct = ifelse(is.na(pufa_g) | is.na(energy_kcal), NA_real_, (pufa_g * 9) / energy_kcal * 100),
    ahei_pufa = dplyr::case_when(
      is.na(pufa_energy_pct) ~ NA_real_,
      pufa_energy_pct <= 2   ~ 0,
      pufa_energy_pct >= 10  ~ 10,
      TRUE ~ (pufa_energy_pct - 2) / (10 - 2) * 10
    )
  )

# ---- 5.9) Alcohol (sex-specific J-curve) ----
# Women: 0.5–1.5 drinks (7–21 g) = 10; Men: 0.5–2.0 drinks (7–28 g) = 10;
# 0 or >2.5 drinks (≥35 g women; ≥49 g men) = 0; linear ramps between.
ahei_alcohol <- ahei_input_mped %>%
  dplyr::mutate(sex = sex_vec) %>%
  dplyr::transmute(
    SEQN,
    ahei_alcohol = dplyr::case_when(
      is.na(alcohol_g) | is.na(sex)                 ~ NA_real_,
      # Women
      sex == 2 & alcohol_g <= 0                     ~ 0,
      sex == 2 & alcohol_g > 0  & alcohol_g < 7     ~ (alcohol_g / 7) * 10,
      sex == 2 & alcohol_g >= 7  & alcohol_g <= 21  ~ 10,
      sex == 2 & alcohol_g > 21 & alcohol_g < 35    ~ ((35 - alcohol_g) / (35 - 21)) * 10,
      sex == 2 & alcohol_g >= 35                     ~ 0,
      # Men
      sex == 1 & alcohol_g <= 0                     ~ 0,
      sex == 1 & alcohol_g > 0  & alcohol_g < 7     ~ (alcohol_g / 7) * 10,
      sex == 1 & alcohol_g >= 7  & alcohol_g <= 28  ~ 10,
      sex == 1 & alcohol_g > 28 & alcohol_g < 49    ~ ((49 - alcohol_g) / (49 - 28)) * 10,
      sex == 1 & alcohol_g >= 49                     ~ 0
    )
  )

# ---- 5.10) Sodium — weighted deciles (lowest decile = 10, highest = 0) ----
# Use energy-adjusted residuals (weighted) to form deciles.
library(survey)
sod_df <- ahei_input_mped %>% dplyr::select(SEQN, WTDRD1, sodium_mg, energy_kcal) %>%
  dplyr::filter(!is.na(WTDRD1) & WTDRD1 > 0 & !is.na(sodium_mg) & !is.na(energy_kcal))

des <- svydesign(ids = ~1, weights = ~WTDRD1, data = sod_df)
fit <- svyglm(sodium_mg ~ energy_kcal, design = des)   # residual method
sod_df$sod_resid <- resid(fit, type = "response")

qs <- svyquantile(~sod_resid, design = des, quantiles = seq(0.1, 0.9, by = 0.1), ci = FALSE, interval.type = "quantile")
qs <- as.numeric(qs)

cut_to_decile <- function(x, dec) {
  brks <- c(-Inf, dec, Inf)
  as.numeric(cut(x, breaks = brks, labels = FALSE))
}

sod_df$sod_decile <- cut_to_decile(sod_df$sod_resid, qs)
sod_sc <- sod_df %>% dplyr::transmute(SEQN, ahei_sodium = (11 - sod_decile) / 10 * 10)

ahei_sodium <- ahei_input_mped %>% dplyr::select(SEQN) %>% dplyr::left_join(sod_sc, by = "SEQN")







# =========================
# Combine components -> exclude any missing component -> total
# =========================
ahei_all <- list(
  ahei_veg, ahei_fruit, ahei_grain, ahei_nutslegumes,
  ahei_meat, ahei_longn3, ahei_pufa, ahei_sodium, ahei_alcohol, ahei_ssb
) %>% purrr::reduce(dplyr::left_join, by = "SEQN")

# Exclude participants with ANY missing component (Patel/Wang approach)
ahei_complete <- ahei_all %>%
  dplyr::filter(dplyr::if_all(dplyr::starts_with("ahei_"), ~ !is.na(.))) %>%
  dplyr::mutate(
    ahei_total = rowSums(dplyr::across(dplyr::starts_with("ahei_")), na.rm = FALSE)
  )

message("Excluded ", nrow(ahei_all) - nrow(ahei_complete), " participants due to ≥1 missing component.")







# Save + QC
readr::write_csv(ahei_complete, file.path(dir$output, "ahei_1999_2004_day1.csv"))
summary(ahei_complete$ahei_total)

ggplot2::ggplot(ahei_complete, ggplot2::aes(x = ahei_total)) +
  ggplot2::geom_histogram(binwidth = 5) +
  ggplot2::labs(title = "AHEI (1999–2004, Day 1) — complete cases", x = "AHEI total", y = "Count") +
  ggplot2::theme_minimal()


# 6) Save + QC -----------------------------------------------------------------

readr::write_csv(ahei_9904, file.path(dir$output, "ahei_1999_2004_day1.csv"))

summary(ahei_9904$ahei_total)

# Quick histogram
ggplot2::ggplot(ahei_9904, ggplot2::aes(x = ahei_total)) +
  ggplot2::geom_histogram(binwidth = 5) +
  ggplot2::labs(title = "AHEI (1999–2004, Day 1)", x = "AHEI total", y = "Count") +
  ggplot2::theme_minimal()


