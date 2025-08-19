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
pkgs <- c("dplyr", "haven", "foreign", "survey", "purrr", "ggplot2", "readr", 
          "stringr", "tidyr", "rlang", "janitor")
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


# 4) Individual foods (1999–2004): read & standardize ---------

# Example: read just one file
path <- file.path(iff_dir, "DR1IFF_C.xpt")
df <- read_xpt(path)

# See the column names
names(df)

#  column might be useful
df$DR1IALCO

df$DR1IFDCD


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

names(iff_all)
iff_all$FOODCODE

#####---------- Helper to read one FNDDS ascii folder into a lookup ----------
read_fndds_lookup <- function(ascii_dir,
                              main_file = "MainFoodDesc.txt",
                              add_file  = "AddFoodDesc.txt") {
  
  main_path <- file.path(ascii_dir, main_file)
  add_path  <- file.path(ascii_dir, add_file)
  
  # Main descriptions (caret-delimited, tilde-quoted, trailing caret)
  main <- read_delim(
    main_path,
    delim = "^",
    quote = "~",
    col_names = c("food_code","start_date","end_date",
                  "main_food_description","main_food_description_upper","_extra"),
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
  ) %>%
    select(-`_extra`) %>%
    transmute(
      FOODCODE = as.numeric(food_code),
      DESC_main = main_food_description
    )
  
  # Additional descriptors (if present); collapse per FOODCODE
  add <- if (file.exists(add_path)) {
    read_delim(
      add_path,
      delim = "^",
      quote = "~",
      col_names = c("food_code","start_date","end_date",
                    "add_food_desc","add_food_desc_upper","_extra"),
      escape_double = FALSE,
      trim_ws = TRUE,
      show_col_types = FALSE
    ) %>%
      select(-`_extra`) %>%
      transmute(
        FOODCODE = as.numeric(food_code),
        DESC_add  = add_food_desc
      ) %>%
      group_by(FOODCODE) %>%
      summarise(DESC_add = paste0(unique(na.omit(DESC_add)), collapse = "; "),
                .groups = "drop")
  } else {
    tibble(FOODCODE = numeric(), DESC_add = character())
  }
  
  # Compose final DESC and return only needed cols
  main %>%
    left_join(add, by = "FOODCODE") %>%
    mutate(
      DESC = coalesce(
        if_else(!is.na(DESC_add) &
                  !str_detect(DESC_main, fixed(DESC_add, ignore_case = TRUE)),
                paste(DESC_main, "-", DESC_add), DESC_main),
        DESC_main
      )
    ) %>%
    transmute(FOODCODE, DESC)
}

# ---------- 4.1) FNDDS 2001–2002 ----------
#### source: 
#### https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fndds-download-databases/

ascii_0102 <- "/Users/dengshuyue/Desktop/SDOH/analysis/data/fndds/FNDDS1_ASCII_unpacked/fndds/ascii"
fndds_lookup_0102 <- read_fndds_lookup(ascii_0102) %>%
  mutate(cycle_src = "2001-2002")

# ---------- 4.2) FNDDS 2003–2004 ----------
ascii_0304 <- "/Users/dengshuyue/Desktop/SDOH/analysis/data/fndds/FNDDS2/ascii"
# Note: filenames in FNDDS2 are lowercase; pass them explicitly
fndds_lookup_0304 <- read_fndds_lookup(
  ascii_dir = ascii_0304,
  main_file = "mainfooddesc.txt",
  add_file  = "addfooddesc.txt"
) %>% mutate(cycle_src = "2003-2004")

# ---------- Combine lookups & dedupe ----------
# If the same FOODCODE exists in both, keep the first (here: 2001–2002). Swap order if you prefer 2003–2004.
fndds_lookup <- bind_rows(fndds_lookup_0102, fndds_lookup_0304) %>%
  distinct(FOODCODE, .keep_all = TRUE)

# ---------- Join to IFF (collision-proof) ----------
iff_joined <- iff_all %>%
  # remove any pre-existing DESC/desc_lc/cycle_src to avoid .x/.y suffixing
  select(-any_of(c("DESC", "desc_lc", "cycle_src"))) %>%
  mutate(FOODCODE = as.numeric(FOODCODE)) %>%
  left_join(fndds_lookup %>% select(FOODCODE, DESC), by = "FOODCODE") %>%
  mutate(
    DESC    = coalesce(DESC, ""),
    desc_lc = stringr::str_to_lower(DESC)  # optional helper; safe because never NA
  )

# ---------- (Optional) quick sanity checks ----------
names(iff_joined)
sum(nchar(iff_joined$DESC) > 0)         # rows with matched descriptions
iff_joined %>% filter(DESC == "") %>% distinct(FOODCODE) %>% head()







# ---------- 4.3) Full/Half/Exclusion logic -> per-person servings from IFF (with debug) ----------

# Base rows used by all debug blocks
iff_base <- iff_joined %>%
  dplyr::transmute(
    SEQN, DAY, cycle, FOODCODE, DESC,
    dl    = stringr::str_to_lower(coalesce(DESC, "")),
    g     = coalesce(GRAMS, 0),
    fl_oz = g / 29.5735
  )


##### DEBUG: Fruit flags & servings (keeps solids that are "packed in juice")-----
fruit_flags <- iff_base %>%
  dplyr::mutate(
    is_juice    = stringr::str_detect(dl, stringr::regex("\\bjuice\\b", TRUE)),
    is_cocktail = stringr::str_detect(dl, stringr::regex("cocktail|\\bdrink\\b|beverage|ade\\b|punch\\b|nectar", TRUE)),
    bev_milkct  = stringr::str_detect(dl, stringr::regex("\\bmilk\\b|coffee|tea|cappuccino|latte|mocha|cocoa|hot\\s*chocolate", TRUE)),
    juice_pack_context = stringr::str_detect(
      dl,
      stringr::regex("juice\\s*pack(ed)?|packed\\s*in\\s*(its\\s*|own\\s*)?juice|in\\s*(its\\s*|own\\s*)?juice", TRUE)
    ),
    fruit_full_raw = stringr::str_detect(
      dl,
      stringr::regex(
        paste0(
          "apple|apricot|banana|berry|blueberry|strawberry|raspberry|blackberry|",
          "citrus|orange|grapefruit|tangerine|clementine|lemon|lime|",
          "melon|watermelon|cantaloupe|honeydew|",
          "pear|peach|plum|nectarine|grape|cherry|mango|papaya|kiwi|pineapple|",
          "dried\\s*(fruit|apple|apricot|date|fig|raisin|prune)|",
          "fruit\\s*(mix|mixture)"
        ), TRUE
      )
    ),
    fruit_half_raw = stringr::str_detect(
      dl,
      stringr::regex("fruit.*(mix|mixture).*(with|&).*(cereal|yogurt|granola|nuts|ice cream|dessert)|baby\\s*food.*fruit.*(mix|mixture)", TRUE)
    ),
    fruit_excl = (
      (is_juice & !juice_pack_context & !is_cocktail & !bev_milkct) |
        stringr::str_detect(
          dl,
          stringr::regex(
            "fruit.*(dessert|pie|cobbler|crisp|sherbet|sorbet|pudding|gelatin|gel\\s*dessert)|\\bgelatin\\b|\\bpudding\\b",
            TRUE
          )
        )
    ),
    fruit_class = dplyr::case_when(
      fruit_excl     ~ "excl",
      fruit_half_raw ~ "half",
      fruit_full_raw ~ "full",
      TRUE           ~ "none"
    ),
    # 100% juice beverage flag for SSB+juice (not counted in fruit cups)
    juice100_flag = is_juice & !juice_pack_context & !is_cocktail & !bev_milkct
  )

# Per-person fruit servings (65 g = 1 serving = 0.5 cup)
fruit_person <- fruit_flags %>%
  dplyr::mutate(fruit_serv_contrib = dplyr::case_when(
    fruit_class == "full" ~ g/65,
    fruit_class == "half" ~ (g/65)*0.5,
    TRUE                  ~ 0
  )) %>%
  dplyr::group_by(SEQN, DAY, cycle) %>%
  dplyr::summarise(
    fruit_serv_iff   = sum(fruit_serv_contrib, na.rm = TRUE),
    fruit_cup_eq_iff = fruit_serv_iff * 0.5,
    fruit_juice_serv = sum((fl_oz/4) * as.numeric(juice100_flag), na.rm = TRUE),
    .groups = "drop"
  )

# (optional checks you can run on demand)
# fruit_flags %>% count(fruit_class)
# summary(fruit_person$fruit_serv_iff)
# summary(fruit_person$fruit_cup_eq_iff)



#####  DEBUG: Vegetable flags & servings ------
veg_flags <- iff_base %>%
  dplyr::mutate(
    # FULL: veg + green salads (exclude pasta/egg/tuna/chicken/macaroni/potato salads)
    veg_full_raw = stringr::str_detect(
      dl,
      stringr::regex(
        paste0(
          "spinach|kale|collard|turnip greens|mustard greens|broccoli|romaine|chard|watercress|",
          "carrot|pumpkin|winter squash|butternut|acorn|sweet potato\\b(?!.*fr(y|ies))|",
          "tomato\\b(?!\\s*(juice|sauce))|",
          "asparagus|green bean|string bean|zucchini|summer squash|cabbage|cauliflower|",
          "eggplant|mushroom|pepper|onion|celery|cucumber|lettuce|greens\\b|",
          "(?=.*\\b(salad|garden\\s*salad|green\\s*salad)\\b)(?!.*\\b(pasta|macaroni|potato|tuna|egg|chicken)\\b)"
        ), TRUE
      )
    ),
    # HALF: soups, baby foods, with sauce, casseroles/mixed dishes, coleslaw
    veg_half_raw = stringr::str_detect(
      dl,
      stringr::regex(
        paste0(
          "vegetable.*soup|veg.*soup|baby\\s*food.*(veg|vegetable)|",
          "vegetable.*(sauce|with\\s*sauce)|coleslaw|cole\\s*slaw|",
          "(casserole|stir[- ]?fry|stew|mixed\\s*dish|pasta\\s*with|rice\\s*with).*(broccoli|carrot|pepper|onion|tomato|veg|vegetable)"
        ), TRUE
      )
    ),
    # EXCLUDE: potatoes & starchy veg, tomato juices/sauces, olives/pickles/relishes
    veg_excl = stringr::str_detect(
      dl,
      stringr::regex(
        "potato|french\\s*fry|hash brown|tater tot|\\bcorn\\b|\\blima bean\\b|\\bhominy\\b|\\bplantain\\b|tomato\\s*(juice|sauce)|\\bolive\\b|\\bpickle\\b|\\brelish\\b",
        TRUE
      )
    ),
    veg_class = dplyr::case_when(
      veg_excl     ~ "excl",
      veg_half_raw ~ "half",
      veg_full_raw ~ "full",
      TRUE         ~ "none"
    )
  )

veg_person <- veg_flags %>%
  dplyr::mutate(veg_serv_contrib = dplyr::case_when(
    veg_class == "full" ~ g/65,
    veg_class == "half" ~ (g/65)*0.5,
    TRUE                ~ 0
  )) %>%
  dplyr::summarise(
    veg_serv_iff   = sum(veg_serv_contrib, na.rm = TRUE),
    veg_cup_eq_iff = veg_serv_iff * 0.5,
    .by = c(SEQN, DAY, cycle)
  )

# (optional checks)
# summary(veg_person$veg_cup_eq_iff)


#### OTHER groups (SSB, meats, nuts/legumes, whole grains, alcohol)------
other_person <- iff_base %>%
  dplyr::mutate(
    # -- SSB & 100% juice --
    bev_core    = stringr::str_detect(dl, stringr::regex("soda|cola|soft\\s*drink|\\bpop\\b|lemonade|fruit\\s*(ade|drink|punch)|sports\\s*drink|energy\\s*drink|sweetened\\s*water|smoothie|frappuccino", TRUE)),
    bev_diet    = stringr::str_detect(dl, stringr::regex("\\bdiet\\b|sugar[- ]?free|unsweetened|zero\\b|low\\s*cal", TRUE)),
    bev_milkct  = stringr::str_detect(dl, stringr::regex("\\bmilk\\b|coffee|tea|cappuccino|latte|mocha|cocoa|hot\\s*chocolate", TRUE)),
    bev_reduced = stringr::str_detect(dl, stringr::regex("reduced\\s*sugar|less\\s*sugar|lower\\s*sugar|50%\\s*less\\s*sugar|\\blight\\b", TRUE)),
    ssb_full_flag = bev_core & !bev_diet & !bev_milkct & !bev_reduced,
    ssb_half_flag = bev_core & !bev_diet & !bev_milkct &  bev_reduced,
    
    is_juice    = stringr::str_detect(dl, stringr::regex("\\bjuice\\b", TRUE)),
    is_cocktail = stringr::str_detect(dl, stringr::regex("cocktail|\\bdrink\\b|beverage|ade\\b|punch\\b|nectar", TRUE)),
    juice_pack_context = stringr::str_detect(dl, stringr::regex("juice\\s*pack(ed)?|packed\\s*in\\s*(its\\s*|own\\s*)?juice|in\\s*(its\\s*|own\\s*)?juice", TRUE)),
    juice100_flag = is_juice & !juice_pack_context & !is_cocktail & !bev_diet & !bev_milkct,
    
    # -- RED / PROCESSED MEAT --
    beef_any  = stringr::str_detect(dl, "\\bbeef\\b|\\bhamburger\\b"),
    beef_core = stringr::str_detect(dl, "beef\\s*(steak|rib|ribs|roast|ground|patty|patties|meatball|sirloin|ribeye|brisket)") |
      stringr::str_detect(dl, "\\bhamburger\\s*patty\\b"),
    beef_mix  = beef_any & stringr::str_detect(dl, "baby\\s*food|gravy|sauce|with\\s*(starch|potato|rice|pasta|noodle|vegetable)|sandwich|burger|frozen|tv\\s*dinner|meal|soup|broth|extract"),
    beef_excl = stringr::str_detect(dl, "beef\\s*bacon|frankfurter|wiener|hot\\s*dog|sausage|luncheon|lunch\\s*meat|bologna|salami|pepperoni|potted\\s*meat|meat\\s*spread|spam|pastrami|corned\\s*beef"),
    
    pork_any  = stringr::str_detect(dl, "\\bpork\\b|\\bham\\b|sparerib|spareribs"),
    pork_core = stringr::str_detect(dl, "pork\\s*(chop|steak|cutlet|roast|loin|tenderloin|sparerib|spareribs)|\\bham\\b"),
    pork_mix  = pork_any & stringr::str_detect(dl, "baby\\s*food|gravy|sauce|with\\s*(starch|potato|rice|pasta|noodle|vegetable)|sandwich|frozen|tv\\s*dinner|meal|soup|broth|extract"),
    pork_excl = stringr::str_detect(dl, "canadian\\s*bacon|\\bbacon\\b|salt\\s*pork|pork\\s*skin|frankfurter|wiener|hot\\s*dog|sausage\\b|luncheon|lunch\\s*meat|bologna|salami|pepperoni|potted\\s*meat|meat\\s*spread|spam"),
    
    beef_full = beef_core & !beef_mix & !beef_excl,
    beef_half = beef_mix  & !beef_excl,
    pork_full = pork_core & !pork_mix & !pork_excl,
    pork_half = pork_mix  & !pork_excl,
    
    proc_core = stringr::str_detect(dl, "canadian\\s*bacon|\\bbacon\\b|salt\\s*pork|pork\\s*skin|frankfurter|wiener|hot\\s*dog|sausage\\b|bratwurst|kielbasa|luncheon|lunch\\s*meat|cold\\s*cut|bologna|salami|pepperoni|potted\\s*meat|meat\\s*spread|liverwurst|pastrami|corned\\s*beef|spam|scrapple"),
    proc_half = stringr::str_detect(dl, "(sandwich|sub|hoagie|roll|bun).*(frankfurter|hot\\s*dog|wiener|luncheon|lunch\\s*meat|potted\\s*meat)"),
    
    # -- NUTS & LEGUMES --
    nls_full = stringr::str_detect(
      dl,
      stringr::regex(
        paste0(
          "\\b(black|pinto|kidney|navy|refried|garbanzo|chickpea|white|great northern|cannellini)\\b.*bean|",
          "\\bdried\\s*bean(s)?\\b|\\bbaked\\s*bean(s)?\\b|",
          "\\blentil(s)?\\b|\\bsplit\\s*pea\\b|\\bgreen\\s*pea(s)?\\b|\\bpea soup\\b(?!.*cream)|",
          "almond|walnut|pecan|peanut\\b(?!\\s*butter)|cashew|pistachio|hazelnut|macadamia|brazil\\s*nut|",
          "sunflower\\s*seed|pumpkin\\s*seed|sesame\\s*seed|chia\\s*seed|flax\\s*seed|hemp\\s*seed"
        ), TRUE
      )
    ),
    nls_half = stringr::str_detect(dl, stringr::regex("(frozen|tv\\s*dinner|plate\\s*meal).*(bean|lentil|pea)|bean.*soup|lentil.*soup|pea.*soup|baby\\s*food.*(bean|lentil|pea)", TRUE)),
    nls_excl = stringr::str_detect(dl, stringr::regex("soy\\b|tofu|tempeh|textured\\s*veg|meat\\s*substitute|veggie\\s*burger|nut\\s*butter|coconut\\s*(milk|beverage)|carob", TRUE)),
    nb_full  = stringr::str_detect(dl, stringr::regex("(peanut|almond|cashew|hazelnut|sunflower)\\s*butter\\b", TRUE)),
    nb_half  = stringr::str_detect(dl, stringr::regex("(peanut|nut)\\s*butter\\s*(and\\s*jelly|&\\s*jelly)|nut\\s*butter\\s*sandwich", TRUE)),
    nb_excl  = stringr::str_detect(dl, stringr::regex("nut\\s*gravy|peanut\\s*sauce", TRUE)),
    tofu_full = stringr::str_detect(dl, stringr::regex("tofu|soy\\s*bean\\b|soybeans\\b|soy\\s*nuts\\b", TRUE)),
    tofu_half = stringr::str_detect(dl, stringr::regex("soy\\s*yogurt|soy\\s*dessert|tofu\\s*soup|tofu.*(mix|dish)", TRUE)),
    soymilk_full = stringr::str_detect(dl, stringr::regex("\\bsoy\\s*(milk|drink|beverage)\\b", TRUE)),
    
    # -- WHOLE GRAINS --
    wg_kw = stringr::str_detect(
      dl,
      stringr::regex(
        paste0(
          "100%\\s*whole|\\bwhole\\s*(wheat|grain|rye|oat)|\\bwholegrain\\b|",
          "oatmeal|rolled\\s*oats|steel[- ]cut\\s*oats|bran\\b|bran\\s*flakes|",
          "brown\\s*rice\\b|popcorn\\b|bulgur|quinoa|farro|barley\\b|",
          "(whole|wholegrain).*\\b(bread|tortilla|pita|pasta|noodle|cracker|cereal)\\b"
        ), TRUE
      )
    ),
    
    # -- Alcohol (QC) --
    beer_full    = stringr::str_detect(dl, "\\b(beer|lager|ale|ipa|stout|porter|pilsner)\\b") & !stringr::str_detect(dl,"non[- ]?alcoholic|near\\s*beer|root\\s*beer|ginger\\s*beer"),
    wine_full    = stringr::str_detect(dl, "\\b(wine|merlot|cabernet|pinot|chardonnay|riesling|zinfandel|sauvignon)\\b") & !stringr::str_detect(dl,"cooking|non[- ]?alcoholic"),
    liq_full     = stringr::str_detect(dl, "whisk(e)?y|bourbon|scotch|vodka|rum|gin|tequila|brandy|cognac|liqueur") & !stringr::str_detect(dl,"non[- ]?alcoholic"),
    alc_cocktail = stringr::str_detect(dl, "cocktail|mixed\\s*drink|margarita|mojito|martini|bloody\\s*mary|long\\s*island|cosmopolitan|daiquiri|pi[nn]a\\s*colada|sangria|spritzer|michelada|shandy")
  ) %>%
  dplyr::summarise(
    # SSB + 100% juice
    ssb_serv_full   = sum((fl_oz/8) * as.numeric(ssb_full_flag), na.rm = TRUE),
    ssb_serv_half   = sum((fl_oz/8) * as.numeric(ssb_half_flag), na.rm = TRUE),
    ssb_serv        = ssb_serv_full + 0.5 * ssb_serv_half,
    fruit_juice_serv= sum((fl_oz/4) * as.numeric(juice100_flag), na.rm = TRUE),
    ssb_juice_serv  = ssb_serv + fruit_juice_serv,
    
    # Red + processed meat
    beef_serv_fw = sum(ifelse(beef_full, g/100, 0), na.rm = TRUE),
    beef_serv_hw = sum(ifelse(beef_half, (g/100)*0.5, 0), na.rm = TRUE),
    pork_serv_fw = sum(ifelse(pork_full, g/100, 0), na.rm = TRUE),
    pork_serv_hw = sum(ifelse(pork_half, (g/100)*0.5, 0), na.rm = TRUE),
    red_serv     = beef_serv_fw + beef_serv_hw + pork_serv_fw + pork_serv_hw,
    proc_serv_fw = sum(ifelse(proc_core, g/100, 0), na.rm = TRUE),
    proc_serv_hw = sum(ifelse(proc_half, (g/100)*0.5, 0), na.rm = TRUE),
    redproc_serv_iff = red_serv + proc_serv_fw + proc_serv_hw,
    
    # Nuts & legumes
    nls_serv_full = sum(ifelse(nls_full & !nls_excl, g/50, 0), na.rm = TRUE),
    nls_serv_half = sum(ifelse(nls_half & !nls_excl, (g/50)*0.5, 0), na.rm = TRUE),
    nb_serv_full  = sum(ifelse(nb_full  & !nb_excl,  g/32, 0), na.rm = TRUE),
    nb_serv_half  = sum(ifelse(nb_half  & !nb_excl,  (g/32)*0.5, 0), na.rm = TRUE),
    tofu_serv_full= sum(ifelse(tofu_full, g/50, 0), na.rm = TRUE),
    tofu_serv_half= sum(ifelse(tofu_half, (g/50)*0.5, 0), na.rm = TRUE),
    soymilk_serv_full = sum(ifelse(soymilk_full, g/50, 0), na.rm = TRUE),
    nuts_legumes_serv_iff = nls_serv_full + nls_serv_half +
      nb_serv_full  + nb_serv_half  +
      tofu_serv_full + tofu_serv_half +
      soymilk_serv_full,
    
    # Whole grains (grams)
    wholegr_g_iff = sum(ifelse(wg_kw, g, 0), na.rm = TRUE),
    
    # Alcohol (QC only)
    alc_beer_serv_full = sum(ifelse(beer_full,  g/340.2,  0), na.rm = TRUE),
    alc_wine_serv_full = sum(ifelse(wine_full,  g/141.75, 0), na.rm = TRUE),
    alc_liq_serv_full  = sum(ifelse(liq_full,   g/42.53,  0), na.rm = TRUE),
    alc_cocktail_half  = sum(ifelse(alc_cocktail, (g/42.53)*0.5, 0), na.rm = TRUE),
    alc_drinks_serv    = alc_beer_serv_full + alc_wine_serv_full + alc_liq_serv_full + alc_cocktail_half,
    
    .by = c(SEQN, DAY, cycle)
  )


####  FINAL: join debug fruit & veg into the combined IFF servings ----
iff_servings <- other_person %>%
  dplyr::left_join(veg_person,   by = c("SEQN","DAY","cycle"), suffix = c("", "_vegdbg")) %>%
  dplyr::left_join(fruit_person, by = c("SEQN","DAY","cycle"), suffix = c("", "_fruitdbg")) %>%
  dplyr::mutate(
    # choose the beverage juice from other_person unless you prefer fruit debug’s version
    fruit_juice_serv = dplyr::coalesce(fruit_juice_serv, fruit_juice_serv_fruitdbg)
  )

# then summarise exactly as you tried
iff_servings_seqn <- iff_servings %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarise(
    dplyr::across(
      c(veg_cup_eq_iff, fruit_cup_eq_iff, ssb_serv, fruit_juice_serv, ssb_juice_serv,
        redproc_serv_iff, nuts_legumes_serv_iff, wholegr_g_iff, alc_drinks_serv),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )



# ========= 5) Build AHEI inputs (MPED + nutrients), Day 1 only =========
### EDIT: replace your existing Section 5 block with this

# Helpers
pick_col <- function(df, candidates) {
  hits <- candidates[candidates %in% names(df)]
  if (length(hits)) hits[1] else NA_character_
}
mean_preserve_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sub_or_keep <- function(total, subtract) ifelse(is.na(total), NA_real_, total - ifelse(is.na(subtract), 0, subtract))

# Standardize DR totals once, then average to person-level
en_col  <- pick_col(dr_9904_mped, c("DRXTKCAL","DR1TKCAL","DR1IKCAL"))
na_col  <- pick_col(dr_9904_mped, c("DRDTSODI","DR1TSODI","DR1ISODI"))
alc_col <- pick_col(dr_9904_mped, c("DRXTALCO","DR1TALCO","DR1IALCO"))
puf_col <- pick_col(dr_9904_mped, c("DRXTPFAT","DR1TPFAT","DR1IPFAT"))
epa_col <- pick_col(dr_9904_mped, c("DRXTP205","DR1TP205","DR1IP205"))
dha_col <- pick_col(dr_9904_mped, c("DRXTP226","DR1TP226","DR1IP226"))
dpa_col <- pick_col(dr_9904_mped, c("DRXTP225","DR1TP225","DR1IP225"))

dr_std <- dr_9904_mped %>%
  transmute(
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
  group_by(SEQN) %>%
  summarise(
    energy_kcal = mean_preserve_na(energy_kcal),
    sodium_mg   = mean_preserve_na(sodium_mg),
    alcohol_g   = mean_preserve_na(alcohol_g),
    pufa_g      = mean_preserve_na(pufa_g),
    epa_g       = mean_preserve_na(epa_g),
    dha_g       = mean_preserve_na(dha_g),
    dpa_g       = mean_preserve_na(dpa_g),
    .groups = "drop"
  )

ahei_input_mped <- mped_9904_person %>%
  filter(!is.na(WTDRD1) & WTDRD1 > 0) %>%
  left_join(nutrients_9904,      by = "SEQN") %>%
  filter(!is.na(energy_kcal) & energy_kcal > 0) %>%
  left_join(iff_servings_seqn,   by = "SEQN") %>%   # <<< EDIT: join to collapsed IFF metrics
  transmute(
    SEQN, WTDRD1, RIAGENDR,
    
    # MPED veg/fruit (fallbacks)
    V_TOTAL   = if ("V_TOTAL"   %in% names(mped_9904_person)) V_TOTAL   else NA_real_,
    V_POTATO  = if ("V_POTATO"  %in% names(mped_9904_person)) V_POTATO  else NA_real_,
    F_TOTAL   = if ("F_TOTAL"   %in% names(mped_9904_person)) F_TOTAL   else NA_real_,
    
    # Prefer IFF full/half logic; fallback to MPED
    veg_cup_eq_final   = coalesce(veg_cup_eq_iff, sub_or_keep(V_TOTAL, V_POTATO)),
    fruit_cup_eq_final = coalesce(fruit_cup_eq_iff, F_TOTAL),
    
    # Whole grains — prefer MPED ounces→grams, else IFF keyword grams
    wholegr_oz_eq = if ("G_WHL" %in% names(mped_9904_person)) G_WHL else NA_real_,
    wholegr_g_mped = ifelse(is.na(wholegr_oz_eq), NA_real_, wholegr_oz_eq * 28.3495),
    wholegr_g_final = coalesce(wholegr_g_mped, wholegr_g_iff),
    
    # Nuts & legumes — MPED construction (fallback)
    nuts_oz_eq        = if ("M_NUTSD" %in% names(mped_9904_person)) M_NUTSD else NA_real_,
    legumes_cup_eq    = if ("LEGUMES" %in% names(mped_9904_person)) LEGUMES else NA_real_,
    nuts_legumes_serv_mped = ifelse(
      is.na(nuts_oz_eq) | is.na(legumes_cup_eq),
      ifelse(is.na(nuts_oz_eq) & is.na(legumes_cup_eq), NA_real_,
             (ifelse(is.na(nuts_oz_eq), 0, nuts_oz_eq) + ifelse(is.na(legumes_cup_eq), 0, legumes_cup_eq) / 0.5)),
      nuts_oz_eq + (legumes_cup_eq / 0.5)
    ),
    nuts_legumes_serv_final = coalesce(nuts_legumes_serv_iff, nuts_legumes_serv_mped),
    
    # Red + processed meat — prefer IFF (full/half), else MPED (oz-eq)
    red_oz_eq   = if ("M_MEAT"  %in% names(mped_9904_person)) M_MEAT  else NA_real_,
    proc_oz_eq  = if ("M_FRANK" %in% names(mped_9904_person)) M_FRANK else NA_real_,
    redproc_oz_eq_mped = ifelse(is.na(red_oz_eq) & is.na(proc_oz_eq), NA_real_,
                                ifelse(is.na(red_oz_eq), 0, red_oz_eq) + ifelse(is.na(proc_oz_eq), 0, proc_oz_eq)),
    redproc_serv_mped  = ifelse(is.na(redproc_oz_eq_mped), NA_real_, redproc_oz_eq_mped / 3.527),
    redproc_serv_final = coalesce(redproc_serv_iff, redproc_serv_mped),
    
    # Beverages (from IFF; DO NOT use A_BEV)
    ssb_serv, fruit_juice_serv, ssb_juice_serv,
    
    # Nutrients for scoring
    energy_kcal, sodium_mg, alcohol_g, pufa_g, epa_g, dha_g
  )

# quick sanity (should be one row per SEQN)
cat("ahei_input_mped rows:", nrow(ahei_input_mped),
    " unique SEQN:", dplyr::n_distinct(ahei_input_mped$SEQN), "\n")




# old ========= 6) AHEI scoring (use *_final variables) =========
### EDIT: replace your current scoring joins with this single-mutate version

lin_pos <- function(x, min0, max10) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (x - min0)/(max10 - min0))) * 10)
lin_rev <- function(x, min10, max0) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (max0 - x)/(max0 - min10))) * 10)

ahei_scores <- ahei_input_mped %>%
  mutate(
    # Veg (0.5 cup = 1 serving; max 5/day)
    veg_serv   = ifelse(is.na(veg_cup_eq_final), NA_real_, veg_cup_eq_final / 0.5),
    ahei_veg   = lin_pos(veg_serv, 0, 5),
    
    # Fruit (0.5 cup = 1 serving; max 4/day)
    fruit_serv = ifelse(is.na(fruit_cup_eq_final), NA_real_, fruit_cup_eq_final / 0.5),
    ahei_fruit = lin_pos(fruit_serv, 0, 4),
    
    # Whole grains (sex-specific caps)
    max_wholegr = case_when(RIAGENDR == 2 ~ 75,
                            RIAGENDR == 1 ~ 90,
                            TRUE          ~ NA_real_),
    ahei_wholegrains = ifelse(is.na(wholegr_g_final) | is.na(max_wholegr),
                              NA_real_, pmin(wholegr_g_final / max_wholegr, 1) * 10),
    
    # SSB + 100% juice
    ahei_ssb = lin_rev(ssb_juice_serv, 0, 1),
    
    # Nuts & legumes
    ahei_nutslegumes = lin_pos(nuts_legumes_serv_final, 0, 1),
    
    # Red + processed meat
    ahei_redprocmeat = lin_rev(redproc_serv_final, 0, 1.5),
    
    # Long-chain n-3
    lc_n3_mg = ifelse(is.na(epa_g) & is.na(dha_g), NA_real_,
                      (ifelse(is.na(epa_g),0,epa_g) + ifelse(is.na(dha_g),0,dha_g)) * 1000),
    ahei_longn3 = lin_pos(lc_n3_mg, 0, 250),
    
    # PUFA % energy
    pufa_energy_pct = ifelse(is.na(pufa_g) | is.na(energy_kcal), NA_real_,
                             (pufa_g * 9) / energy_kcal * 100),
    ahei_pufa = case_when(
      is.na(pufa_energy_pct) ~ NA_real_,
      pufa_energy_pct <= 2   ~ 0,
      pufa_energy_pct >= 10  ~ 10,
      TRUE ~ (pufa_energy_pct - 2) / (10 - 2) * 10
    ),
    
    # Alcohol (sex-specific J-curve; uses grams from DR totals)
    ahei_alcohol = case_when(
      is.na(alcohol_g) | is.na(RIAGENDR) ~ NA_real_,
      RIAGENDR == 2 & alcohol_g <= 0                    ~ 0,
      RIAGENDR == 2 & alcohol_g > 0  & alcohol_g < 7    ~ (alcohol_g / 7) * 10,
      RIAGENDR == 2 & alcohol_g >= 7  & alcohol_g <= 21 ~ 10,
      RIAGENDR == 2 & alcohol_g > 21 & alcohol_g < 35   ~ ((35 - alcohol_g) / (35 - 21)) * 10,
      RIAGENDR == 2 & alcohol_g >= 35                   ~ 0,
      RIAGENDR == 1 & alcohol_g <= 0                    ~ 0,
      RIAGENDR == 1 & alcohol_g > 0  & alcohol_g < 7    ~ (alcohol_g / 7) * 10,
      RIAGENDR == 1 & alcohol_g >= 7  & alcohol_g <= 28 ~ 10,
      RIAGENDR == 1 & alcohol_g > 28 & alcohol_g < 49   ~ ((49 - alcohol_g) / (49 - 28)) * 10,
      RIAGENDR == 1 & alcohol_g >= 49                   ~ 0
    )
  ) %>%
  select(SEQN, starts_with("ahei_"))

# Complete cases & total
ahei_complete <- ahei_scores %>%
  filter(if_all(starts_with("ahei_"), ~ !is.na(.))) %>%
  mutate(ahei_total = rowSums(across(starts_with("ahei_")), na.rm = FALSE))



message("Excluded ", nrow(ahei_all) - nrow(ahei_complete), " participants due to ≥1 missing component.")


quick_dup_check <- function(df, nm) {
  dups <- sum(duplicated(df$SEQN))
  cat(nm, ": rows=", nrow(df), " unique SEQN=", dplyr::n_distinct(df$SEQN),
      " dupes=", dups, "\n")
}
quick_dup_check(ahei_input_mped, "ahei_input_mped")




# ========= Means for each AHEI section =========
# ahei_complete: one row per SEQN, all ahei_* present + ahei_total
# ahei_input_mped: has SEQN, WTDRD1 (and RIAGENDR if you want by-sex later)

# Attach day-1 weight (and sex if you want to stratify later)
ahei_comp <- ahei_complete %>%
  left_join(ahei_input_mped %>% select(SEQN, WTDRD1, RIAGENDR), by = "SEQN")

# Components to summarize (include total at the end)
comp_cols <- names(ahei_comp)[grepl("^ahei_", names(ahei_comp))]
comp_cols <- union(setdiff(comp_cols, "ahei_total"), "ahei_total")

# Weighted mean helper (ignores NA in either x or w)
wmean <- function(x, w) {
  i <- !is.na(x) & !is.na(w)
  if (!any(i)) return(NA_real_)
  sum(x[i] * w[i]) / sum(w[i])
}

# Unweighted means
means_unw <- ahei_comp %>%
  summarise(across(all_of(comp_cols), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "component", values_to = "mean_unweighted")

# Weighted means (WTDRD1)
means_w1 <- tibble(
  component = comp_cols,
  mean_weighted = map_dbl(comp_cols, ~ wmean(ahei_comp[[.x]], ahei_comp$WTDRD1))
)

# Weighted means with 6-year combined weight (WTDRD1/3)
wt6 <- ahei_comp$WTDRD1 / 3
means_w6 <- tibble(
  component = comp_cols,
  mean_weighted_6yr = map_dbl(comp_cols, ~ wmean(ahei_comp[[.x]], wt6))
)

# Merge and display
ahei_means <- means_unw %>%
  left_join(means_w1, by = "component") %>%
  left_join(means_w6, by = "component") %>%
  arrange(component)

print(ahei_means, n = nrow(ahei_means))

# ---- Optional: by sex (1=Male, 2=Female)
# ahei_means_by_sex <- ahei_comp %>%
#   group_by(RIAGENDR) %>%
#   summarise(
#     across(all_of(comp_cols), ~ wmean(.x, WTDRD1), .names = "wmean_{.col}")
#   )
# print(ahei_means_by_sex)


# ========= 6) AHEI scoring — one component at a time =========

lin_pos <- function(x, min0, max10) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (x - min0)/(max10 - min0))) * 10)
lin_rev <- function(x, min10, max0) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (max0 - x)/(max0 - min10))) * 10)

# Weighted mean helper (ignores NA in either x or w)
wmean <- function(x, w) {
  i <- !is.na(x) & !is.na(w)
  if (!any(i)) return(NA_real_)
  sum(x[i] * w[i]) / sum(w[i])
}

# Pretty reporter for a single component
report_component <- function(df, score_col, wt_col = "WTDRD1", title = NULL) {
  sc <- df[[score_col]]
  wt <- df[[wt_col]]
  wt6 <- wt / 3  # 6-yr pooled
  
  if (!is.null(title)) cat("\n---", title, "---\n")
  print(summary(sc))
  cat("Mean (unweighted):", mean(sc, na.rm = TRUE), "\n")
  cat("Mean (WTDRD1):    ", wmean(sc, wt), "\n")
  cat("Mean (WTDRD1/3):  ", wmean(sc, wt6), "\n")
}

# ---------------- 6.1 Vegetables ----------------
ahei_veg_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1,
    veg_serv   = ifelse(is.na(veg_cup_eq_final), NA_real_, veg_cup_eq_final ),  
    ahei_veg   = lin_pos(veg_serv, 0, 5)
  )
report_component(ahei_veg_tbl, "ahei_veg", title = "Vegetables")

# ---------------- 6.2 Fruit ----------------
ahei_fruit_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1,
    fruit_serv = ifelse(is.na(fruit_cup_eq_final), NA_real_, fruit_cup_eq_final / 0.5),  # 0.5 cup = 1 serving
    ahei_fruit = lin_pos(fruit_serv, 0, 4)
  )
report_component(ahei_fruit_tbl, "ahei_fruit", title = "Fruit")

# ---------------- 6.3 Whole grains ----------------
ahei_grain_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, RIAGENDR,
    max_wholegr = case_when(RIAGENDR == 2 ~ 75,   # women
                            RIAGENDR == 1 ~ 90,   # men
                            TRUE          ~ NA_real_),
    ahei_wholegrains = ifelse(is.na(wholegr_g_final) | is.na(max_wholegr),
                              NA_real_,
                              pmin(wholegr_g_final / max_wholegr, 1) * 10)
  )
report_component(ahei_grain_tbl, "ahei_wholegrains", title = "Whole grains")

# ---------------- 6.4 SSB + 100% fruit juice ----------------
ahei_ssb_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1,
    ahei_ssb = lin_rev(ssb_juice_serv, 0, 1)   # 0–≥1 servings/day (8 oz SSB, 4 oz 100% juice)
  )
report_component(ahei_ssb_tbl, "ahei_ssb", title = "SSB + 100% juice")

# ---------------- 6.5 Nuts & legumes ----------------
ahei_nutsleg_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1,
    ahei_nutslegumes = lin_pos(nuts_legumes_serv_final, 0, 1)  # servings/day
  )
report_component(ahei_nutsleg_tbl, "ahei_nutslegumes", title = "Nuts & legumes")

# ---------------- 6.6 Red + processed meat ----------------
ahei_meat_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1,
    ahei_redprocmeat = lin_rev(redproc_serv_final, 0, 1.5)    # servings/day
  )
report_component(ahei_meat_tbl, "ahei_redprocmeat", title = "Red + processed meat")

# ---------------- 6.7 Long-chain n-3 (EPA + DHA) ----------------
ahei_longn3_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, epa_g, dha_g,
    lc_n3_mg = ifelse(is.na(epa_g) & is.na(dha_g),
                      NA_real_,
                      (ifelse(is.na(epa_g), 0, epa_g) + ifelse(is.na(dha_g), 0, dha_g)) * 1000),
    ahei_longn3 = lin_pos(lc_n3_mg, 0, 250)
  )
report_component(ahei_longn3_tbl, "ahei_longn3", title = "Long-chain n-3 (EPA+DHA)")

# ---------------- 6.8 PUFA % energy ----------------
ahei_pufa_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, pufa_g, energy_kcal,
    pufa_energy_pct = ifelse(is.na(pufa_g) | is.na(energy_kcal), NA_real_, (pufa_g * 9) / energy_kcal * 100),
    ahei_pufa = case_when(
      is.na(pufa_energy_pct) ~ NA_real_,
      pufa_energy_pct <= 2   ~ 0,
      pufa_energy_pct >= 10  ~ 10,
      TRUE ~ (pufa_energy_pct - 2) / (10 - 2) * 10
    )
  )
report_component(ahei_pufa_tbl, "ahei_pufa", title = "PUFA % energy")

# ---------------- 6.9 Alcohol ----------------
ahei_alcohol_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, RIAGENDR, alcohol_g,
    ahei_alcohol = case_when(
      is.na(alcohol_g) | is.na(RIAGENDR) ~ NA_real_,
      # Women
      RIAGENDR == 2 & alcohol_g <= 0                    ~ 0,
      RIAGENDR == 2 & alcohol_g > 0  & alcohol_g < 7    ~ (alcohol_g / 7) * 10,
      RIAGENDR == 2 & alcohol_g >= 7  & alcohol_g <= 21 ~ 10,
      RIAGENDR == 2 & alcohol_g > 21 & alcohol_g < 35   ~ ((35 - alcohol_g) / (35 - 21)) * 10,
      RIAGENDR == 2 & alcohol_g >= 35                   ~ 0,
      # Men
      RIAGENDR == 1 & alcohol_g <= 0                    ~ 0,
      RIAGENDR == 1 & alcohol_g > 0  & alcohol_g < 7    ~ (alcohol_g / 7) * 10,
      RIAGENDR == 1 & alcohol_g >= 7  & alcohol_g <= 28 ~ 10,
      RIAGENDR == 1 & alcohol_g > 28 & alcohol_g < 49   ~ ((49 - alcohol_g) / (49 - 28)) * 10,
      RIAGENDR == 1 & alcohol_g >= 49                   ~ 0
    )
  )
report_component(ahei_alcohol_tbl, "ahei_alcohol", title = "Alcohol")

ahei_input_mped$alcohol_g




# ---- 6.10) Sodium — energy-adjusted & outlier-handled (AHEI decile scoring) ----
# Approach:
# 1) Compute sodium density (mg / 1000 kcal).
# 2) Filter implausible energy days, then winsorize sodium density at weighted 0.5% / 99.5%.
# 3) Compute weighted deciles on the winsorized density.
# 4) Assign 10..0 points by decile (lower sodium density = better).
#    (Also computes a continuous 10↘0 score between Q10 and Q90 for sensitivity checks.)

# --- helper: weighted quantiles (no external packages)
wtd_quantile <- function(x, w, probs = seq(0.1, 0.9, 0.1)) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

# --- choose weights (Day 1)
wt_vec <- if ("WTDRD1" %in% names(ahei_input_mped)) ahei_input_mped$WTDRD1 else rep(1, nrow(ahei_input_mped))

# --- base with density; drop obvious bad energy days (adjust bounds if desired)
sod_base <- ahei_input_mped %>%
  transmute(
    SEQN, sodium_mg, energy_kcal, WTDRD1 = wt_vec
  ) %>%
  filter(!is.na(sodium_mg), !is.na(energy_kcal), energy_kcal > 0) %>%
  # Exclude implausible Day-1 energy recalls (tune these if your SOP differs)
  filter(energy_kcal >= 500, energy_kcal <= 6000) %>%
  mutate(
    sod_den = sodium_mg / (energy_kcal / 1000)    # mg per 1000 kcal
  )

# --- winsorize sodium density at weighted 0.5% and 99.5% to tame outliers
TRIM_P <- 0.005
trim_lo <- wtd_quantile(sod_base$sod_den, sod_base$WTDRD1, probs = TRIM_P)
trim_hi <- wtd_quantile(sod_base$sod_den, sod_base$WTDRD1, probs = 1 - TRIM_P)
sod_base <- sod_base %>%
  mutate(sod_den_w = pmin(pmax(sod_den, trim_lo), trim_hi))

cat("Sodium density (mg/1000 kcal) winsorization bounds:\n")
print(round(c(lo = trim_lo, hi = trim_hi), 2))

# --- weighted deciles on winsorized density
sod_dec <- wtd_quantile(sod_base$sod_den_w, sod_base$WTDRD1, probs = seq(0.1, 0.9, 0.1))
names(sod_dec) <- paste0("Q", seq(10, 90, 10))
cat("AHEI sodium density deciles (mg/1000 kcal, weighted, winsorized):\n")
print(round(sod_dec, 2))

# --- make a lookup tibble for joining back
sod_scores <- sod_base %>%
  mutate(
    # Discrete decile scoring (AHEI-2010 style): lowest decile = 10, highest decile = 0
    ahei_sodium_dec = case_when(
      sod_den_w <= sod_dec[1] ~ 10,
      sod_den_w <= sod_dec[2] ~ 9,
      sod_den_w <= sod_dec[3] ~ 8,
      sod_den_w <= sod_dec[4] ~ 7,
      sod_den_w <= sod_dec[5] ~ 6,
      sod_den_w <= sod_dec[6] ~ 5,
      sod_den_w <= sod_dec[7] ~ 4,
      sod_den_w <= sod_dec[8] ~ 3,
      sod_den_w <= sod_dec[9] ~ 2,
      TRUE                    ~ 0
    ),
    # Optional: continuous scoring between Q10 (10 pts) and Q90 (0 pts)
    ahei_sodium_cont = case_when(
      is.na(sod_den_w)                 ~ NA_real_,
      sod_den_w <= sod_dec[1]          ~ 10,
      sod_den_w >= sod_dec[9]          ~ 0,
      TRUE ~ 10 - ( (sod_den_w - sod_dec[1]) / (sod_dec[9] - sod_dec[1]) ) * 10
    )
  ) %>%
  select(SEQN, ahei_sodium = ahei_sodium_dec, ahei_sodium_cont)

# --- final output used in your component joins
ahei_sodium_tbl <- sod_scores %>% select(SEQN, ahei_sodium)

# QC
summary(sod_base$sod_den)     # raw density
summary(sod_base$sod_den_w)   # winsorized density
summary(ahei_sodium_tbl$ahei_sodium)


# ========= Join all components & compute totals =========
ahei_scores <- list(
  ahei_veg_tbl %>% select(SEQN, ahei_veg),
  ahei_fruit_tbl %>% select(SEQN, ahei_fruit),
  ahei_grain_tbl %>% select(SEQN, ahei_wholegrains),
  ahei_nutsleg_tbl %>% select(SEQN, ahei_nutslegumes),
  ahei_meat_tbl %>% select(SEQN, ahei_redprocmeat),
  ahei_longn3_tbl %>% select(SEQN, ahei_longn3),
  ahei_pufa_tbl %>% select(SEQN, ahei_pufa),
  ahei_sodium_tbl %>% select(SEQN, ahei_sodium),
  ahei_ssb_tbl %>% select(SEQN, ahei_ssb),
  ahei_alcohol_tbl %>% select(SEQN, ahei_alcohol)
) %>% reduce(left_join, by = "SEQN")


# Complete cases & total
ahei_complete <- ahei_scores %>%
  filter(if_all(starts_with("ahei_"), ~ !is.na(.))) %>%
  mutate(ahei_total = rowSums(across(starts_with("ahei_")), na.rm = FALSE))

cat("\nParticipants with complete AHEI components: ", nrow(ahei_complete), "\n")

# Optional: quick overall summaries
print(summary(ahei_complete$ahei_total))





################ !!!!!!!! ------
# 6) AHEI scoring — energy-adjusted (per 1000 kcal) =========

lin_pos <- function(x, min0, max10) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (x - min0)/(max10 - min0))) * 10)
lin_rev <- function(x, min10, max0) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (max0 - x)/(max0 - min10))) * 10)

# Weighted mean helper (ignores NA in either x or w)
wmean <- function(x, w) {
  i <- !is.na(x) & !is.na(w)
  if (!any(i)) return(NA_real_)
  sum(x[i] * w[i]) / sum(w[i])
}

# Reporter for a single component
report_component <- function(df, score_col, wt_col = "WTDRD1", title = NULL) {
  sc <- df[[score_col]]
  wt <- df[[wt_col]]
  wt6 <- wt / 3  # 6-yr pooled
  if (!is.null(title)) cat("\n---", title, "---\n")
  print(summary(sc))
  cat("Mean (unweighted):", mean(sc, na.rm = TRUE), "\n")
  cat("Mean (WTDRD1):    ", wmean(sc, wt), "\n")
  cat("Mean (WTDRD1/3):  ", wmean(sc, wt6), "\n")
}

# ----- Energy-density helper & reference -----
E_REF <- 2000
per_1000 <- function(x, kcal) ifelse(is.na(x) | is.na(kcal) | kcal <= 0, NA_real_, x / (kcal / 1000))

# ================= Veg (cups/1000 kcal; max at 1.25 cups/1000) =================
ahei_veg_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, energy_kcal,
    veg_cups_per_1000 = per_1000(veg_cup_eq_final, energy_kcal),
    ahei_veg = ifelse(
      is.na(veg_cups_per_1000), NA_real_,
      pmin(veg_cups_per_1000 / (2.5 / (E_REF/1000)), 1) * 10  # 2.5 cups/day @ 2000kcal = 1.25 cups/1000
    )
  )
report_component(ahei_veg_tbl, "ahei_veg", title = "Vegetables (per 1000 kcal)")

# ================= Fruit (cups/1000 kcal; max at 1.0 cups/1000) =================
ahei_fruit_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, energy_kcal,
    fruit_cups_per_1000 = per_1000(fruit_cup_eq_final, energy_kcal),
    ahei_fruit = ifelse(
      is.na(fruit_cups_per_1000), NA_real_,
      pmin(fruit_cups_per_1000 / (2.0 / (E_REF/1000)), 1) * 10  # 2.0 cups/day @ 2000kcal = 1.0 cups/1000
    )
  )
report_component(ahei_fruit_tbl, "ahei_fruit", title = "Fruit (per 1000 kcal)")

# ================= Whole grains (g/1000 kcal; women 37.5, men 45) =================
ahei_grain_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, RIAGENDR, energy_kcal,
    wholegr_g_per_1000 = per_1000(wholegr_g_final, energy_kcal),
    max_wholegr_per_1000 = case_when(
      RIAGENDR == 2 ~ 75 / (E_REF/1000),  # 37.5 g/1000
      RIAGENDR == 1 ~ 90 / (E_REF/1000),  # 45.0 g/1000
      TRUE          ~ NA_real_
    ),
    ahei_wholegrains = ifelse(
      is.na(wholegr_g_per_1000) | is.na(max_wholegr_per_1000),
      NA_real_,
      pmin(wholegr_g_per_1000 / max_wholegr_per_1000, 1) * 10
    )
  )
report_component(ahei_grain_tbl, "ahei_wholegrains", title = "Whole grains (per 1000 kcal)")

# ================= SSB + 100% juice (servings/1000 kcal; max penalty at ≥0.5) =================
ahei_ssb_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, energy_kcal,
    ssb_juice_per_1000 = per_1000(ssb_juice_serv, energy_kcal),
    ahei_ssb = lin_rev(ssb_juice_per_1000, 0, 1 / (E_REF/1000))  # 1.0/day @ 2000kcal => 0.5/1000
  )
report_component(ahei_ssb_tbl, "ahei_ssb", title = "SSB + 100% juice (per 1000 kcal)")

# ================= Nuts & legumes (servings/1000 kcal; max at 0.5) =================
ahei_nutsleg_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, energy_kcal,
    nuts_leg_per_1000 = per_1000(nuts_legumes_serv_final, energy_kcal),
    ahei_nutslegumes = lin_pos(nuts_leg_per_1000, 0, 1 / (E_REF/1000) / 2)  # 1/day @ 2000kcal = 0.5/1000
  )
report_component(ahei_nutsleg_tbl, "ahei_nutslegumes", title = "Nuts & legumes (per 1000 kcal)")

# ================= Red + processed meat (servings/1000 kcal; min at 0.75) =================
# 🔥 this is higher need check 
ahei_meat_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, energy_kcal,
    redproc_per_1000 = per_1000(redproc_serv_final, energy_kcal),
    ahei_redprocmeat = lin_rev(redproc_per_1000, 0, 1.5 / (E_REF/1000))  # 1.5/day @ 2000kcal = 0.75/1000
  )
report_component(ahei_meat_tbl, "ahei_redprocmeat", title = "Red + processed meat (per 1000 kcal)")


# ================= Long-chain n-3 (EPA + DHA) =================
# AHEI standard uses absolute mg/day (kept here). If you prefer density, uncomment the *per_1000* lines.
ahei_longn3_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, epa_g, dha_g, energy_kcal,
    lc_n3_mg = ifelse(is.na(epa_g) & is.na(dha_g),
                      NA_real_,
                      (ifelse(is.na(epa_g), 0, epa_g) + ifelse(is.na(dha_g), 0, dha_g)) * 1000),
    # lc_n3_mg_per_1000 = per_1000(lc_n3_mg, energy_kcal),                 # OPTIONAL density
    ahei_longn3 = lin_pos(lc_n3_mg, 0, 250)                                 # standard
    # ahei_longn3 = lin_pos(lc_n3_mg_per_1000, 0, 250 / (E_REF/1000))      # OPTIONAL density (125 mg/1000)
  )
report_component(ahei_longn3_tbl, "ahei_longn3", title = "Long-chain n-3 (EPA+DHA)")

# ================= PUFA % energy (already energy-adjusted) =================
ahei_pufa_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, pufa_g, energy_kcal,
    pufa_energy_pct = ifelse(is.na(pufa_g) | is.na(energy_kcal), NA_real_, (pufa_g * 9) / energy_kcal * 100),
    ahei_pufa = case_when(
      is.na(pufa_energy_pct) ~ NA_real_,
      pufa_energy_pct <= 2   ~ 0,
      pufa_energy_pct >= 10  ~ 10,
      TRUE ~ (pufa_energy_pct - 2) / (10 - 2) * 10
    )
  )
report_component(ahei_pufa_tbl, "ahei_pufa", title = "PUFA % energy")

# ================= Alcohol (AHEI J-curve with 2.5 for nondrinkers) =================
ahei_alcohol_tbl <- ahei_input_mped %>%
  transmute(
    SEQN, WTDRD1, RIAGENDR, alcohol_g,
    ahei_alcohol = dplyr::case_when(
      is.na(alcohol_g) | is.na(RIAGENDR) ~ NA_real_,
      
      # ---------- Women (RIAGENDR == 2) ----------
      # Nondrinkers get 2.5
      RIAGENDR == 2 & alcohol_g <= 0                     ~ 2.5,
      # Ramp up from 2.5 at 0 g to 10 at 7 g (0.5 drink)
      RIAGENDR == 2 & alcohol_g > 0  & alcohol_g < 7     ~ 2.5 + (alcohol_g / 7) * (10 - 2.5),
      # Optimal window 7–21 g (0.5–1.5 drinks)
      RIAGENDR == 2 & alcohol_g >= 7  & alcohol_g <= 21  ~ 10,
      # Ramp down from 10 at 21 g to 0 at 35 g (≥2.5 drinks)
      RIAGENDR == 2 & alcohol_g > 21 & alcohol_g < 35    ~ ((35 - alcohol_g) / (35 - 21)) * 10,
      # Heavy drinking
      RIAGENDR == 2 & alcohol_g >= 35                    ~ 0,
      
      # ---------- Men (RIAGENDR == 1) ----------
      # Nondrinkers get 2.5
      RIAGENDR == 1 & alcohol_g <= 0                     ~ 2.5,
      # Ramp up from 2.5 at 0 g to 10 at 7 g (0.5 drink)
      RIAGENDR == 1 & alcohol_g > 0  & alcohol_g < 7     ~ 2.5 + (alcohol_g / 7) * (10 - 2.5),
      # Optimal window 7–28 g (0.5–2.0 drinks)
      RIAGENDR == 1 & alcohol_g >= 7  & alcohol_g <= 28  ~ 10,
      # Ramp down from 10 at 28 g to 0 at 49 g (≥2.5 drinks)
      RIAGENDR == 1 & alcohol_g > 28 & alcohol_g < 49    ~ ((49 - alcohol_g) / (49 - 28)) * 10,
      # Heavy drinking
      RIAGENDR == 1 & alcohol_g >= 49                    ~ 0
    )
  )

# (Optional) quick summary with your helper
report_component(ahei_alcohol_tbl, "ahei_alcohol", title = "Alcohol (AHEI w/ 2.5 for nondrinkers)")




# ---- 6.10) Sodium — energy-adjusted & outlier-handled (AHEI decile scoring) ----
# Approach:
# 1) Compute sodium density (mg / 1000 kcal).
# 2) Filter implausible energy days, then winsorize sodium density at weighted 0.5% / 99.5%.
# 3) Compute weighted deciles on the winsorized density.
# 4) Assign 10..0 points by decile (lower sodium density = better).
#    (Also computes a continuous 10↘0 score between Q10 and Q90 for sensitivity checks.)

# --- helper: weighted quantiles (no external packages)
wtd_quantile <- function(x, w, probs = seq(0.1, 0.9, 0.1)) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}

# --- choose weights (Day 1)
wt_vec <- if ("WTDRD1" %in% names(ahei_input_mped)) ahei_input_mped$WTDRD1 else rep(1, nrow(ahei_input_mped))

# --- base with density; drop obvious bad energy days (adjust bounds if desired)
sod_base <- ahei_input_mped %>%
  transmute(
    SEQN, sodium_mg, energy_kcal, WTDRD1 = wt_vec
  ) %>%
  filter(!is.na(sodium_mg), !is.na(energy_kcal), energy_kcal > 0) %>%
  # Exclude implausible Day-1 energy recalls (tune these if your SOP differs)
  filter(energy_kcal >= 500, energy_kcal <= 6000) %>%
  mutate(
    sod_den = sodium_mg / (energy_kcal / 1000)    # mg per 1000 kcal
  )

# --- winsorize sodium density at weighted 0.5% and 99.5% to tame outliers
TRIM_P <- 0.005
trim_lo <- wtd_quantile(sod_base$sod_den, sod_base$WTDRD1, probs = TRIM_P)
trim_hi <- wtd_quantile(sod_base$sod_den, sod_base$WTDRD1, probs = 1 - TRIM_P)
sod_base <- sod_base %>%
  mutate(sod_den_w = pmin(pmax(sod_den, trim_lo), trim_hi))

cat("Sodium density (mg/1000 kcal) winsorization bounds:\n")
print(round(c(lo = trim_lo, hi = trim_hi), 2))

# --- weighted deciles on winsorized density
sod_dec <- wtd_quantile(sod_base$sod_den_w, sod_base$WTDRD1, probs = seq(0.1, 0.9, 0.1))
names(sod_dec) <- paste0("Q", seq(10, 90, 10))
cat("AHEI sodium density deciles (mg/1000 kcal, weighted, winsorized):\n")
print(round(sod_dec, 2))

# --- make a lookup tibble for joining back
sod_scores <- sod_base %>%
  mutate(
    # Discrete decile scoring (AHEI-2010 style): lowest decile = 10, highest decile = 0
    ahei_sodium_dec = case_when(
      sod_den_w <= sod_dec[1] ~ 10,
      sod_den_w <= sod_dec[2] ~ 9,
      sod_den_w <= sod_dec[3] ~ 8,
      sod_den_w <= sod_dec[4] ~ 7,
      sod_den_w <= sod_dec[5] ~ 6,
      sod_den_w <= sod_dec[6] ~ 5,
      sod_den_w <= sod_dec[7] ~ 4,
      sod_den_w <= sod_dec[8] ~ 3,
      sod_den_w <= sod_dec[9] ~ 2,
      TRUE                    ~ 0
    ),
    # Optional: continuous scoring between Q10 (10 pts) and Q90 (0 pts)
    ahei_sodium_cont = case_when(
      is.na(sod_den_w)                 ~ NA_real_,
      sod_den_w <= sod_dec[1]          ~ 10,
      sod_den_w >= sod_dec[9]          ~ 0,
      TRUE ~ 10 - ( (sod_den_w - sod_dec[1]) / (sod_dec[9] - sod_dec[1]) ) * 10
    )
  ) %>%
  select(SEQN, ahei_sodium = ahei_sodium_dec, ahei_sodium_cont)

# --- final output used in your component joins
ahei_sodium_tbl <- sod_scores %>% select(SEQN, ahei_sodium)

# QC
summary(sod_base$sod_den)     # raw density
summary(sod_base$sod_den_w)   # winsorized density
summary(ahei_sodium_tbl$ahei_sodium)

# ---- If you want a single data frame of scores (energy-adjusted where applicable) ----
ahei_scores <- list(
  ahei_veg_tbl %>% select(SEQN, ahei_veg),
  ahei_fruit_tbl %>% select(SEQN, ahei_fruit),
  ahei_grain_tbl %>% select(SEQN, ahei_wholegrains),
  ahei_ssb_tbl %>% select(SEQN, ahei_ssb),
  ahei_nutsleg_tbl %>% select(SEQN, ahei_nutslegumes),
  ahei_meat_tbl %>% select(SEQN, ahei_redprocmeat),
  ahei_longn3_tbl %>% select(SEQN, ahei_longn3),
  ahei_pufa_tbl %>% select(SEQN, ahei_pufa),
  ahei_alcohol_tbl %>% select(SEQN, ahei_alcohol),
  ahei_sodium_tbl %>% select(SEQN, ahei_sodium)
) %>% purrr::reduce(dplyr::left_join, by = "SEQN")



# Complete cases & total
ahei_complete <- ahei_scores %>%
  filter(if_all(starts_with("ahei_"), ~ !is.na(.))) %>%
  mutate(ahei_total = rowSums(across(starts_with("ahei_")), na.rm = FALSE))

cat("\nParticipants with complete AHEI components: ", nrow(ahei_complete), "\n")

# Optional: quick overall summaries
print(summary(ahei_complete$ahei_total))


#### print mean for all groups ====
# ---- helper: weighted mean ----
wmean <- function(x, w) {
  i <- !is.na(x) & !is.na(w)
  if (!any(i)) return(NA_real_)
  sum(x[i] * w[i]) / sum(w[i])
}

# ---- which component columns to summarize ----
comp_cols <- names(ahei_complete)[grepl("^ahei_", names(ahei_complete))]
comp_cols <- setdiff(comp_cols, "ahei_total")  # handle total separately

# ---- attach Day-1 weights ----
ahei_comp_w <- ahei_complete %>%
  dplyr::left_join(ahei_input_mped %>% dplyr::select(SEQN, WTDRD1), by = "SEQN")

# ---- compute means for each component ----
means_unw <- colMeans(ahei_comp_w[, comp_cols, drop = FALSE], na.rm = TRUE)
means_w1  <- sapply(comp_cols, function(cc) wmean(ahei_comp_w[[cc]], ahei_comp_w$WTDRD1))
means_w6  <- sapply(comp_cols, function(cc) wmean(ahei_comp_w[[cc]], ahei_comp_w$WTDRD1 / 3))

mean_table <- tibble::tibble(component = comp_cols) %>%
  dplyr::mutate(
    mean_unweighted = as.numeric(means_unw[component]),
    mean_wt_day1    = as.numeric(means_w1[component]),
    mean_wt_6yr     = as.numeric(means_w6[component])
  ) %>%
  dplyr::arrange(component)

# ---- add total row ----
mean_total <- tibble::tibble(
  component        = "ahei_total",
  mean_unweighted  = mean(ahei_comp_w$ahei_total, na.rm = TRUE),
  mean_wt_day1     = wmean(ahei_comp_w$ahei_total, ahei_comp_w$WTDRD1),
  mean_wt_6yr      = wmean(ahei_comp_w$ahei_total, ahei_comp_w$WTDRD1 / 3)
)

mean_table <- dplyr::bind_rows(mean_table, mean_total)

print(mean_table, n = nrow(mean_table))


# 7) Save + QC -----------------------------------------------------------------

readr::write_csv(ahei_complete, file.path(dir$output, "ahei_1999_2004_day1.csv"))

summary(ahei_complete$ahei_total)

summary(ahei_complete$ahei_total)

# Quick histogram
ggplot2::ggplot(ahei_complete, ggplot2::aes(x = ahei_total)) +
  ggplot2::geom_histogram(binwidth = 5) +
  ggplot2::labs(title = "AHEI (1999–2004, Day 1)", x = "AHEI total", y = "Count") +
  ggplot2::theme_minimal()


