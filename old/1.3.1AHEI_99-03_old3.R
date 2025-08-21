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


# ---------- 4.3) Derive SSB + 100% fruit juice from IFF ----------

# ---------- 4.3) Full/Half/Exclusion logic -> per-person servings from IFF ----------

iff_servings <- iff_joined %>%
  mutate(
    dl    = str_to_lower(coalesce(DESC, "")),
    g     = coalesce(GRAMS, 0),
    fl_oz = g / 29.5735
  ) %>%
  mutate(
    # =======================
    # VEGETABLES (65 g = 1 serving = 0.5 cup)
    # Full: dark-green, deep-yellow, tomatoes (raw/cooked), other veg (raw/cooked)
    veg_full = str_detect(
      dl,
      regex(
        paste0(
          # dark-green
          "spinach|kale|collard|turnip greens|mustard greens|broccoli|romaine|chard|watercress|",
          # deep-yellow
          "carrot|pumpkin|winter squash|butternut|acorn|sweet potato\\b(?!.*fr(y|ies))|",
          # tomatoes raw/cooked (not juice/sauce)
          "tomato\\b(?!\\s*(juice|sauce))|",
          # other common veg
          "asparagus|green bean|string bean|zucchini|summer squash|cabbage|cauliflower|",
          "eggplant|mushroom|pepper|onion|celery|cucumber|lettuce|greens\\b"
        ),
        ignore_case = TRUE
      )
    ),
    # Half: tomato mixtures/sandwiches, vegetables with sauces, vegetable soups, baby foods
    veg_half = str_detect(
      dl,
      regex(
        paste0(
          "tomato.*(mix|sandwich)|vegetable.*(sauce|with\\s*sauce)|",
          "vegetable.*soup|veg.*soup|baby\\s*food.*(veg|vegetable)"
        ),
        ignore_case = TRUE
      )
    ),
    # Exclude: potatoes/starchy veg, tomato juices/sauces, olives, pickles, relishes
    veg_excl = str_detect(
      dl,
      regex(
        paste0(
          "potato|french\\s*fry|hash brown|tater tot|",
          "\\bcorn\\b|\\bgreen pea\\b|\\blima bean\\b|\\bhominy\\b|\\bplantain\\b|",
          "tomato\\s*(juice|sauce)|\\bolive\\b|\\bpickle\\b|\\brelish\\b"
        ),
        ignore_case = TRUE
      )
    ),
    
    # =======================
    # FRUIT (65 g = 1 serving = 0.5 cup)
    # Full: citrus, dried fruits, other fruits, berries, fruit mixtures
    fruit_full = str_detect(
      dl,
      regex(
        paste0(
          "apple|apricot|banana|berry|blueberry|strawberry|raspberry|blackberry|",
          "citrus|orange|grapefruit|tangerine|clementine|lemon|lime|",
          "melon|watermelon|cantaloupe|honeydew|",
          "pear|peach|plum|nectarine|grape|cherry|mango|papaya|kiwi|pineapple|",
          "dried\\s*(fruit|apple|apricot|date|fig|raisin|prune)|",
          "fruit\\s*(mix|mixture)"
        ),
        ignore_case = TRUE
      )
    ),
    # Half: mixtures with non-fruits; baby food mixtures
    fruit_half = str_detect(
      dl,
      regex(
        paste0(
          "fruit.*(mix|mixture).*(with|&).*(cereal|yogurt|granola|nuts|ice cream|dessert)\\b|",
          "baby\\s*food.*fruit.*(mix|mixture)"
        ),
        ignore_case = TRUE
      )
    ),
    # Exclude: juices; fruit desserts & fruit-flavored puddings
    fruit_excl = str_detect(
      dl,
      regex(
        paste0(
          "\\bjuice\\b|",
          "fruit.*(dessert|pie|cobbler|crisp|sherbet|sorbet|pudding|gelatin|gel dessert)|",
          "\\bgelatin\\b|\\bpudding\\b"
        ),
        ignore_case = TRUE
      )
    ),
    
    # =======================
    # SSB & 100% JUICE
    bev_core    = str_detect(dl, regex("soda|cola|soft\\s*drink|\\bpop\\b|lemonade|fruit\\s*(ade|drink|punch)|sports\\s*drink|energy\\s*drink|sweetened\\s*water|smoothie|frappuccino", TRUE)),
    bev_diet    = str_detect(dl, regex("\\bdiet\\b|sugar[- ]?free|unsweetened|zero\\b|low\\s*cal", TRUE)),
    bev_milkct  = str_detect(dl, regex("\\bmilk\\b|coffee|tea|cappuccino|latte|mocha|cocoa|hot\\s*chocolate", TRUE)),
    bev_reduced = str_detect(dl, regex("reduced\\s*sugar|less\\s*sugar|lower\\s*sugar|50%\\s*less\\s*sugar|\\blight\\b", TRUE)),
    
    ssb_full_flag = bev_core & !bev_diet & !bev_milkct & !bev_reduced,
    ssb_half_flag = bev_core & !bev_diet & !bev_milkct &  bev_reduced,
    
    is_juice       = str_detect(dl, regex("\\bjuice\\b", TRUE)),
    is_cocktail    = str_detect(dl, regex("cocktail|\\bdrink\\b|beverage|ade\\b|punch\\b|nectar", TRUE)),
    juice100_flag  = is_juice & !is_cocktail & !bev_diet & !bev_milkct,
    
    # =======================
    # RED / PROCESSED MEAT
    beef_any  = str_detect(dl, "\\bbeef\\b|\\bhamburger\\b"),
    beef_core = str_detect(dl, "beef\\s*(steak|rib|ribs|roast|ground|patty|patties|meatball|sirloin|ribeye|brisket)") |
      str_detect(dl, "\\bhamburger\\s*patty\\b"),
    beef_mix  = beef_any & str_detect(dl, "baby\\s*food|gravy|sauce|with\\s*(starch|potato|rice|pasta|noodle|vegetable)|sandwich|burger|frozen|tv\\s*dinner|meal|soup|broth|extract"),
    beef_excl = str_detect(dl, "beef\\s*bacon|frankfurter|wiener|hot\\s*dog|sausage|luncheon|lunch\\s*meat|bologna|salami|pepperoni|potted\\s*meat|meat\\s*spread|spam|pastrami|corned\\s*beef"),
    
    pork_any  = str_detect(dl, "\\bpork\\b|\\bham\\b|sparerib|spareribs"),
    pork_core = str_detect(dl, "pork\\s*(chop|steak|cutlet|roast|loin|tenderloin|sparerib|spareribs)|\\bham\\b"),
    pork_mix  = pork_any & str_detect(dl, "baby\\s*food|gravy|sauce|with\\s*(starch|potato|rice|pasta|noodle|vegetable)|sandwich|frozen|tv\\s*dinner|meal|soup|broth|extract"),
    pork_excl = str_detect(dl, "canadian\\s*bacon|\\bbacon\\b|salt\\s*pork|pork\\s*skin|frankfurter|wiener|hot\\s*dog|sausage\\b|luncheon|lunch\\s*meat|bologna|salami|pepperoni|potted\\s*meat|meat\\s*spread|spam"),
    
    beef_full = beef_core & !beef_mix & !beef_excl,
    beef_half = beef_mix  & !beef_excl,
    pork_full = pork_core & !pork_mix & !pork_excl,
    pork_half = pork_mix  & !pork_excl,
    
    proc_core = str_detect(dl, "canadian\\s*bacon|\\bbacon\\b|salt\\s*pork|pork\\s*skin|frankfurter|wiener|hot\\s*dog|sausage\\b|bratwurst|kielbasa|luncheon|lunch\\s*meat|cold\\s*cut|bologna|salami|pepperoni|potted\\s*meat|meat\\s*spread|liverwurst|pastrami|corned\\s*beef|spam|scrapple"),
    proc_half = str_detect(dl, "(sandwich|sub|hoagie|roll|bun).*(frankfurter|hot\\s*dog|wiener|luncheon|lunch\\s*meat|potted\\s*meat)"),
    
    # =======================
    # NUTS & LEGUMES (serving sizes per your spec)
    # ---- Nuts/Legumes/Seeds (50 g per serving)
    nls_full = str_detect(
      dl,
      regex(
        paste0(
          # beans/peas/lentils plain
          "\\b(black|pinto|kidney|navy|refried|garbanzo|chickpea|white|great northern|cannellini)\\b.*bean|",
          "\\bdried\\s*bean\\b|\\bdried\\s*beans\\b|\\bbaked\\s*bean\\b|\\bbaked\\s*beans\\b|",
          "\\blentil\\b|\\blentils\\b|\\bsplit\\s*pea\\b|\\bgreen\\s*pea(s)?\\b|\\bpea soup\\b(?!.*cream)|",
          # nuts & seeds
          "almond|walnut|pecan|peanut\\b(?!\\s*butter)|cashew|pistachio|hazelnut|macadamia|brazil\\s*nut|",
          "sunflower\\s*seed|pumpkin\\s*seed|sesame\\s*seed|chia\\s*seed|flax\\s*seed|hemp\\s*seed"
        ),
        TRUE
      )
    ),
    # Half: frozen plate meals, soups with legumes, legumes baby food
    nls_half = str_detect(
      dl,
      regex(
        "(frozen|tv\\s*dinner|plate\\s*meal).*(bean|lentil|pea)|bean.*soup|lentil.*soup|pea.*soup|baby\\s*food.*(bean|lentil|pea)",
        TRUE
      )
    ),
    # Exclude: soy-derived products, meat substitutes, nut butters, coconut beverages, carob products
    nls_excl = str_detect(
      dl,
      regex("soy\\b|tofu|tempeh|textured\\s*veg|meat\\s*substitute|veggie\\s*burger|nut\\s*butter|coconut\\s*(milk|beverage)|carob", TRUE)
    ),
    
    # ---- Nut butters (32 g per serving)
    nb_full = str_detect(dl, regex("(peanut|almond|cashew|hazelnut|sunflower)\\s*butter\\b", TRUE)),
    nb_half = str_detect(dl, regex("(peanut|nut)\\s*butter\\s*(and\\s*jelly|&\\s*jelly)|nut\\s*butter\\s*sandwich", TRUE)),
    nb_excl = str_detect(dl, regex("nut\\s*gravy|peanut\\s*sauce", TRUE)),
    
    # ---- Tofu / soybeans / soy nuts (50 g per serving)
    tofu_full = str_detect(dl, regex("tofu|soy\\s*bean\\b|soybeans\\b|soy\\s*nuts\\b", TRUE)),
    tofu_half = str_detect(dl, regex("soy\\s*yogurt|soy\\s*dessert|tofu\\s*soup|tofu.*(mix|dish)", TRUE)),
    
    # ---- Soy milk (50 g per serving; includes low-fat)
    soymilk_full = str_detect(dl, regex("\\bsoy\\s*(milk|drink|beverage)\\b", TRUE)),
    
    # =======================
    # WHOLE GRAINS (AHEI uses grams/day; here we do keyword flag; optional upgrade to 10:1 later)
    wg_kw = str_detect(
      dl,
      regex(
        paste0(
          "100%\\s*whole|\\bwhole\\s*(wheat|grain|rye|oat)|\\bwholegrain\\b|",
          "oatmeal|rolled\\s*oats|steel[- ]cut\\s*oats|bran\\b|bran\\s*flakes|",
          "brown\\s*rice\\b|popcorn\\b|bulgur|quinoa|farro|barley\\b|",
          "(whole|wholegrain).*\\b(bread|tortilla|pita|pasta|noodle|cracker|cereal)\\b"
        ),
        TRUE
      )
    )
    # If you later build a nutrient-based 10:1 table (by FOODCODE), add: wg_flag <- wg_kw | wg_ratio_flag
  ) %>%
  summarise(
    # ---------- Vegetables ----------
    veg_serv_full = sum(if_else(veg_full & !veg_excl, g/65,        0), na.rm = TRUE),
    veg_serv_half = sum(if_else(veg_half & !veg_excl, (g/65)*0.5,  0), na.rm = TRUE),
    veg_serv_iff  = veg_serv_full + veg_serv_half,
    veg_cup_eq_iff = veg_serv_iff * 0.5,
    
    # ---------- Fruit ----------
    fruit_serv_full = sum(if_else(fruit_full & !fruit_excl, g/65,        0), na.rm = TRUE),
    fruit_serv_half = sum(if_else(fruit_half & !fruit_excl, (g/65)*0.5,  0), na.rm = TRUE),
    fruit_serv_iff  = fruit_serv_full + fruit_serv_half,
    fruit_cup_eq_iff = fruit_serv_iff * 0.5,
    
    # ---------- SSB + 100% JUICE ----------
    ssb_serv_full = sum((fl_oz/8) * as.numeric(ssb_full_flag), na.rm = TRUE),
    ssb_serv_half = sum((fl_oz/8) * as.numeric(ssb_half_flag), na.rm = TRUE),
    ssb_serv      = ssb_serv_full + 0.5 * ssb_serv_half,
    fruit_juice_serv = sum((fl_oz/4) * as.numeric(juice100_flag), na.rm = TRUE),
    ssb_juice_serv   = ssb_serv + fruit_juice_serv,
    
    # ---------- RED + PROCESSED MEAT ----------
    beef_serv_fw = sum(if_else(beef_full, g/100,          0), na.rm = TRUE),
    beef_serv_hw = sum(if_else(beef_half, (g/100)*0.5,    0), na.rm = TRUE),
    pork_serv_fw = sum(if_else(pork_full, g/100,          0), na.rm = TRUE),
    pork_serv_hw = sum(if_else(pork_half, (g/100)*0.5,    0), na.rm = TRUE),
    red_serv     = beef_serv_fw + beef_serv_hw + pork_serv_fw + pork_serv_hw,
    
    proc_serv_fw = sum(if_else(proc_core, g/100,          0), na.rm = TRUE),
    proc_serv_hw = sum(if_else(proc_half, (g/100)*0.5,    0), na.rm = TRUE),
    redproc_serv_iff = red_serv + proc_serv_fw + proc_serv_hw,
    
    # ---------- NUTS & LEGUMES ----------
    # Nuts/Legumes/Seeds (exclude nls_excl; 50 g per serving)
    nls_serv_full = sum(if_else(nls_full & !nls_excl, g/50,        0), na.rm = TRUE),
    nls_serv_half = sum(if_else(nls_half & !nls_excl, (g/50)*0.5,  0), na.rm = TRUE),
    
    # Nut butters (32 g per serving; exclude gravy/sauce)
    nb_serv_full  = sum(if_else(nb_full & !nb_excl, g/32,        0), na.rm = TRUE),
    nb_serv_half  = sum(if_else(nb_half & !nb_excl, (g/32)*0.5,  0), na.rm = TRUE),
    
    # Tofu / soybeans / soy nuts (50 g per serving)
    tofu_serv_full = sum(if_else(tofu_full, g/50,        0), na.rm = TRUE),
    tofu_serv_half = sum(if_else(tofu_half, (g/50)*0.5,  0), na.rm = TRUE),
    
    # Soy milk (50 g per serving)
    soymilk_serv_full = sum(if_else(soymilk_full, g/50, 0), na.rm = TRUE),
    
    # Total nuts+legumes AHEI servings from IFF
    nuts_legumes_serv_iff = nls_serv_full + nls_serv_half +
      nb_serv_full  + nb_serv_half  +
      tofu_serv_full + tofu_serv_half +
      soymilk_serv_full,
    
    # ---------- WHOLE GRAINS (grams) ----------
    wholegr_g_iff = sum(if_else(wg_kw, g, 0), na.rm = TRUE),
    
    # ---------- Alcohol (QC only; AHEI scoring still uses alcohol_g from DR) ----------
    alc_noalco   = FALSE, # already filtered in flags; left for clarity
    beer_full    = str_detect(dl, "\\b(beer|lager|ale|ipa|stout|porter|pilsner)\\b") & !str_detect(dl,"non[- ]?alcoholic|near\\s*beer|root\\s*beer|ginger\\s*beer"),
    wine_full    = str_detect(dl, "\\b(wine|merlot|cabernet|pinot|chardonnay|riesling|zinfandel|sauvignon)\\b") & !str_detect(dl,"cooking|non[- ]?alcoholic"),
    liq_full     = str_detect(dl, "whisk(e)?y|bourbon|scotch|vodka|rum|gin|tequila|brandy|cognac|liqueur") & !str_detect(dl,"non[- ]?alcoholic"),
    alc_cocktail = str_detect(dl, "cocktail|mixed\\s*drink|margarita|mojito|martini|bloody\\s*mary|long\\s*island|cosmopolitan|daiquiri|pi[nn]a\\s*colada|sangria|spritzer|michelada|shandy"),
    
    alc_beer_serv_full = sum(if_else(beer_full,  g/340.2,   0), na.rm = TRUE),
    alc_wine_serv_full = sum(if_else(wine_full,  g/141.75,  0), na.rm = TRUE),
    alc_liq_serv_full  = sum(if_else(liq_full,   g/42.53,   0), na.rm = TRUE),
    alc_cocktail_half  = sum(if_else(alc_cocktail, (g/42.53)*0.5, 0), na.rm = TRUE),
    alc_drinks_serv    = alc_beer_serv_full + alc_wine_serv_full + alc_liq_serv_full + alc_cocktail_half,
    
    .by = c(SEQN, DAY, cycle)
  )

iff_servings$wholegr_g_iff


# 5) Build AHEI inputs (MPED + nutrients), Day 1 only — NO coalesce -------


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

dr_9904_mped$DRDTSODI
dr_std$sodium_mg


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


# MPED food groups -> working inputs (NO coalesce)

# Keep valid Day-1 weighted recalls , join nutrients, and require kcal present.

names(mped_9904_person)
mped_9904_person$F_OTHER
mped_9904_person$F_TOTAL

mped_9904_person$A_BEV

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
    
    # --- WHOLE FRUIT: F_CITMLB + F_OTHER (do NOT use F_TOTAL)
    # F_CITMLB  = if ("F_CITMLB"  %in% names(mped_9904_person))   F_CITMLB  else NA_real_,
    # F_OTHER   = if ("F_OTHER"   %in% names(mped_9904_person))   F_OTHER   else NA_real_,
    # fruit_cup_eq = ifelse(is.na(F_CITMLB) | is.na(F_OTHER),
    #                      NA_real_,
    #                      F_CITMLB + F_OTHER),
    
    
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
    # ssb_serv         = NA_real_,
    # fruit_juice_serv = NA_real_,     # set later from DRXIFF if available
    # ssb_juice_serv   = ifelse(is.na(ssb_serv) & is.na(fruit_juice_serv), NA_real_,
    #                           ifelse(is.na(ssb_serv), 0, ssb_serv) + ifelse(is.na(fruit_juice_serv), 0, fruit_juice_serv)),
    
    A_BEV    = if ("A_BEV"   %in% names(mped_9904_person))   A_BEV   else NA_real_,
    ssb_serv = A_BEV,
    
    # --- Nutrients (carried through)
    energy_kcal, sodium_mg, alcohol_g, pufa_g, epa_g, dha_g
  )




ahei_input_mped <- mped_9904_person %>%
  filter(!is.na(WTDRD1) & WTDRD1 > 0) %>%
  left_join(nutrients_9904, by = "SEQN") %>%
  filter(!is.na(energy_kcal) & energy_kcal > 0) %>%
  left_join(
    iff_bevs %>% select(SEQN, DAY, ssb_serv, fruit_juice_serv, ssb_juice_serv),
    by = c("SEQN")  # day is 1 for your build; if you carry DAY in mped_9904_person, use by=c("SEQN","DAY")
  ) %>%
  transmute(
    SEQN, WTDRD1, RIAGENDR,
    V_TOTAL   = if ("V_TOTAL"  %in% names(mped_9904_person)) V_TOTAL  else NA_real_,
    V_POTATO  = if ("V_POTATO" %in% names(mped_9904_person)) V_POTATO else NA_real_,
    veg_cup_eq = sub_or_keep(V_TOTAL, V_POTATO),
    F_TOTAL     = if ("F_TOTAL" %in% names(mped_9904_person)) F_TOTAL else NA_real_,
    fruit_cup_eq= F_TOTAL,
    wholegr_oz_eq = if ("G_WHL" %in% names(mped_9904_person)) G_WHL else NA_real_,
    wholegr_g     = ifelse(is.na(wholegr_oz_eq), NA_real_, wholegr_oz_eq * 28.3495),
    nuts_oz_eq        = if ("M_NUTSD" %in% names(mped_9904_person)) M_NUTSD else NA_real_,
    legumes_cup_eq    = if ("LEGUMES" %in% names(mped_9904_person)) LEGUMES else NA_real_,
    nuts_legumes_serv = ifelse(is.na(nuts_oz_eq) | is.na(legumes_cup_eq),
                               ifelse(is.na(nuts_oz_eq) & is.na(legumes_cup_eq), NA_real_,
                                      (ifelse(is.na(nuts_oz_eq), 0, nuts_oz_eq) + ifelse(is.na(legumes_cup_eq), 0, legumes_cup_eq) / 0.5)),
                               nuts_oz_eq + (legumes_cup_eq / 0.5)),
    red_oz_eq   = if ("M_MEAT"  %in% names(mped_9904_person)) M_MEAT  else NA_real_,
    proc_oz_eq  = if ("M_FRANK" %in% names(mped_9904_person)) M_FRANK else NA_real_,
    redproc_oz_eq = ifelse(is.na(red_oz_eq) & is.na(proc_oz_eq), NA_real_,
                           ifelse(is.na(red_oz_eq), 0, red_oz_eq) + ifelse(is.na(proc_oz_eq), 0, proc_oz_eq)),
    redproc_serv  = ifelse(is.na(redproc_oz_eq), NA_real_, redproc_oz_eq / 3.527),
    
    # >>> use SSB + 100% juice from IFF
    ssb_serv         = ssb_serv,
    fruit_juice_serv = fruit_juice_serv,
    ssb_juice_serv   = ssb_juice_serv,
    
    energy_kcal, sodium_mg, alcohol_g, pufa_g, epa_g, dha_g
  )


# Sanity messages
cat("AHEI input rows after valid Day-1 & non-missing kcal: ", nrow(ahei_input_mped), "\n")

# mean(ahei_input_mped$veg_cup_eq, na.rm = TRUE) # optional peek


# 6) AHEI component scoring — keep NA as NA (no NA->0)-------


# ---- helpers ----
lin_pos <- function(x, min0, max10) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (x - min0)/(max10 - min0))) * 10)
lin_rev <- function(x, min10, max0) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (max0 - x)/(max0 - min10))) * 10)

# If you have DEMO (RIAGENDR: 1=Male, 2=Female) merged in ahei_input_mped, this will be used.
# If RIAGENDR is missing, whole-grain score becomes NA to enforce complete-case exclusion.
sex_vec <- if ("RIAGENDR" %in% names(ahei_input_mped)) ahei_input_mped$RIAGENDR else rep(NA_real_, nrow(ahei_input_mped))

# ---- 6.1) Vegetables (servings/day; 0 -> 0 pts; 5 -> 10 pts) ----


# Choose reference energy for density scaling (AHEI work commonly uses 2000 kcal)
E_REF <- 2000

# helper: x per 1000 kcal, NA if x or kcal is NA
per_1000 <- function(x, kcal) ifelse(is.na(x) | is.na(kcal), NA_real_, x / (kcal / 1000))

# density cutpoints scaled from per-day AHEI targets
VEG_CUTOFF_PER1000   <- 5 / (E_REF / 1000)   # 5/day at 2000 kcal -> 2.5/1000
FRUIT_CUTOFF_PER1000 <- 4 / (E_REF / 1000)   # 4/day at 2000 kcal -> 2.0/1000

# Vegetables (past density spec): max at ≥5 cups/1000 kcal
# 1 serving = 0.5 cup; ensure veg_cup_eq is non-starchy veg
ahei_veg <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    veg_serv = ifelse(is.na(veg_cup_eq), NA_real_, veg_cup_eq / 0.5),
    ahei_veg = lin_pos(veg_serv, 0, 5)
  )

summary(ahei_veg)

## adjusted Vegetables (MAX at ≥ 2.5 cups per 1000 kcal)
ahei_veg <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    veg_per_1000 = per_1000(veg_cup_eq, energy_kcal),
    # optional sanity filter for outliers
    veg_per_1000 = ifelse(!is.na(veg_per_1000) & veg_per_1000 > 10, NA_real_, veg_per_1000),
    ahei_veg = ifelse(
      is.na(veg_per_1000), NA_real_,
      pmin(veg_per_1000 / VEG_CUTOFF_PER1000, 1) * 10
    )
  )

summary(ahei_veg)

# ---- 6.2) Fruit (servings/day; 0 -> 0 pts; 4 -> 10 pts) ----
# 1 serving = 0.5 cup; exclude juice upstream if you have it
ahei_fruit <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    fruit_serv = ifelse(is.na(fruit_cup_eq), NA_real_, fruit_cup_eq / 0.5),
    ahei_fruit = lin_pos(fruit_serv, 0, 4)
  )

summary(ahei_fruit)


## adjusted Fruit (MAX at ≥ 2.0 cups per 1000 kcal) — unchanged target, made explicit-----
# ahei_fruit <- ahei_input_mped %>%
#  dplyr::transmute(
#    SEQN,
#    fruit_per_1000 = per_1000(fruit_cup_eq, energy_kcal),
#    fruit_per_1000 = ifelse(!is.na(fruit_per_1000) & fruit_per_1000 > 10, NA_real_, fruit_per_1000),
#    ahei_fruit = ifelse(
#      is.na(fruit_per_1000), NA_real_,
#      pmin(fruit_per_1000 / FRUIT_CUTOFF_PER1000, 1) * 10
#    )
#  )

summary(ahei_fruit)


# ---- 6.3) Whole grains (g/day; women 75g -> 10, men 90g -> 10) ----
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

summary(ahei_grain)

# ---- 6.4) SSB + 100% fruit juice (servings/day; 0 -> 10 pts; ≥1 -> 0 pts) ----
ahei_ssb <- ahei_input_mped %>%
  transmute(SEQN, ahei_ssb = lin_rev(ssb_juice_serv, 0, 1))

summary(ahei_ssb)

# ---- 5.5) Nuts & legumes (servings/day; 0 -> 0 pts; 1 -> 10 pts) ----
ahei_nutslegumes <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    ahei_nutslegumes = lin_pos(nuts_legumes_serv, 0, 1)
  )

summary(ahei_nutslegumes)

# ---- 5.6) Red + processed meat (servings/day; 0 -> 10 pts; ≥1.5 -> 0 pts) ----
ahei_meat <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    ahei_redprocmeat = lin_rev(redproc_serv, 0, 1.5)
  )

summary(ahei_meat)

# ---- 5.7) Long-chain n-3 (EPA + DHA), mg/day; 0 -> 0 pts; 250 mg -> 10 pts ----
ahei_longn3 <- ahei_input_mped %>%
  dplyr::transmute(
    SEQN,
    lc_n3_mg = ifelse(is.na(epa_g) & is.na(dha_g),
                      NA_real_,
                      (ifelse(is.na(epa_g), 0, epa_g) + ifelse(is.na(dha_g), 0, dha_g)) * 1000),
    ahei_longn3 = lin_pos(lc_n3_mg, 0, 250)
  )

summary(ahei_longn3)

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

summary(ahei_pufa)

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

summary(ahei_alcohol)

# ---- 5.10) Sodium —----
# standard AHEI-2010 sodium scoring (10 at ≤1100 mg/1000 kcal, 0 at ≥2000 mg/1000 kcal, 
# linear in between) and is what Wang et al. use

ahei_sodium <- ahei_input_mped %>%
  dplyr::select(SEQN, sodium_mg, energy_kcal) %>%
  # exclude if either piece is missing or energy is 0/neg
  dplyr::filter(!is.na(sodium_mg), !is.na(energy_kcal), energy_kcal > 0) %>%
  dplyr::mutate(
    sodium_per_1000kcal = sodium_mg / (energy_kcal / 1000),
    ahei_sodium = dplyr::case_when(
      sodium_per_1000kcal <= 1100 ~ 10,
      sodium_per_1000kcal >= 2000 ~ 0,
      TRUE ~ (2000 - sodium_per_1000kcal) / (2000 - 1100) * 10
    )
  ) %>%
  dplyr::select(SEQN, ahei_sodium)

summary(ahei_sodium)

# =========================
# Combine components -> exclude any missing component -> total
# =========================
ahei_all <- list(
  ahei_veg, ahei_fruit, ahei_grain, ahei_nutslegumes,
  ahei_meat, ahei_longn3, ahei_pufa, ahei_sodium, ahei_alcohol #,  ahei_ssb
) %>% purrr::reduce(dplyr::left_join, by = "SEQN")

# Exclude participants with ANY missing component (Patel/Wang approach)
ahei_complete <- ahei_all %>%
  dplyr::filter(dplyr::if_all(dplyr::starts_with("ahei_"), ~ !is.na(.))) %>%
  dplyr::mutate(
    ahei_total = rowSums(dplyr::across(dplyr::starts_with("ahei_")), na.rm = FALSE)
  )

message("Excluded ", nrow(ahei_all) - nrow(ahei_complete), " participants due to ≥1 missing component.")

summary(ahei_all)





# 7) Save + QC -----------------------------------------------------------------

readr::write_csv(ahei_complete, file.path(dir$output, "ahei_1999_2004_day1.csv"))

summary(ahei_complete$ahei_total)

# Quick histogram
ggplot2::ggplot(ahei_complete, ggplot2::aes(x = ahei_total)) +
  ggplot2::geom_histogram(binwidth = 5) +
  ggplot2::labs(title = "AHEI (1999–2004, Day 1)", x = "AHEI total", y = "Count") +
  ggplot2::theme_minimal()


