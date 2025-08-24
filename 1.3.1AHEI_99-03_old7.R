# ================================================================
# AHEI 1999–2004 — MPED IFF (item-level) with mixture 0.5,
# juice exclusions, beverage block, and pyr_tot fallback.
# Uses your original reading flow and directory layout.
# ================================================================

# 1) Paths, packages, helpers ------------------------------------

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
invisible(lapply(dir[c("data","output","code","fped","nhanes")], dir.create,
                 showWarnings = FALSE, recursive = TRUE))

# R session defaults
options(stringsAsFactors = FALSE, scipen = 999,
        repos = c(CRAN = "https://cloud.r-project.org"))
set.seed(148)

# Packages
pkgs <- c("dplyr","haven","foreign","survey","purrr","ggplot2",
          "readr","stringr","tidyr","rlang","janitor","tibble")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, quiet = TRUE)
suppressPackageStartupMessages(invisible(lapply(pkgs, require, character.only = TRUE)))

cat("Project root:", dir$root, "\nData dir:", dir$data, "\nFPED/MPED dir:", dir$fped, "\n\n")

# ---- small helpers
`%||%` <- function(a,b) if (!is.null(a)) a else b
wmean <- function(x, w) { i <- !is.na(x) & !is.na(w); if (!any(i)) return(NA_real_); sum(x[i]*w[i]) / sum(w[i]) }
mean_preserve_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

safe_read_sas <- function(path) {
  if (!file.exists(path)) { message("Missing file: ", path); return(NULL) }
  tryCatch(haven::read_sas(path), error = function(e){ message("read_sas failed: ", path); NULL })
}
safe_read_xpt <- function(path) {
  if (!file.exists(path)) { message("Missing file: ", path); return(NULL) }
  tryCatch(haven::read_xpt(path), error = function(e){ message("read_xpt failed: ", path); NULL })
}
pick_name <- function(nms, candidates) { hits <- candidates[candidates %in% nms]; if (length(hits)) hits[1] else NA_character_ }
safe_vec <- function(df, candidates, default = NA_real_) {
  nm <- pick_name(names(df), candidates)
  if (is.na(nm)) rep(default, nrow(df)) else df[[nm]]
}
infer_day_from_name <- function(path) {
  b <- tolower(basename(path))
  if (grepl("dr1|day1|_d1\\b", b)) return(1)
  if (grepl("dr2|day2|_d2\\b", b)) return(2)
  1
}
report_component <- function(df, score_col, title = NULL) {
  tmp <- df %>%
    dplyr::select(SEQN, dplyr::all_of(score_col)) %>%
    dplyr::left_join(nutrients_9904 %>% dplyr::select(SEQN, WTDRD1), by = "SEQN") %>%
    dplyr::mutate(WTDRD1 = as.numeric(WTDRD1))
  sc  <- tmp[[score_col]]
  wt  <- tmp$WTDRD1
  wt6 <- wt / 3
  if (!is.null(title)) cat("\n---", title, "---\n")
  print(summary(sc))
  cat("Mean (unweighted):", mean(sc, na.rm = TRUE), "\n")
  cat("Mean (WTDRD1):    ", wmean(sc, wt), "\n")
  cat("Mean (WTDRD1/3):  ", wmean(sc, wt6), "\n")
}

# 2) MPED MyPyramid totals (1999–2004) ---------------------------

# Paths
pyr_tot_9902_path    <- file.path(dir$fped, "pyr_tot.sas7bdat")  # MPED v1.0 totals (1999–2002)
pyr_tot_d1_0304_path <- file.path(dir$fped, "MyPyrEquivDB2", "data", "intakes", "pyr_tot_d1.sas7bdat")
pyr_tot_d2_0304_path <- file.path(dir$fped, "MyPyrEquivDB2", "data", "intakes", "pyr_tot_d2.sas7bdat")

# Read
pyr_tot_9902    <- safe_read_sas(pyr_tot_9902_path)
pyr_tot_d1_0304 <- safe_read_sas(pyr_tot_d1_0304_path)
pyr_tot_d2_0304 <- safe_read_sas(pyr_tot_d2_0304_path)

# Tag DAY for 2003–2004
if (!is.null(pyr_tot_d1_0304)) pyr_tot_d1_0304 <- dplyr::mutate(pyr_tot_d1_0304, DAY = 1)
if (!is.null(pyr_tot_d2_0304)) pyr_tot_d2_0304 <- dplyr::mutate(pyr_tot_d2_0304, DAY = 2)

# Harmonize DAY for 1999–2002 and force Day 1 for join
pyr_tot_9902_d1 <- if (!is.null(pyr_tot_9902)) {
  tmp <- pyr_tot_9902
  nm  <- names(tmp)
  if ("DAYNO" %in% nm && !("DAY" %in% nm)) names(tmp)[nm == "DAYNO"] <- "DAY"
  if (!("DAY" %in% names(tmp))) tmp <- dplyr::mutate(tmp, DAY = 1)
  dplyr::filter(tmp, DAY == 1)
} else NULL

# Keep Day 1 for 2003–2004
pyr_tot_0304_d1 <- if (!is.null(pyr_tot_d1_0304)) dplyr::mutate(pyr_tot_d1_0304, DAY = 1) else NULL

cat(
  "Loaded MPED files:\n",
  "- 1999–2002 totals: ", ifelse(is.null(pyr_tot_9902), "NO", "YES"), "\n",
  "- 2003–2004 Day1:   ", ifelse(is.null(pyr_tot_d1_0304), "NO", "YES"), "\n",
  "- 2003–2004 Day2:   ", ifelse(is.null(pyr_tot_d2_0304), "NO", "YES"), "\n\n",
  sep = ""
)

pyr_tot_0304_d1$SEQN
# 3) NHANES Day-1 totals (1999–2004) + nutrients -----------------

dr_paths <- tibble::tribble(
  ~label,       ~dr1,                                        ~dr2,
  "1999-2000",  file.path(dir$nhanes, "DRXTOT.XPT"),          NA_character_,
  "2001-2002",  file.path(dir$nhanes, "DRXTOT_B.XPT"),        NA_character_,
  "2003-2004",  file.path(dir$nhanes, "DR1TOT_C.XPT"),        NA_character_
)

read_dr <- function(dr_path, day_num){
  if (is.na(dr_path)) return(NULL)
  df <- safe_read_xpt(dr_path)
  if (is.null(df)) return(NULL)
  dplyr::mutate(df, DAY = day_num)
}

dr_9900 <- read_dr(dr_paths$dr1[dr_paths$label=="1999-2000"], 1)
dr_0102 <- read_dr(dr_paths$dr1[dr_paths$label=="2001-2002"], 1)
dr_0304 <- read_dr(dr_paths$dr1[dr_paths$label=="2003-2004"], 1)
dr_all_d1 <- dplyr::bind_rows(dr_9900, dr_0102, dr_0304)
if (nrow(dr_all_d1) == 0) stop("No DR Day-1 rows found.")

# DEMO (gender for alcohol/whole grains sex-specific)
demo_paths <- tibble::tribble(
  ~label,       ~demo,
  "1999-2000",  file.path(dir$nhanes, "DEMO.XPT"),
  "2001-2002",  file.path(dir$nhanes, "DEMO_B.XPT"),
  "2003-2004",  file.path(dir$nhanes, "DEMO_C.XPT")
)
demo_all <- demo_paths$demo %>%
  purrr::map(safe_read_xpt) %>% dplyr::bind_rows() %>% dplyr::select(SEQN, RIAGENDR)

demo_all$SEQN
# Standardize nutrient columns (Day 1)
pick_first <- function(nms, cand) { hit <- cand[cand %in% nms]; if (length(hit)) hit[1] else NA_character_ }
en_col  <- pick_first(names(dr_all_d1), c("DRXTKCAL","DR1TKCAL","DR1IKCAL"))
na_col  <- pick_first(names(dr_all_d1), c("DRDTSODI","DR1TSODI","DR1ISODI"))
alc_col <- pick_first(names(dr_all_d1), c("DRXTALCO","DR1TALCO","DR1IALCO"))
puf_col <- pick_first(names(dr_all_d1), c("DRXTPFAT","DR1TPFAT","DR1IPFAT"))
epa_col <- pick_first(names(dr_all_d1), c("DRXTP205","DR1TP205","DR1IP205"))
dha_col <- pick_first(names(dr_all_d1), c("DRXTP226","DR1TP226","DR1IP226"))
dpa_col <- pick_first(names(dr_all_d1), c("DRXTP225","DR1TP225","DR1IP225"))
w_col   <- pick_first(names(dr_all_d1), c("WTDRD1"))

dr_std <- dr_all_d1 %>%
  dplyr::transmute(
    SEQN,
    DAY         = 1,
    WTDRD1      = if (!is.na(w_col))   .data[[w_col]]   else NA_real_,
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
    WTDRD1      = mean_preserve_na(WTDRD1),
    energy_kcal = mean_preserve_na(energy_kcal),
    sodium_mg   = mean_preserve_na(sodium_mg),
    alcohol_g   = mean_preserve_na(alcohol_g),
    pufa_g      = mean_preserve_na(pufa_g),
    epa_g       = mean_preserve_na(epa_g),
    dha_g       = mean_preserve_na(dha_g),
    dpa_g       = mean_preserve_na(dpa_g),
    .groups = "drop"
  ) %>%
  dplyr::left_join(demo_all, by = "SEQN")

nutrients_9904$SEQN
# 3b) Join DR by cycle to MPED totals (Day 1) ---------------------

join_mped <- function(dr_df, mped_df) {
  if (is.null(dr_df) || !nrow(dr_df)) return(tibble())
  if (is.null(mped_df) || !nrow(mped_df)) return(dr_df)  # keep DR if MPED missing
  if ("DAY" %in% names(mped_df)) dplyr::left_join(dr_df, mped_df, by = c("SEQN","DAY"))
  else                           dplyr::left_join(dr_df, mped_df, by = "SEQN")
}


dr_9902_mped <- join_mped(dplyr::bind_rows(dr_9900, dr_0102), pyr_tot_9902_d1)
dr_0304_mped <- join_mped(dr_0304,                              pyr_tot_0304_d1)

dr_9904_mped <- dplyr::bind_rows(
  dr_9902_mped %||% dplyr::tibble(),
  dr_0304_mped %||% dplyr::tibble()
)


cat("DR Day1 rows:", nrow(dplyr::bind_rows(dr_9900, dr_0102, dr_0304)),
    "| after MPED join:", nrow(dr_9904_mped), "x", ncol(dr_9904_mped), "\n")

mped_9904_person <- dr_9904_mped %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

cat("Persons in MPED 1999–2004 (Day 1 only): ", nrow(mped_9904_person), "\n")

mped_9904_person$SEQN
# 4) Individual foods (IFF) + FNDDS DESC (1999–2004) --------------

iff_files <- tibble::tribble(
  ~cycle,        ~file,
  "1999-2000",   file.path(dir$nhanes, "DRXIFF.xpt"),
  "2001-2002",   file.path(dir$nhanes, "DRXIFF_B.xpt"),
  "2003-2004",   file.path(dir$nhanes, "DR1IFF_C.xpt")
)

read_iff_min <- function(path, cycle_label) {
  df <- safe_read_xpt(path); if (is.null(df)) return(tibble())
  nms <- names(df)
  seqn_nm  <- pick_name(nms, c("SEQN"))
  iline_nm <- pick_name(nms, c("DR1ILINE","DRXILINE","ILINE","LINENUM"))
  fcode_nm <- pick_name(nms, c("DR1IFDCD","DRDIFDCD","DRXIFDCD","FOODCODE","IFCODE"))
  grams_nm <- pick_name(nms, c("DR1IGRMS","DRXIGRMS","GRAMS"))
  wweia_nm <- pick_name(nms, c("DR1IWWEIA","DRXIWWEIA","WWEIACAT","WWEIACODE"))
  desc_nm  <- pick_name(nms, c("DR1IFDCDTX","DRXIFDCDTX","DR1I_FDNAME","FDNAM","FD_NAME"))
  n <- nrow(df)
  tibble(
    SEQN      = if (!is.na(seqn_nm))   as.numeric(df[[seqn_nm]])  else rep(NA_real_, n),
    ILINE     = if (!is.na(iline_nm))  as.numeric(df[[iline_nm]]) else rep(NA_real_, n),
    FOODCODE  = if (!is.na(fcode_nm))  as.numeric(df[[fcode_nm]]) else rep(NA_real_, n),
    GRAMS     = if (!is.na(grams_nm))  as.numeric(df[[grams_nm]]) else rep(NA_real_, n),
    WWEIA_CAT = if (!is.na(wweia_nm))  as.numeric(df[[wweia_nm]]) else rep(NA_real_, n),
    DESC      = if (!is.na(desc_nm))   as.character(df[[desc_nm]]) else rep(NA_character_, n),
    DAY       = 1,
    cycle     = cycle_label
  )
}

iff_all <- purrr::map2_dfr(iff_files$file, iff_files$cycle, read_iff_min) %>%
  dplyr::mutate(
    FOODCODE = as.numeric(FOODCODE),
    GRAMS    = as.numeric(GRAMS),
    DESC     = dplyr::coalesce(DESC, ""),
    dl       = stringr::str_to_lower(DESC)
  )
cat("IFF rows stacked (1999–2004): ", nrow(iff_all), "\n")

iff_all$SEQN
# ---- FNDDS ASCII lookups
read_fndds_lookup <- function(ascii_dir, main_file="MainFoodDesc.txt", add_file="AddFoodDesc.txt") {
  main_path <- file.path(ascii_dir, main_file)
  add_path  <- file.path(ascii_dir, add_file)
  if (!file.exists(main_path)) return(tibble::tibble(FOODCODE = numeric(), DESC = character()))
  main <- readr::read_delim(
    main_path, delim="^", quote="~",
    col_names=c("food_code","start_date","end_date","main_food_description",
                "main_food_description_upper","_extra"),
    escape_double=FALSE, trim_ws=TRUE, show_col_types=FALSE
  ) %>% dplyr::select(-`_extra`) %>%
    dplyr::transmute(FOODCODE = suppressWarnings(as.numeric(food_code)),
                     DESC_main = main_food_description)
  add <- if (file.exists(add_path)) {
    readr::read_delim(
      add_path, delim="^", quote="~",
      col_names=c("food_code","start_date","end_date","add_food_desc",
                  "add_food_desc_upper","_extra"),
      escape_double=FALSE, trim_ws=TRUE, show_col_types=FALSE
    ) %>% dplyr::select(-`_extra`) %>%
      dplyr::transmute(FOODCODE = suppressWarnings(as.numeric(food_code)),
                       DESC_add = add_food_desc) %>%
      dplyr::group_by(FOODCODE) %>%
      dplyr::summarise(DESC_add = paste0(unique(stats::na.omit(DESC_add)), collapse="; "), .groups="drop")
  } else tibble::tibble(FOODCODE = numeric(), DESC_add = character())
  main %>% dplyr::left_join(add, by="FOODCODE") %>%
    dplyr::mutate(
      DESC = dplyr::coalesce(
        dplyr::if_else(!is.na(DESC_add) &
                         !stringr::str_detect(DESC_main, stringr::fixed(DESC_add, ignore_case = TRUE)),
                       paste(DESC_main, "-", DESC_add), DESC_main),
        DESC_main)) %>%
    dplyr::transmute(FOODCODE, DESC)
}

ascii_0102 <- "/Users/dengshuyue/Desktop/SDOH/analysis/data/fndds/FNDDS1_ASCII_unpacked/fndds/ascii"
ascii_0304 <- "/Users/dengshuyue/Desktop/SDOH/analysis/data/fndds/FNDDS2/ascii"

fndds_lookup_0102 <- read_fndds_lookup(ascii_0102) %>% dplyr::mutate(cycle_src = "2001-2002")
fndds_lookup_0304 <- read_fndds_lookup(ascii_0304, main_file="mainfooddesc.txt", add_file="addfooddesc.txt") %>%
  dplyr::mutate(cycle_src = "2003-2004")

fndds_lookup <- dplyr::bind_rows(fndds_lookup_0102, fndds_lookup_0304) %>%
  dplyr::distinct(FOODCODE, .keep_all = TRUE) %>%
  dplyr::select(FOODCODE, DESC)

iff_joined <- iff_all %>%
  dplyr::select(-dplyr::any_of(c("DESC","desc_lc","cycle_src"))) %>%
  dplyr::mutate(FOODCODE = suppressWarnings(as.numeric(FOODCODE))) %>%
  dplyr::left_join(fndds_lookup %>% dplyr::select(FOODCODE, DESC), by = "FOODCODE") %>%
  dplyr::mutate(DESC = dplyr::coalesce(DESC, ""),
                dl   = stringr::str_to_lower(DESC))


##!!!!!!-----

mped_iff_joined <- iff_joined  # alias used later

mped_iff_joined$SEQN

# 5) MPED IFF-level equivalents (item level) ----------------------

# ---------- NHANES IFF (correct: dir$nhanes; DRXIFF*.xpt) ----------
iff_dir <- dir$nhanes

iff_files <- tibble::tribble(
  ~cycle,        ~file,
  "1999-2000",   file.path(iff_dir, "DRXIFF.xpt"),
  "2001-2002",   file.path(iff_dir, "DRXIFF_B.xpt"),
  "2003-2004",   file.path(iff_dir, "DR1IFF_C.xpt")
)

iff_all <- purrr::map2_dfr(iff_files$file, iff_files$cycle, read_iff_min)

# (keep your FNDDS lookup code as-is, then:)
iff_joined <- iff_all %>%
  dplyr::select(-dplyr::any_of(c("DESC", "desc_lc"))) %>%
  dplyr::mutate(FOODCODE = as.numeric(FOODCODE)) %>%
  dplyr::left_join(fndds_lookup %>% dplyr::select(FOODCODE, DESC), by = "FOODCODE") %>%
  dplyr::mutate(
    DESC    = dplyr::coalesce(DESC, ""),
    dl      = stringr::str_to_lower(DESC)
  )

iff_joined$SEQN

# ---------- MPED per-item equivalents (pyr_iff*) ----------

#### HERE IS THE ISSUE !!!!! fix in next version ------
# These are NOT the NHANES IFF; they usually live under dir$fped.
pyr_iff_files <- list.files(
  dir$fped, pattern = "(?i)^pyr_iff.*\\.(sas7bdat|xpt)$",
  recursive = TRUE, full.names = TRUE
)

if (!length(pyr_iff_files)) {
  message("No MPED per-item files (pyr_iff*) found under: ", dir$fped,
          " — item-level equivalents will be unavailable; will use totals fallback only.")
  pyr_iff_all <- tibble::tibble()
} else {
  read_pyr_iff <- function(path){
    df <- safe_read_sas(path)
    if (is.null(df)) return(tibble())
    day_guess <- infer_day_from_name(path)
    tibble::tibble(
      SEQN     = as.numeric(safe_vec(df, c("SEQN"), NA_real_)),
      DAY      = as.numeric(safe_vec(df, c("DAY","DAYNO","DAYN"), day_guess)),
      ILINE    = as.numeric(safe_vec(df, c("ILINE","DR1ILINE","DRXILINE","LINENUM"), NA_real_)),
      FOODCODE = as.numeric(safe_vec(df, c("FOODCODE","DR1IFDCD","DRDIFDCD","IFCODE"), NA_real_)),
      V_TOTAL  = as.numeric(safe_vec(df, c("V_TOTAL","V_TOTAL_IFF","V_TOT"), NA_real_)),
      V_POTATO = as.numeric(safe_vec(df, c("V_POTATO","V_POTAT"), 0)),
      F_TOTAL  = as.numeric(safe_vec(df, c("F_TOTAL","F_TOTAL_IFF","F_TOT"), NA_real_)),
      G_WHL    = as.numeric(safe_vec(df, c("G_WHL","G_WHOLE","G_WG"), NA_real_)),   # oz-eq/item
      LEGUMES  = as.numeric(safe_vec(df, c("LEGUMES"), NA_real_)),                  # cup-eq/item
      M_NUTSD  = as.numeric(safe_vec(df, c("M_NUTSD","M_NUTS"), NA_real_)),         # oz-eq/item
      M_MEAT   = as.numeric(safe_vec(df, c("M_MEAT","MEAT_RED","M_RED"), NA_real_)),# oz-eq/item
      M_FRANK  = as.numeric(safe_vec(df, c("M_FRANK","MEAT_PROC","M_PROC"), NA_real_))
    ) %>% dplyr::filter(!is.na(SEQN))
  }
  pyr_iff_all <- purrr::map_dfr(pyr_iff_files, read_pyr_iff)
}

pyr_iff_all$SEQN

# ---------- Join choice ----------
have_il_iff  <- "ILINE" %in% names(iff_joined)   && any(!is.na(iff_joined$ILINE))
have_il_mped <- "ILINE" %in% names(pyr_iff_all) && any(!is.na(pyr_iff_all$ILINE))
if (nrow(pyr_iff_all)) {
  # Preferred: MPED per-item rows, add IFF DESC/GRAMS for mixture flags
  if (have_il_iff && have_il_mped) {
    mped_iff_joined <- pyr_iff_all %>%
      dplyr::left_join(
        iff_joined %>% dplyr::select(SEQN, DAY, ILINE, FOODCODE, GRAMS, DESC, dl),
        by = c("SEQN","DAY","ILINE","FOODCODE")
      )
  } else {
    message("ILINE unavailable; joining on (SEQN, DAY, FOODCODE).")
    mped_iff_joined <- pyr_iff_all %>%
      dplyr::left_join(
        iff_joined %>% dplyr::select(SEQN, DAY, FOODCODE, GRAMS, DESC, dl),
        by = c("SEQN","DAY","FOODCODE")
      )
  }
} else {
  # No per-item MPED: keep IFF for beverage flags and mixture tags; item-level equivalents will be NA
  mped_iff_joined <- iff_joined
}

mped_iff_joined$SEQN

# 5b) Mixture & beverage flags; item-level totals -----------------

# mixture flag
mixture_flag <- function(desc) {
  x <- stringr::str_to_lower(desc %||% "")
  stringr::str_detect(
    x,
    stringr::regex(
      paste0(
        "soup|stew|casserole|chili|",
        "pizza|calzone|stromboli|",
        "(sandwich|wrap|pita|sub|hoagie|burger)(?!.*veggie)|",
        "burrito|taco|quesadilla|enchilada|nacho|fajita|empanada|",
        "fried\\s*rice|stir[- ]?fry|lo\\s*mein|chow\\s*mein|",
        "pasta\\s*with|noodles?\\s*with|rice\\s*with|",
        "mixed\\s*dish|tv\\s*dinner|frozen\\s*meal|plate\\s*meal"
      ),
      ignore_case = TRUE
    )
  )
}

# beverage flags
vegjuice_pat       <- "(?i)(\\bv8\\b|vegetable\\s*juice|veg\\s*juice|tomato\\s*juice|carrot\\s*juice|beet\\s*juice|clamato|low[- ]sodium\\s*vegetable\\s*juice)"
is_juice_pat       <- "\\bjuice\\b"
is_cocktail_pat    <- "cocktail|\\bdrink\\b|beverage|ade\\b|punch\\b|nectar"
milk_coffee_pat    <- "\\bmilk\\b|coffee|tea|cappuccino|latte|mocha|cocoa|hot\\s*chocolate"
juice_pack_ctx_pat <- "juice\\s*pack(ed)?|packed\\s*in\\s*(its\\s*|own\\s*)?juice|in\\s*(its\\s*|own\\s*)?juice"

mped_iff_w <- mped_iff_joined %>%
  dplyr::mutate(
    dl     = stringr::str_to_lower(DESC %||% ""),
    is_mix = mixture_flag(DESC),
    w_mix  = dplyr::if_else(is_mix, 0.5, 1.0),
    
    vegjuice_flag = stringr::str_detect(dl, vegjuice_pat),
    is_juice      = stringr::str_detect(dl, stringr::regex(is_juice_pat, TRUE)),
    is_cocktail   = stringr::str_detect(dl, stringr::regex(is_cocktail_pat, TRUE)),
    bev_milkct    = stringr::str_detect(dl, stringr::regex(milk_coffee_pat, TRUE)),
    juice_pack_context = stringr::str_detect(dl, stringr::regex(juice_pack_ctx_pat, TRUE)),
    juice100_flag = is_juice & !juice_pack_context & !is_cocktail & !bev_milkct
  )

# item-level contributions (Day 1 only)
mped_item <- mped_iff_w %>%
  dplyr::filter(is.na(DAY) | DAY == 1) %>%
  dplyr::transmute(
    SEQN, DAY,
    veg_cup_item        = dplyr::if_else(vegjuice_flag, 0,
                                         (dplyr::coalesce(V_TOTAL,0) - dplyr::coalesce(V_POTATO,0))) * w_mix,
    fruit_cup_item      = dplyr::if_else(juice100_flag, 0, dplyr::coalesce(F_TOTAL,0)) * w_mix,
    wholegr_g_item      = dplyr::coalesce(G_WHL,0) * 28.3495 * w_mix,
    nuts_serv_item      = dplyr::coalesce(M_NUTSD,0) * w_mix,
    legumes_serv_item   = dplyr::coalesce(LEGUMES,0) / 0.5 * w_mix,
    nutsleg_serv_item   = nuts_serv_item + legumes_serv_item,
    redproc_serv_item   = (dplyr::coalesce(M_MEAT,0) + dplyr::coalesce(M_FRANK,0)) / 3.527 * w_mix
  ) %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarise(
    veg_cup_eq_item        = sum(veg_cup_item,   na.rm=TRUE),
    fruit_cup_eq_item      = sum(fruit_cup_item, na.rm=TRUE),
    wholegr_g_item         = sum(wholegr_g_item, na.rm=TRUE),
    nuts_legumes_serv_item = sum(nutsleg_serv_item, na.rm=TRUE),
    redproc_serv_item      = sum(redproc_serv_item, na.rm=TRUE),
    .groups="drop"
  )
item_mixed_seqn <- mped_item

# 6) Beverages (SSB + 100% juice) from IFF ------------------------

ssb_person <- iff_joined %>%
  dplyr::mutate(
    dl    = stringr::str_to_lower(DESC %||% ""),
    fl_oz = dplyr::coalesce(GRAMS, 0)/29.5735,
    bev_core    = stringr::str_detect(dl, stringr::regex("soda|cola|soft\\s*drink|\\bpop\\b|lemonade|fruit\\s*(ade|drink|punch)|sports\\s*drink|energy\\s*drink|sweetened\\s*water|smoothie|frappuccino", TRUE)),
    bev_diet    = stringr::str_detect(dl, stringr::regex("\\bdiet\\b|sugar[- ]?free|unsweetened|zero\\b|low\\s*cal", TRUE)),
    bev_milkct  = stringr::str_detect(dl, stringr::regex("\\bmilk\\b|coffee|tea|cappuccino|latte|mocha|cocoa|hot\\s*chocolate", TRUE)),
    bev_reduced = stringr::str_detect(dl, stringr::regex("reduced\\s*sugar|less\\s*sugar|lower\\s*sugar|50%\\s*less\\s*sugar|\\blight\\b", TRUE)),
    ssb_full_flag = bev_core & !bev_diet & !bev_milkct & !bev_reduced,
    ssb_half_flag = bev_core & !bev_diet & !bev_milkct &  bev_reduced,
    is_juice    = stringr::str_detect(dl, stringr::regex("\\bjuice\\b", TRUE)),
    is_cocktail = stringr::str_detect(dl, stringr::regex("cocktail|\\bdrink\\b|beverage|ade\\b|punch\\b|nectar", TRUE)),
    juice_pack_context = stringr::str_detect(dl, stringr::regex("juice\\s*pack(ed)?|packed\\s*in\\s*(its\\s*|own\\s*)?juice|in\\s*(its\\s*|own\\s*)?juice", TRUE)),
    juice100_flag = is_juice & !juice_pack_context & !is_cocktail & !bev_diet & !bev_milkct
  ) %>%
  dplyr::summarise(
    ssb_serv_full     = sum((fl_oz/8) * as.numeric(ssb_full_flag), na.rm = TRUE),
    ssb_serv_half     = sum((fl_oz/8) * as.numeric(ssb_half_flag), na.rm = TRUE),
    ssb_serv          = ssb_serv_full + 0.5 * ssb_serv_half,
    fruit_juice_serv  = sum((fl_oz/4) * as.numeric(juice100_flag), na.rm = TRUE),
    ssb_juice_serv    = ssb_serv + fruit_juice_serv,
    .by = c(SEQN, DAY)
  ) %>%
  dplyr::group_by(SEQN) %>%
  dplyr::summarise(
    ssb_serv = sum(ssb_serv, na.rm=TRUE),
    fruit_juice_serv = sum(fruit_juice_serv, na.rm=TRUE),
    ssb_juice_serv = sum(ssb_juice_serv, na.rm=TRUE),
    .groups="drop"
  )
beverage_seqn <- ssb_person

# 7) Fallback from MPED day totals (pyr_tot*) ---------------------

infer_day_from_name  # ensure exists
tot_files <- list.files(dir$fped, pattern="pyr_tot.*\\.sas7bdat$", recursive=TRUE, full.names=TRUE)
if (!length(tot_files)) message("No pyr_tot*.sas7bdat files found under: ", dir$fped)

read_pyr_tot <- function(path){
  df <- safe_read_sas(path)
  if (is.null(df)) return(tibble())
  day_guess <- infer_day_from_name(path)
  tibble(
    SEQN     = as.numeric(safe_vec(df, c("SEQN"), NA_real_)),
    DAY      = as.numeric(safe_vec(df, c("DAY","DAYNO","DAYN"), day_guess)),
    V_TOTAL  = as.numeric(safe_vec(df, c("V_TOTAL","V_TOT"), NA_real_)),
    V_POTATO = as.numeric(safe_vec(df, c("V_POTATO","V_POTAT"), 0)),
    F_TOTAL  = as.numeric(safe_vec(df, c("F_TOTAL","F_TOT"), NA_real_)),
    G_WHL    = as.numeric(safe_vec(df, c("G_WHL","G_WHOLE","G_WG"), NA_real_)),   # oz-eq/day
    LEGUMES  = as.numeric(safe_vec(df, c("LEGUMES"), NA_real_)),                  # cup-eq/day
    M_NUTSD  = as.numeric(safe_vec(df, c("M_NUTSD","M_NUTS"), NA_real_)),         # oz-eq/day
    M_MEAT   = as.numeric(safe_vec(df, c("M_MEAT","MEAT_RED","M_RED"), NA_real_)),# oz-eq/day
    M_FRANK  = as.numeric(safe_vec(df, c("M_FRANK","MEAT_PROC","M_PROC"), NA_real_))
  ) %>% dplyr::filter(!is.na(SEQN))
}

pyr_tot_all <- if (length(tot_files)) purrr::map_dfr(tot_files, read_pyr_tot) else tibble()

pyr_tot_all_dedup <- pyr_tot_all %>% dplyr::filter(DAY %in% c(1, NA)) %>%
  dplyr::group_by(SEQN, DAY) %>%
  dplyr::summarise(dplyr::across(
    c(V_TOTAL, V_POTATO, F_TOTAL, G_WHL, LEGUMES, M_NUTSD, M_MEAT, M_FRANK),
    ~ if (all(is.na(.))) NA_real_ else mean(., na.rm = TRUE)
  ), .groups = "drop")

pyr_tot_d1 <- pyr_tot_all_dedup %>%
  dplyr::mutate(DAY = dplyr::coalesce(DAY, 1)) %>%
  dplyr::filter(DAY == 1) %>%
  dplyr::select(-DAY)

pyr_tot_fallback <- pyr_tot_d1 %>%
  dplyr::transmute(
    SEQN,
    veg_cup_eq_fb    = pmax(V_TOTAL - V_POTATO, 0),   # cups/day
    fruit_cup_eq_fb  = F_TOTAL,                       # cups/day
    wholegr_g_fb     = G_WHL * 28.3495,               # grams/day
    nuts_leg_serv_fb = M_NUTSD + (LEGUMES / 0.5),     # servings/day
    redproc_serv_fb  = (M_MEAT + M_FRANK) / 3.527     # servings/day
  )

# 8) Build AHEI input (prefer item-level; fallback to totals) -----
ahei_input <- nutrients_9904 %>%
  dplyr::filter(!is.na(WTDRD1) & WTDRD1 > 0, !is.na(energy_kcal) & energy_kcal > 0) %>%
  dplyr::left_join(item_mixed_seqn,  by = "SEQN") %>%
  dplyr::left_join(pyr_tot_fallback, by = "SEQN") %>%
  dplyr::left_join(beverage_seqn,    by = "SEQN") %>%
  dplyr::transmute(
    SEQN, WTDRD1, RIAGENDR,
    veg_cup_eq_final        = dplyr::coalesce(veg_cup_eq_item,        veg_cup_eq_fb),
    fruit_cup_eq_final      = dplyr::coalesce(fruit_cup_eq_item,      fruit_cup_eq_fb),
    wholegr_g_final         = dplyr::coalesce(wholegr_g_item,         wholegr_g_fb),
    nuts_legumes_serv_final = dplyr::coalesce(nuts_legumes_serv_item, nuts_leg_serv_fb),
    redproc_serv_final      = dplyr::coalesce(redproc_serv_item,      redproc_serv_fb),
    ssb_serv, fruit_juice_serv, ssb_juice_serv,
    energy_kcal, sodium_mg, alcohol_g, pufa_g, epa_g, dha_g
  )

summary(item_mixed_seqn$veg_cup_eq_item)
summary(ahei_input$SEQN)

# 9) AHEI scoring -------------------------------------------------

lin_pos  <- function(x, min0, max10) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (x - min0)/(max10 - min0))) * 10)
lin_rev  <- function(x, min10, max0) ifelse(is.na(x), NA_real_, pmax(0, pmin(1, (max0 - x)/(max0 - min10))) * 10)
per_1000 <- function(x, kcal) ifelse(is.na(x) | is.na(kcal) | kcal <= 0, NA_real_, x / (kcal / 1000))
E_REF <- 2000

# Vegetables
VEG_SERV_MAX   <- 5
SERVING_CUP    <- 0.5
VEG_CUPS_MAX   <- VEG_SERV_MAX * SERVING_CUP        # 2.5 cups/day
VEG_PER1000MAX <- VEG_CUPS_MAX / (E_REF/1000)       # 1.25 cups/1000

ahei_veg_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, energy_kcal,
                   veg_cups_per_1000 = per_1000(veg_cup_eq_final, energy_kcal),
                   ahei_veg = ifelse(is.na(veg_cups_per_1000), NA_real_,
                                     pmin(veg_cups_per_1000/VEG_PER1000MAX, 1) * 10))
report_component(ahei_veg_tbl, "ahei_veg", "Vegetables (per 1000 kcal)")

# Fruit
ahei_fruit_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, energy_kcal,
                   fruit_cups_per_1000 = per_1000(fruit_cup_eq_final, energy_kcal),
                   ahei_fruit = ifelse(is.na(fruit_cups_per_1000), NA_real_,
                                       pmin(fruit_cups_per_1000 / (2.0/(E_REF/1000)), 1) * 10))
report_component(ahei_fruit_tbl, "ahei_fruit", "Fruit (per 1000 kcal)")

# Whole grains (sex-specific max)
ahei_grain_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, RIAGENDR, energy_kcal,
                   wholegr_g_per_1000 = per_1000(wholegr_g_final, energy_kcal),
                   max_wholegr_per_1000 = dplyr::case_when(
                     RIAGENDR == 2 ~ 75/(E_REF/1000),  # 37.5 g/1000
                     RIAGENDR == 1 ~ 90/(E_REF/1000),  # 45.0 g/1000
                     TRUE ~ NA_real_
                   ),
                   ahei_wholegrains = ifelse(is.na(wholegr_g_per_1000) | is.na(max_wholegr_per_1000),
                                             NA_real_,
                                             pmin(wholegr_g_per_1000/max_wholegr_per_1000, 1) * 10))
report_component(ahei_grain_tbl, "ahei_wholegrains", "Whole grains (per 1000 kcal)")

# SSB + 100% juice (reverse)
ahei_ssb_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, energy_kcal,
                   ssb_juice_per_1000 = per_1000(ssb_juice_serv, energy_kcal),
                   ahei_ssb = lin_rev(ssb_juice_per_1000, 0, 1/(E_REF/1000)))
report_component(ahei_ssb_tbl, "ahei_ssb", "SSB + 100% juice (per 1000 kcal)")

# Nuts & legumes
ahei_nutsleg_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, energy_kcal,
                   nuts_leg_per_1000 = per_1000(nuts_legumes_serv_final, energy_kcal),
                   ahei_nutslegumes = lin_pos(nuts_leg_per_1000, 0, 1/(E_REF/1000)/2))
report_component(ahei_nutsleg_tbl, "ahei_nutslegumes", "Nuts & legumes (per 1000 kcal)")

# Red + processed meat (reverse)
MEAT_MAX_PER1000 <- 1.5/(E_REF/1000)  # = 0.75 servings/1000 at 2000 kcal
ahei_meat_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, energy_kcal,
                   redproc_per_1000 = per_1000(redproc_serv_final, energy_kcal),
                   ahei_redprocmeat = lin_rev(redproc_per_1000, 0, MEAT_MAX_PER1000))
report_component(ahei_meat_tbl, "ahei_redprocmeat", "Red + processed meat (per 1000 kcal)")

# Long-chain n-3 (EPA + DHA)
ahei_longn3_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, epa_g, dha_g,
                   lc_n3_mg = ifelse(is.na(epa_g) & is.na(dha_g), NA_real_,
                                     (dplyr::coalesce(epa_g,0) + dplyr::coalesce(dha_g,0)) * 1000),
                   ahei_longn3 = lin_pos(lc_n3_mg, 0, 250))
report_component(ahei_longn3_tbl, "ahei_longn3", "Long-chain n-3 (EPA+DHA)")

# PUFA % energy
ahei_pufa_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, pufa_g, energy_kcal,
                   pufa_energy_pct = ifelse(is.na(pufa_g) | is.na(energy_kcal), NA_real_, (pufa_g*9)/energy_kcal*100),
                   ahei_pufa = dplyr::case_when(
                     is.na(pufa_energy_pct) ~ NA_real_,
                     pufa_energy_pct <= 2   ~ 0,
                     pufa_energy_pct >= 10  ~ 10,
                     TRUE ~ (pufa_energy_pct - 2)/(10 - 2) * 10
                   ))
report_component(ahei_pufa_tbl, "ahei_pufa", "PUFA % energy")

# Alcohol (J-shape with 2.5 for non-drinkers)
ahei_alcohol_tbl <- ahei_input %>%
  dplyr::transmute(SEQN, WTDRD1, RIAGENDR, alcohol_g,
                   ahei_alcohol = dplyr::case_when(
                     is.na(alcohol_g) | is.na(RIAGENDR) ~ NA_real_,
                     RIAGENDR == 2 & alcohol_g <= 0 ~ 2.5,
                     RIAGENDR == 2 & alcohol_g > 0  & alcohol_g < 7  ~ 2.5 + (alcohol_g/7)*(10-2.5),
                     RIAGENDR == 2 & alcohol_g >= 7 & alcohol_g <= 21 ~ 10,
                     RIAGENDR == 2 & alcohol_g > 21 & alcohol_g < 35  ~ ((35 - alcohol_g)/(35 - 21))*10,
                     RIAGENDR == 2 & alcohol_g >= 35 ~ 0,
                     RIAGENDR == 1 & alcohol_g <= 0 ~ 2.5,
                     RIAGENDR == 1 & alcohol_g > 0  & alcohol_g < 7   ~ 2.5 + (alcohol_g/7)*(10-2.5),
                     RIAGENDR == 1 & alcohol_g >= 7 & alcohol_g <= 28  ~ 10,
                     RIAGENDR == 1 & alcohol_g > 28 & alcohol_g < 49   ~ ((49 - alcohol_g)/(49 - 28))*10,
                     RIAGENDR == 1 & alcohol_g >= 49 ~ 0
                   ))
report_component(ahei_alcohol_tbl, "ahei_alcohol", "Alcohol (AHEI)")

# Sodium — decile scoring on mg/1000 kcal (winsorized)
wtd_quantile <- function(x, w, probs = seq(0.1, 0.9, 0.1)) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]; w <- w[ok]; o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w); sapply(probs, function(p) x[which(cw >= p)[1]])
}
wt_vec <- if ("WTDRD1" %in% names(ahei_input)) ahei_input$WTDRD1 else rep(1, nrow(ahei_input))
sod_base <- ahei_input %>%
  dplyr::transmute(SEQN, sodium_mg, energy_kcal, WTDRD1 = wt_vec) %>%
  dplyr::filter(!is.na(sodium_mg), !is.na(energy_kcal), energy_kcal > 0) %>%
  dplyr::filter(energy_kcal >= 500, energy_kcal <= 6000) %>%
  dplyr::mutate(sod_den = sodium_mg / (energy_kcal/1000))

TRIM_P <- 0.005
trim_lo <- wtd_quantile(sod_base$sod_den, sod_base$WTDRD1, probs = TRIM_P)
trim_hi <- wtd_quantile(sod_base$sod_den, sod_base$WTDRD1, probs = 1 - TRIM_P)
sod_base <- sod_base %>% dplyr::mutate(sod_den_w = pmin(pmax(sod_den, trim_lo), trim_hi))
sod_dec <- wtd_quantile(sod_base$sod_den_w, sod_base$WTDRD1, probs = seq(0.1, 0.9, 0.1))
names(sod_dec) <- paste0("Q", seq(10, 90, 10))
cat("Sodium density winsorization bounds:", round(c(lo=trim_lo, hi=trim_hi),2), "\n")
cat("Sodium deciles:", round(sod_dec,2), "\n")

ahei_sodium_tbl <- sod_base %>%
  dplyr::mutate(
    ahei_sodium = dplyr::case_when(
      sod_den_w <= sod_dec[1] ~ 10,
      sod_den_w <= sod_dec[2] ~ 9,
      sod_den_w <= sod_dec[3] ~ 8,
      sod_den_w <= sod_dec[4] ~ 7,
      sod_den_w <= sod_dec[5] ~ 6,
      sod_den_w <= sod_dec[6] ~ 5,
      sod_den_w <= sod_dec[7] ~ 4,
      sod_den_w <= sod_dec[8] ~ 3,
      sod_den_w <= sod_dec[9] ~ 2,
      TRUE ~ 0
    )
  ) %>% dplyr::select(SEQN, ahei_sodium)

# 10) Combine components, require complete, summarise -------------

ahei_scores <- list(
  ahei_veg_tbl %>% dplyr::select(SEQN, ahei_veg),
  ahei_fruit_tbl %>% dplyr::select(SEQN, ahei_fruit),
  ahei_grain_tbl %>% dplyr::select(SEQN, ahei_wholegrains),
  ahei_ssb_tbl %>% dplyr::select(SEQN, ahei_ssb),
  ahei_nutsleg_tbl %>% dplyr::select(SEQN, ahei_nutslegumes),
  ahei_meat_tbl %>% dplyr::select(SEQN, ahei_redprocmeat),
  ahei_longn3_tbl %>% dplyr::select(SEQN, ahei_longn3),
  ahei_pufa_tbl %>% dplyr::select(SEQN, ahei_pufa),
  ahei_alcohol_tbl %>% dplyr::select(SEQN, ahei_alcohol),
  ahei_sodium_tbl %>% dplyr::select(SEQN, ahei_sodium)
) %>% purrr::reduce(dplyr::left_join, by="SEQN")

ahei_complete <- ahei_scores %>%
  dplyr::filter(if_all(starts_with("ahei_"), ~ !is.na(.))) %>%
  dplyr::mutate(ahei_total = rowSums(dplyr::across(starts_with("ahei_")), na.rm = FALSE))

cat("\nParticipants with complete AHEI components: ", nrow(ahei_complete), "\n")
print(summary(ahei_complete$ahei_total))

# Weighted means table
ahei_comp_w <- ahei_complete %>% dplyr::left_join(nutrients_9904 %>% dplyr::select(SEQN, WTDRD1), by="SEQN")
comp_cols <- names(ahei_complete)[grepl("^ahei_", names(ahei_complete))]
comp_cols <- setdiff(comp_cols, "ahei_total")

means_unw <- colMeans(ahei_comp_w[, comp_cols, drop=FALSE], na.rm=TRUE)
means_w1  <- sapply(comp_cols, function(cc) wmean(ahei_comp_w[[cc]], ahei_comp_w$WTDRD1))
means_w6  <- sapply(comp_cols, function(cc) wmean(ahei_comp_w[[cc]], ahei_comp_w$WTDRD1/3))

mean_table <- tibble::tibble(component = comp_cols) %>%
  dplyr::mutate(
    mean_unweighted = as.numeric(means_unw[component]),
    mean_wt_day1    = as.numeric(means_w1[component]),
    mean_wt_6yr     = as.numeric(means_w6[component])
  ) %>%
  dplyr::arrange(component)

mean_total <- tibble::tibble(
  component        = "ahei_total",
  mean_unweighted  = mean(ahei_comp_w$ahei_total, na.rm=TRUE),
  mean_wt_day1     = wmean(ahei_comp_w$ahei_total, ahei_comp_w$WTDRD1),
  mean_wt_6yr      = wmean(ahei_comp_w$ahei_total, ahei_comp_w$WTDRD1/3)
)
mean_table <- dplyr::bind_rows(mean_table, mean_total)
print(mean_table, n = nrow(mean_table))

# 11) Save + quick histogram --------------------------------------

readr::write_csv(ahei_complete, file.path(dir$output, "ahei_1999_2004_day1_item+fallback.csv"))

ggplot2::ggplot(ahei_complete, ggplot2::aes(x = ahei_total)) +
  ggplot2::geom_histogram(binwidth = 5) +
  ggplot2::labs(title = "AHEI (1999–2004, Day 1) — MPED IFF + mixture 0.5 (+ totals fallback)",
                x = "AHEI total", y = "Count") +
  ggplot2::theme_minimal()

# -------------------- diagnostics (optional) ---------------------
cat("\nOverlap checks:\n")
cat("  IFF SEQN  ∩ nutrients SEQN : ",
    length(intersect(iff_joined$SEQN, nutrients_9904$SEQN)), "\n")
cat("  IFF-eq SEQN ∩ nutrients SEQN: ",
    length(intersect(pyr_iff_all$SEQN, nutrients_9904$SEQN)), "\n")
cat("  Item SEQN  ∩ Fallback SEQN  : ",
    length(intersect(item_mixed_seqn$SEQN, pyr_tot_fallback$SEQN)), "\n")






