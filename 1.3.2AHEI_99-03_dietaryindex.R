
# SECTION I ================================================================
# AHEI (1999–2004) using dietaryindex::AHEI_NHANES_MPED
# - Reads MPED per-100g equivalences from .sas7bdat or fixed-width .txt
# - Dummy WJFRT to satisfy internals
# - MPED-safe SSB handling (SSB_code = numeric(0))
# - Runs 99–00, 01–02 (Day 1 only), 03–04 (Day 1 + Day 2 if present)
# - Saves raw outputs + unweighted & survey-weighted summaries

# 0) Packages -----------------------------------------------------
if (!requireNamespace("dietaryindex", quietly = TRUE)) install.packages("dietaryindex")
suppressPackageStartupMessages({
  library(dietaryindex)
  library(dplyr)
  library(readr)
  library(haven)
  library(survey)
  library(tidyr)
  library(purrr)
})

# (Optional) If you hit reshape bugs, uncomment to install dev build:
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("jamesjiadazhan/dietaryindex")

# 1) Paths --------------------------------------------------------
root <- "/Users/dengshuyue/Desktop/SDOH/analysis"
paths <- list(
  root   = root,
  data   = file.path(root, "data"),
  output = file.path(root, "output"),
  nhanes = file.path(root, "data", "nhanes_deit"),
  fped   = file.path(root, "data", "fped")
)
invisible(lapply(paths, dir.create, showWarnings = FALSE, recursive = TRUE))

dir <- list(
  root   = getwd(),
  data   = file.path(getwd(), "data"),
  output = file.path(getwd(), "output"),
  code   = file.path(getwd(), "code"),
  fped   = file.path(getwd(), "data", "fped"),
  nhanes = file.path(getwd(), "data", "nhanes_deit")
)

# 2) Helpers ------------------------------------------------------

# Read MPED per-100g equivalences as df from either .sas7bdat or fixed-width .txt
read_equiv_any <- function(path_char_or_df) {
  if (is.data.frame(path_char_or_df)) {
    df <- path_char_or_df
  } else if (is.character(path_char_or_df) && file.exists(path_char_or_df)) {
    if (grepl("\\.sas7bdat$", path_char_or_df, ignore.case = TRUE)) {
      df <- haven::read_sas(path_char_or_df)
    } else if (grepl("\\.txt$", path_char_or_df, ignore.case = TRUE)) {
      widths <- c(8, 1, 6, rep(8, 32))
      col_names <- c(
        "DRDIFDCD","EQUIVFLAG","DRDIMC",
        "G_TOTAL","G_WHL","G_NWHL",
        "V_TOTAL","V_DRKGR","V_DPYEL","V_POTATO","V_STARCY","V_TOMATO","V_OTHER",
        "F_TOTAL","F_CITMLB","F_OTHER",
        "D_TOTAL","D_MILK","D_YOGURT","D_CHEESE",
        "M_MPF","M_MEAT","M_ORGAN","M_FRANK","M_POULT","M_FISH_HI","M_FISH_LO","M_EGG","M_SOY","M_NUTSD","LEGUMES",
        "DISCFAT_OIL","DISCFAT_SOL","ADD_SUG","A_BEV"
      )
      df <- utils::read.fwf(path_char_or_df, widths = widths, col.names = col_names)
    } else {
      stop("Unsupported MPED file type for: ", path_char_or_df)
    }
  } else {
    stop("MPED path not found or invalid.")
  }
  if (!("FOODCODE" %in% names(df))) df <- df %>% mutate(FOODCODE = as.numeric(DRDIFDCD))
  if (!("MODCODE"  %in% names(df)) && "DRDIMC" %in% names(df)) df <- df %>% mutate(MODCODE = as.numeric(DRDIMC))
  df
}

# Build minimal WJFRT expected by the function (FOODCODE, WHOLEFRT, FRTJUICE)
build_wjfrt_dummy_from_equiv <- function(equiv_df) {
  equiv_df %>%
    transmute(FOODCODE = as.numeric(FOODCODE)) %>%
    distinct() %>%
    mutate(WHOLEFRT = 0, FRTJUICE = 0)
}

# Unweighted + Survey-weighted summaries from AHEI output
summarize_ahei <- function(ahei_df, has_day2, f_dr1tot, f_dr2tot, f_demo, out_prefix, out_dir) {
  ahei_cols <- c(
    "AHEI_ALL","AHEI_NOETOH","AHEI_VEG","AHEI_FRT","AHEI_WGRAIN",
    "AHEI_NUTSLEG","AHEI_N3FAT","AHEI_PUFA","AHEI_SSB_FRTJ",
    "AHEI_REDPROC_MEAT","AHEI_SODIUM","AHEI_ALCOHOL"
  )
  ahei_cols <- intersect(ahei_cols, names(ahei_df))
  
  # Unweighted
  unw <- ahei_df %>%
    select(SEQN, all_of(ahei_cols)) %>%
    pivot_longer(-SEQN, names_to = "component", values_to = "value") %>%
    group_by(component) %>%
    summarise(
      n      = sum(!is.na(value)),
      mean   = mean(value, na.rm = TRUE),
      sd     = sd(value, na.rm = TRUE),
      p25    = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
      median = median(value, na.rm = TRUE),
      p75    = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
      min    = min(value, na.rm = TRUE),
      max    = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(match(component, ahei_cols))
  write_csv(unw, file.path(out_dir, paste0(out_prefix, "_unweighted_summary.csv")))
  
  # Weighted (WTDR2D if Day1+2; else WTDRD1), with design vars from DEMO
  if (has_day2) {
    nutr <- haven::read_xpt(f_dr2tot) %>% transmute(SEQN, WT = WTDR2D)
  } else {
    nutr <- haven::read_xpt(f_dr1tot) %>% transmute(SEQN, WT = WTDRD1)
  }
  demo <- haven::read_xpt(f_demo) %>% select(SEQN, SDMVPSU, SDMVSTRA)
  
  analytic <- ahei_df %>%
    select(SEQN, all_of(ahei_cols)) %>%
    inner_join(nutr, by = "SEQN") %>%
    inner_join(demo, by = "SEQN") %>%
    filter(!is.na(WT) & WT > 0)
  
  des <- survey::svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT, data = analytic, nest = TRUE)
  
  wmeans <- map_dfr(ahei_cols, function(v) {
    est <- try(svymean(as.formula(paste0("~", v)), design = des, na.rm = TRUE), silent = TRUE)
    if (inherits(est, "try-error")) {
      tibble(component = v, mean = NA_real_, SE = NA_real_, n = sum(!is.na(analytic[[v]])))
    } else {
      tibble(component = v, mean = as.numeric(coef(est)), SE = as.numeric(SE(est)), n = sum(!is.na(analytic[[v]])))
    }
  }) %>% arrange(match(component, ahei_cols))
  write_csv(wmeans, file.path(out_dir, paste0(out_prefix, "_weighted_means.csv")))
}

# Run one cycle end-to-end
run_cycle <- function(cycle, files, paths) {
  message("=== Cycle ", cycle, " ===")
  f_equiv_path <- file.path(paths$fped,   files$equiv)   # .sas7bdat or .txt
  f_demo       <- file.path(paths$nhanes, files$demo)
  f_dr1tot     <- file.path(paths$nhanes, files$dr1tot)
  f_dr1iff     <- file.path(paths$nhanes, files$dr1iff)
  f_dr2tot     <- if (!is.null(files$dr2tot)) file.path(paths$nhanes, files$dr2tot) else NULL
  f_dr2iff     <- if (!is.null(files$dr2iff)) file.path(paths$nhanes, files$dr2iff) else NULL
  
  exists_vec <- c(
    equiv = file.exists(f_equiv_path),
    demo  = file.exists(f_demo),
    dr1tot= file.exists(f_dr1tot),
    dr1iff= file.exists(f_dr1iff),
    dr2tot= ifelse(is.null(f_dr2tot), FALSE, file.exists(f_dr2tot)),
    dr2iff= ifelse(is.null(f_dr2iff), FALSE, file.exists(f_dr2iff))
  )
  print(exists_vec)
  
  if (!all(exists_vec[c("equiv","demo","dr1tot","dr1iff")])) {
    warning("Skipping ", cycle, " — missing required day 1 files.")
    return(invisible(NULL))
  }
  
  has_day2 <- (!is.null(f_dr2tot) && !is.null(f_dr2iff) &&
                 file.exists(f_dr2tot) && file.exists(f_dr2iff))
  
  # Read MPED equivalences (sas7bdat or txt) and build dummy WJFRT
  equiv_df   <- read_equiv_any(f_equiv_path)
  WJFRT_dummy <- build_wjfrt_dummy_from_equiv(equiv_df)
  mped_ssb_codes <- numeric(0)
  
  # Run dietaryindex
  if (has_day2) {
    message("Running AHEI_NHANES_MPED (Day 1 + Day 2)…")
    res <- AHEI_NHANES_MPED(
      MPED_PER_100_GRAM_PATH = equiv_df,
      WJFRT                  = WJFRT_dummy,
      NUTRIENT_PATH          = f_dr1tot,
      NUTRIENT_IND_PATH      = f_dr1iff,
      DEMO_PATH              = f_demo,
      NUTRIENT_PATH2         = f_dr2tot,
      NUTRIENT_IND_PATH2     = f_dr2iff,
      SSB_code               = mped_ssb_codes
    )
  } else {
    message("Running AHEI_NHANES_MPED (Day 1 only)…")
    res <- AHEI_NHANES_MPED(
      MPED_PER_100_GRAM_PATH = equiv_df,
      WJFRT                  = WJFRT_dummy,
      NUTRIENT_PATH          = f_dr1tot,
      NUTRIENT_IND_PATH      = f_dr1iff,
      DEMO_PATH              = f_demo,
      NUTRIENT_PATH2         = NULL,
      NUTRIENT_IND_PATH2     = NULL,
      SSB_code               = mped_ssb_codes
    )
  }
  
  # Save raw + summaries
  out_base <- paste0("ahei_mped_", gsub("[^0-9-]", "", cycle))
  if (!has_day2) out_base <- paste0(out_base, "_d1only")
  out_csv <- file.path(paths$output, paste0(out_base, ".csv"))
  write_csv(res, out_csv)
  message("Saved: ", out_csv)
  
  summarize_ahei(
    ahei_df   = res,
    has_day2  = has_day2,
    f_dr1tot  = f_dr1tot,
    f_dr2tot  = if (has_day2) f_dr2tot else f_dr1tot, # placeholder; not used when has_day2=FALSE
    f_demo    = f_demo,
    out_prefix= out_base,
    out_dir   = paths$output
  )
  
  invisible(res)
}

# 3) Cycle file maps ----------------------------------------------
# Adjust names if your local filenames differ (extensions are case-insensitive on most systems).
cycles <- list(
  "1999-2000" = list(
    equiv  = "equiv9400.txt",     # your 99–00 MPED per-100g file
    demo   = "DEMO.XPT",
    dr1tot = "DRXTOT.XPT",        # change to DR1TOT.XPT if that's what you have
    dr1iff = "DRXIFF.XPT",        # change to DR1IFF.XPT if that's what you have
    dr2tot = NULL,                # no Day 2 in public release
    dr2iff = NULL
  ),
  "2001-2002" = list(
    equiv  = "equiv0102.txt",     # txt is fine; read.fwf path supported
    demo   = "DEMO_B.XPT",
    dr1tot = "DRXTOT_B.XPT",      # <-- per your note
    dr1iff = "DRXIFF_B.XPT",      # <-- per your note
    dr2tot = NULL,                # no Day 2 in public release
    dr2iff = NULL
  ),
  "2003-2004" = list(
    equiv  = "equiv0304.sas7bdat",
    demo   = "DEMO_C.XPT",
    dr1tot = "DR1TOT_C.XPT",
    dr1iff = "DR1IFF_C.XPT",
    dr2tot = "DR2TOT_C.XPT",      # Day 2 starts here
    dr2iff = "DR2IFF_C.XPT"
  )
)

# 4) Run all cycles -----------------------------------------------
results <- lapply(names(cycles), function(cy) run_cycle(cy, cycles[[cy]], paths))
names(results) <- names(cycles)

# Done. Outputs written to: file.path(paths$output, *)

# merge to one long table with a cycle label
ahei_9904 <- bind_rows(results, .id = "cycle") %>%
  mutate(cycle = factor(cycle, levels = c("1999-2000","2001-2002","2003-2004")))

# quick checks
dim(ahei_9904)
table(ahei_9904$cycle)


# (optional) write it out
# write_csv(ahei_9904, "/Users/dengshuyue/Desktop/SDOH/analysis/output/ahei_1999_2004.csv")



# SECTION II ================================================================
# AHEI 1999–2004 (MPED-based) with:
#  - WJFRT built from joined_fcode.csv (FOODCODE → WHOLEFRT/FRTJUICE)
#  - SSB+100% juice recomputed from item-level fcode + kcal/1000
#  - Totals recomputed
# Produces: ahei_9904 and writes CSV

## 0) Packages ----------------------------------------------------
if (!requireNamespace("dietaryindex", quietly = TRUE)) install.packages("dietaryindex")
suppressPackageStartupMessages({
  library(dietaryindex)
  library(dplyr)
  library(readr)
  library(haven)
  library(survey)
  library(tidyr)
  library(purrr)
  library(stringr)
})

## 1) Paths -------------------------------------------------------
# Adjust root to your project if needed
root <- "/Users/dengshuyue/Desktop/SDOH/analysis"
paths <- list(
  root   = root,
  data   = file.path(root, "data"),
  output = file.path(root, "output"),
  nhanes = file.path(root, "data", "nhanes_deit"),
  fped   = file.path(root, "data", "fped")
)
invisible(lapply(paths, dir.create, showWarnings = FALSE, recursive = TRUE))

## 2) Cycle file maps ---------------------------------------------
cycles <- list(
  "1999-2000" = list(
    equiv  = "equiv9400.txt",
    demo   = "DEMO.XPT",
    dr1tot = "DRXTOT.XPT",     # or DR1TOT.XPT on your machine
    dr1iff = "DRXIFF.XPT",
    dr2tot = NULL,
    dr2iff = NULL
  ),
  "2001-2002" = list(
    equiv  = "equiv0102.txt",
    demo   = "DEMO_B.XPT",
    dr1tot = "DRXTOT_B.XPT",
    dr1iff = "DRXIFF_B.XPT",
    dr2tot = NULL,
    dr2iff = NULL
  ),
  "2003-2004" = list(
    equiv  = "equiv0304.sas7bdat",
    demo   = "DEMO_C.XPT",
    dr1tot = "DR1TOT_C.XPT",
    dr1iff = "DR1IFF_C.XPT",
    dr2tot = "DR2TOT_C.XPT",
    dr2iff = "DR2IFF_C.XPT"
  )
)

## 3) Load item-level joined_fcode and build WJFRT ----------------
# joined_fcode.csv must have: SEQN, DAY, ILINE, FOODCODE, GRAMS, cycle, dl, DESC
fcode_path <- file.path(paths$output, "joined_fcode.csv")
stopifnot(file.exists(fcode_path))
fcode <- readr::read_csv(fcode_path, show_col_types = FALSE)

# A canonical description per FOODCODE (most frequent non-missing), then flags
fcode_desc <- fcode %>%
  mutate(DESC = coalesce(DESC, "")) %>%
  filter(!is.na(FOODCODE)) %>%
  count(FOODCODE, DESC, sort = FALSE) %>%
  group_by(FOODCODE) %>%
  slice_max(n, with_ties = FALSE) %>%
  ungroup() %>%
  select(FOODCODE, DESC)

wjfrt <- fcode_desc %>%
  mutate(
    dl = str_to_lower(DESC),
    
    # fruit juice (100%): exclude cocktails/nectars/ade/“drink”, veg juices, diet/unsweetened, milk/coffee contexts
    is_juice        = str_detect(dl, "\\bjuice\\b"),
    is_cocktail     = str_detect(dl, "cocktail|\\bdrink\\b|beverage|\\bade\\b|punch\\b|nectar"),
    juice_ctx_pack  = str_detect(dl, "juice\\s*pack(ed)?|packed\\s*in\\s*(its\\s*|own\\s*)?juice|in\\s*(its\\s*|own\\s*)?juice"),
    veg_juice       = str_detect(dl, "(\\bv8\\b|vegetable\\s*juice|veg\\s*juice|tomato\\s*juice|carrot\\s*juice|beet\\s*juice|clamato|bloody\\s*mary)"),
    diet_unsweet    = str_detect(dl, "\\bdiet\\b|sugar[- ]?free|unsweetened|zero\\b|low\\s*cal"),
    milk_coffee     = str_detect(dl, "\\bmilk\\b|coffee|tea|cappuccino|latte|mocha|cocoa|hot\\s*chocolate"),
    
    FRTJUICE_flag   = is_juice & !is_cocktail & !juice_ctx_pack & !veg_juice & !diet_unsweet & !milk_coffee,
    
    # whole fruit (fresh/frozen/canned/dried); exclude spreads/sauces & juices
    fruit_word = str_detect(
      dl,
      paste0(
        "\\b(apple|banana|orange|grapefruit|tangerine|mandarin|",
        "grape|raisin|currant|prune|date|fig|pear|peach|nectarine|plum|apricot|",
        "cherry|pineapple|mango|papaya|kiwi|melon|cantaloupe|honeydew|watermelon|",
        "berry|berries|strawberry|blueberry|raspberry|blackberry|cranberry|",
        "pomegranate|guava|passion\\s*fruit|dragon\\s*fruit|lychee|persimmon)\\b"
      )
    ),
    not_spread_sauce = !str_detect(dl, "jam|jelly|preserve|fruit\\s*sauce|apple\\s*sauce|applesauce|pie\\s*fill|fruit\\s*butter\\b|sorbet|ice\\s*pop"),
    WHOLEFRT_flag    = fruit_word & not_spread_sauce & !FRTJUICE_flag
  ) %>%
  transmute(
    FOODCODE = as.numeric(FOODCODE),
    WHOLEFRT = as.integer(WHOLEFRT_flag),
    FRTJUICE = as.integer(FRTJUICE_flag)
  ) %>%
  distinct(FOODCODE, .keep_all = TRUE)

# Make sure types & exclusivity are clean
wjfrt_use <- wjfrt %>%
  transmute(
    FOODCODE = as.numeric(FOODCODE),
    WHOLEFRT = as.integer(WHOLEFRT > 0),
    FRTJUICE = as.integer(FRTJUICE > 0)
  ) %>%
  distinct(FOODCODE, .keep_all = TRUE)
stopifnot(all(c("FOODCODE","WHOLEFRT","FRTJUICE") %in% names(wjfrt_use)))
stopifnot(sum(wjfrt_use$WHOLEFRT == 1 & wjfrt_use$FRTJUICE == 1, na.rm = TRUE) == 0)

## 4) Helpers -----------------------------------------------------
# Read MPED equivalences from .sas7bdat or fixed-width txt (dietaryindex accepts a df)
read_equiv_any <- function(path_char_or_df) {
  if (is.data.frame(path_char_or_df)) return(path_char_or_df)
  stopifnot(is.character(path_char_or_df), file.exists(path_char_or_df))
  if (grepl("\\.sas7bdat$", path_char_or_df, ignore.case = TRUE)) {
    df <- haven::read_sas(path_char_or_df)
  } else if (grepl("\\.txt$", path_char_or_df, ignore.case = TRUE)) {
    widths <- c(8, 1, 6, rep(8, 32))
    col_names <- c(
      "DRDIFDCD","EQUIVFLAG","DRDIMC",
      "G_TOTAL","G_WHL","G_NWHL",
      "V_TOTAL","V_DRKGR","V_DPYEL","V_POTATO","V_STARCY","V_TOMATO","V_OTHER",
      "F_TOTAL","F_CITMLB","F_OTHER",
      "D_TOTAL","D_MILK","D_YOGURT","D_CHEESE",
      "M_MPF","M_MEAT","M_ORGAN","M_FRANK","M_POULT","M_FISH_HI","M_FISH_LO","M_EGG","M_SOY","M_NUTSD","LEGUMES",
      "DISCFAT_OIL","DISCFAT_SOL","ADD_SUG","A_BEV"
    )
    df <- utils::read.fwf(path_char_or_df, widths = widths, col.names = col_names)
  } else stop("Unsupported MPED file type: ", path_char_or_df)
  if (!("FOODCODE" %in% names(df))) df <- df %>% mutate(FOODCODE = as.numeric(DRDIFDCD))
  if (!("MODCODE"  %in% names(df)) && "DRDIMC" %in% names(df)) df <- df %>% mutate(MODCODE = as.numeric(DRDIMC))
  df
}

# energy helpers
per_1000 <- function(x, kcal) x / (kcal / 1000)
pick_first <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit)) hit[1] else NA_character_
}
get_energy <- function(files, paths) {
  f1 <- file.path(paths$nhanes, files$dr1tot)
  d1 <- haven::read_xpt(f1)
  e1 <- pick_first(names(d1), c("DR1TKCAL","DRXTKCAL","DR1IKCAL"))
  out <- d1 %>% transmute(SEQN, TKCAL1 = as.numeric(.data[[e1]]))
  if (!is.null(files$dr2tot)) {
    f2 <- file.path(paths$nhanes, files$dr2tot)
    d2 <- haven::read_xpt(f2)
    e2 <- pick_first(names(d2), c("DR2TKCAL","DRXTKCAL","DR2IKCAL"))
    out <- out %>%
      left_join(d2 %>% transmute(SEQN, TKCAL2 = as.numeric(.data[[e2]])), by = "SEQN") %>%
      mutate(energy_kcal = rowMeans(cbind(TKCAL1, TKCAL2), na.rm = TRUE)) %>%
      select(SEQN, energy_kcal)
  } else {
    out <- out %>% rename(energy_kcal = TKCAL1)
  }
  out
}

# beverages directly from fcode (handles missing dl by falling back to DESC)
bev_from_fcode <- function(fcode_df) {
  x <- fcode_df %>%
    mutate(
      dl    = ifelse(is.na(dl) | dl == "", str_to_lower(coalesce(DESC, "")), dl),
      fl_oz = coalesce(GRAMS, 0) / 29.5735
    )
  
  bev_core <- str_detect(
    x$dl,
    paste0(
      "(soft\\s*drink|carbonated\\s*(beverage|drink)|soda|cola|root\\s*beer|ginger\\s*ale|",
      "lemonade|orangeade|grapeade|\\b[a-z]+\\s+drink\\b|fruit\\s*drink|fruit\\s*punch|",
      "sports\\s*drink|energy\\s*drink|sweet(?:ened)?\\s*(tea|iced\\s*tea)|vitamin\\s*water|",
      "sweetened\\s*water)"
    )
  )
  bev_diet    <- str_detect(x$dl, "\\bdiet\\b|sugar[- ]?free|unsweetened|zero\\b|low\\s*cal")
  milk_coffee <- str_detect(x$dl, "\\bmilk\\b|coffee|cappuccino|latte|mocha|cocoa|hot\\s*chocolate")
  bev_reduced <- str_detect(x$dl, "reduced\\s*sugar|less\\s*sugar|lower\\s*sugar|50%\\s*less\\s*sugar|\\blight\\b")
  
  is_juice    <- str_detect(x$dl, "\\bjuice\\b")
  is_cocktail <- str_detect(x$dl, "cocktail|\\bdrink\\b|beverage|\\bade\\b|punch\\b|nectar")
  juice_ctx   <- str_detect(x$dl, "juice\\s*pack(ed)?|packed\\s*in\\s*(its\\s*|own\\s*)?juice|in\\s*(its\\s*|own\\s*)?juice")
  
  ssb_full <- bev_core & !bev_diet & !milk_coffee & !bev_reduced
  ssb_half <- bev_core & !bev_diet & !milk_coffee &  bev_reduced
  juice100 <- is_juice & !juice_ctx & !is_cocktail & !bev_diet & !milk_coffee
  
  per_day <- x %>%
    group_by(SEQN, cycle, DAY) %>%
    summarise(
      ssb_serv_full    = sum(fl_oz[ssb_full], na.rm = TRUE) / 8,
      ssb_serv_half    = sum(fl_oz[ssb_half], na.rm = TRUE) / 8,
      SSB_SERV         = ssb_serv_full + 0.5 * ssb_serv_half,
      FRUIT_JUICE_SERV = sum(fl_oz[juice100], na.rm = TRUE) / 4,
      SSB_JUICE_SERV   = SSB_SERV + FRUIT_JUICE_SERV,
      .groups = "drop_last"
    ) %>% ungroup()
  
  per_day %>%
    group_by(SEQN, cycle) %>%
    summarise(
      SSB_SERV_D1         = SSB_SERV[DAY == 1]              %>% { if (length(.)==0) NA_real_ else . },
      FRUIT_JUICE_SERV_D1 = FRUIT_JUICE_SERV[DAY == 1]      %>% { if (length(.)==0) NA_real_ else . },
      SSB_JUICE_SERV_D1   = SSB_JUICE_SERV[DAY == 1]        %>% { if (length(.)==0) NA_real_ else . },
      SSB_SERV_D2         = SSB_SERV[DAY == 2]              %>% { if (length(.)==0) NA_real_ else . },
      FRUIT_JUICE_SERV_D2 = FRUIT_JUICE_SERV[DAY == 2]      %>% { if (length(.)==0) NA_real_ else . },
      SSB_JUICE_SERV_D2   = SSB_JUICE_SERV[DAY == 2]        %>% { if (length(.)==0) NA_real_ else . },
      SSB_SERV_AVG         = mean(SSB_SERV, na.rm = TRUE),
      FRUIT_JUICE_SERV_AVG = mean(FRUIT_JUICE_SERV, na.rm = TRUE),
      SSB_JUICE_SERV_AVG   = mean(SSB_JUICE_SERV, na.rm = TRUE),
      .groups = "drop"
    )
}

score_ssb_ahei <- function(serv_per_1000, max_serv = 1.0/(2000/1000)) {
  pmax(0, pmin(1, (max_serv - serv_per_1000)/max_serv)) * 10
}

## 5) Run dietaryindex per cycle using your WJFRT ------------------
run_cycle <- function(cycle, files, paths, wjfrt_map) {
  message("=== Cycle ", cycle, " ===")
  f_equiv <- file.path(paths$fped,   files$equiv)
  f_demo  <- file.path(paths$nhanes, files$demo)
  f_dr1t  <- file.path(paths$nhanes, files$dr1tot)
  f_dr1i  <- file.path(paths$nhanes, files$dr1iff)
  f_dr2t  <- if (!is.null(files$dr2tot)) file.path(paths$nhanes, files$dr2tot) else NULL
  f_dr2i  <- if (!is.null(files$dr2iff)) file.path(paths$nhanes, files$dr2iff) else NULL
  
  stopifnot(file.exists(f_equiv), file.exists(f_demo), file.exists(f_dr1t), file.exists(f_dr1i))
  has_day2 <- (!is.null(f_dr2t) && !is.null(f_dr2i) && file.exists(f_dr2t) && file.exists(f_dr2i))
  
  equiv_df <- read_equiv_any(f_equiv)
  
  if (has_day2) {
    res <- AHEI_NHANES_MPED(
      MPED_PER_100_GRAM_PATH = equiv_df,
      WJFRT                  = wjfrt_map,
      NUTRIENT_PATH          = f_dr1t,
      NUTRIENT_IND_PATH      = f_dr1i,
      DEMO_PATH              = f_demo,
      NUTRIENT_PATH2         = f_dr2t,
      NUTRIENT_IND_PATH2     = f_dr2i,
      SSB_code               = numeric(0)
    )
  } else {
    res <- AHEI_NHANES_MPED(
      MPED_PER_100_GRAM_PATH = equiv_df,
      WJFRT                  = wjfrt_map,
      NUTRIENT_PATH          = f_dr1t,
      NUTRIENT_IND_PATH      = f_dr1i,
      DEMO_PATH              = f_demo,
      NUTRIENT_PATH2         = NULL,
      NUTRIENT_IND_PATH2     = NULL,
      SSB_code               = numeric(0)
    )
  }
  
  # save raw per-cycle (optional)
  out_base <- paste0("ahei_mped_", gsub("[^0-9-]", "", cycle), if (!has_day2) "_d1only" else "")
  readr::write_csv(res, file.path(paths$output, paste0(out_base, ".csv")))
  res
}

results <- lapply(names(cycles), function(cy) run_cycle(cy, cycles[[cy]], paths, wjfrt_use))
names(results) <- names(cycles)

# Merge to one table with cycle label
ahei_9904 <- bind_rows(results, .id = "cycle") %>%
  mutate(cycle = factor(cycle, levels = c("1999-2000","2001-2002","2003-2004")))

## 6) Fix SSB+100% juice using fcode + kcal/1000 -------------------
ssb_fix_all <- map_dfr(names(cycles), function(cy) {
  files  <- cycles[[cy]]
  energy <- get_energy(files, paths) %>% mutate(cycle = cy)
  bev    <- bev_from_fcode(fcode %>% filter(cycle == cy))
  
  bev %>%
    right_join(energy, by = c("SEQN","cycle")) %>%
    mutate(
      ssb_per1000        = per_1000(SSB_JUICE_SERV_AVG, energy_kcal),
      AHEI_SSB_FRTJ_FIX  = score_ssb_ahei(ssb_per1000)
    ) %>%
    select(cycle, SEQN, AHEI_SSB_FRTJ_FIX)
})

ahei_9904 <- ahei_9904 %>%
  left_join(ssb_fix_all, by = c("cycle","SEQN")) %>%
  mutate(AHEI_SSB_FRTJ = coalesce(AHEI_SSB_FRTJ_FIX, AHEI_SSB_FRTJ)) %>%
  select(-AHEI_SSB_FRTJ_FIX)

## 7) Recompute totals (with/without alcohol) ---------------------
comp_cols <- c(
  "AHEI_VEG","AHEI_FRT","AHEI_WGRAIN","AHEI_NUTSLEG",
  "AHEI_N3FAT","AHEI_PUFA","AHEI_SSB_FRTJ","AHEI_REDPROC_MEAT",
  "AHEI_SODIUM","AHEI_ALCOHOL"
)
comp_cols <- intersect(comp_cols, names(ahei_9904))

ahei_9904 <- ahei_9904 %>%
  rowwise() %>%
  mutate(
    AHEI_NOETOH_RECOMP = sum(c_across(all_of(setdiff(comp_cols, "AHEI_ALCOHOL"))), na.rm = TRUE),
    AHEI_ALL_RECOMP    = AHEI_NOETOH_RECOMP + coalesce(AHEI_ALCOHOL, 0)
  ) %>%
  ungroup()

## 8) Quick sanity & save -----------------------------------------
print(summary(ahei_9904$AHEI_FRT))
print(summary(ahei_9904$AHEI_SSB_FRTJ))
print(summary(ahei_9904$AHEI_ALL_RECOMP))

out_csv <- file.path(paths$output, "ahei_9904_wjfrt_ssbfix.csv")
readr::write_csv(ahei_9904, out_csv)
message("Saved: ", out_csv)

# ========= AHEI summary table (overall + by cycle) =========

# 1) pick all AHEI* columns that exist
ahei_cols <- names(ahei_9904)
ahei_cols <- ahei_cols[startsWith(ahei_cols, "AHEI_")]

# 2) helper to compute stats
summarise_stats <- function(df, id_cols = c("cycle")) {
  df %>%
    select(all_of(c(id_cols, ahei_cols))) %>%
    pivot_longer(cols = all_of(ahei_cols),
                 names_to = "component", values_to = "value") %>%
    group_by(across(all_of(id_cols)), component, .drop = FALSE) %>%
    summarise(
      n        = sum(!is.na(value)),
      missing  = sum(is.na(value)),
      mean     = mean(value, na.rm = TRUE),
      sd       = sd(value, na.rm = TRUE),
      p25      = as.numeric(quantile(value, 0.25, na.rm = TRUE)),
      median   = median(value, na.rm = TRUE),
      p75      = as.numeric(quantile(value, 0.75, na.rm = TRUE)),
      min      = suppressWarnings(min(value, na.rm = TRUE)),
      max      = suppressWarnings(max(value, na.rm = TRUE)),
      .groups  = "drop"
    ) %>%
    arrange(component, across(all_of(id_cols)))
}

# 3) build both views
ahei_summary_by_cycle <- summarise_stats(ahei_9904, id_cols = "cycle")
ahei_summary_overall  <- summarise_stats(ahei_9904, id_cols = character()) %>%
  mutate(cycle = "Overall") %>%
  relocate(cycle)

# 4) combine
ahei_summary <- bind_rows(ahei_summary_by_cycle, ahei_summary_overall) %>%
  mutate(component = factor(component, levels = sort(unique(component))))

# 5) preview & save
print(ahei_summary, n = 50)
# write_csv(ahei_summary, file.path(paths$output, "ahei_9904_summary_unweighted.csv"))













# SECTION III FPED era (2005–2018) ----------------------------------------------------
## AHEI_NHANES_FPED needs only the IND-level FPED + DR1IFF/DR2IFF paths.
## If you supply both Day 1 and Day 2, it returns the combined result automatically.
## See: https://jamesjiadazhan.github.io/dietaryindex_manual/reference/AHEI_NHANES_FPED.html


fped_cfg <- list(
  `2005-2006` = list(demo=file.path(dir$nhanes,"DEMO_D.XPT"),
                     fped1=file.path(dir$fped,"FPED_DR1IFF_D.SAS7BDAT"),
                     dr1i =file.path(dir$nhanes,"DR1IFF_D.XPT"),
                     fped2=file.path(dir$fped,"FPED_DR2IFF_D.SAS7BDAT"),
                     dr2i =file.path(dir$nhanes,"DR2IFF_D.XPT")),
  `2007-2008` = list(demo=file.path(dir$nhanes,"DEMO_E.XPT"),
                     fped1=file.path(dir$fped,"FPED_DR1IFF_E.SAS7BDAT"),
                     dr1i =file.path(dir$nhanes,"DR1IFF_E.XPT"),
                     fped2=file.path(dir$fped,"FPED_DR2IFF_E.SAS7BDAT"),
                     dr2i =file.path(dir$nhanes,"DR2IFF_E.XPT")),
  `2009-2010` = list(demo=file.path(dir$nhanes,"DEMO_F.XPT"),
                     fped1=file.path(dir$fped,"FPED_DR1IFF_F.SAS7BDAT"),
                     dr1i =file.path(dir$nhanes,"DR1IFF_F.XPT"),
                     fped2=file.path(dir$fped,"FPED_DR2IFF_F.SAS7BDAT"),
                     dr2i =file.path(dir$nhanes,"DR2IFF_F.XPT")),
  `2011-2012` = list(demo=file.path(dir$nhanes,"DEMO_G.XPT"),
                     fped1=file.path(dir$fped,"FPED_DR1IFF_G.SAS7BDAT"),
                     dr1i =file.path(dir$nhanes,"DR1IFF_G.XPT"),
                     fped2=file.path(dir$fped,"FPED_DR2IFF_G.SAS7BDAT"),
                     dr2i =file.path(dir$nhanes,"DR2IFF_G.XPT")),
  `2013-2014` = list(demo=file.path(dir$nhanes,"DEMO_H.XPT"),
                     fped1=file.path(dir$fped,"FPED_DR1IFF_H.SAS7BDAT"),
                     dr1i =file.path(dir$nhanes,"DR1IFF_H.XPT"),
                     fped2=file.path(dir$fped,"FPED_DR2IFF_H.SAS7BDAT"),
                     dr2i =file.path(dir$nhanes,"DR2IFF_H.XPT")),
  `2015-2016` = list(demo=file.path(dir$nhanes,"DEMO_I.XPT"),
                     fped1=file.path(dir$fped,"FPED_DR1IFF_I.SAS7BDAT"),
                     dr1i =file.path(dir$nhanes,"DR1IFF_I.XPT"),
                     fped2=file.path(dir$fped,"FPED_DR2IFF_I.SAS7BDAT"),
                     dr2i =file.path(dir$nhanes,"DR2IFF_I.XPT")),
  `2017-2018` = list(demo=file.path(dir$nhanes,"DEMO_J.XPT"),
                     fped1=file.path(dir$fped,"FPED_DR1IFF_J.SAS7BDAT"),
                     dr1i =file.path(dir$nhanes,"DR1IFF_J.XPT"),
                     fped2=file.path(dir$fped,"FPED_DR2IFF_J.SAS7BDAT"),
                     dr2i =file.path(dir$nhanes,"DR2IFF_J.XPT"))
)

run_fped_cycle <- function(cfg, label) {
  message("FPED ", label, " …")
  di <- dietaryindex::AHEI_NHANES_FPED(
    FPED_IND_PATH      = cfg$fped1,
    NUTRIENT_IND_PATH  = cfg$dr1i,
    FPED_IND_PATH2     = cfg$fped2,
    NUTRIENT_IND_PATH2 = cfg$dr2i
    # SSB_code = NULL  # (optional) pass your own FNDDS beverage codes
  ) %>% janitor::clean_names()
  
  di_adult <- drop_children_and_preg(di, cfg$demo)
  
  out <- file.path(dir$output, paste0("ahei_", gsub("-", "", label), "_fped.csv"))
  readr::write_csv(di_adult, out)
  invisibly(di_adult)
}

fped_list <- imap(fped_cfg, run_fped_cycle)
ahei_fped_0518 <- bind_rows(fped_list)

# Combine & quick sanity checks -------------------------------------------
ahei_all <- bind_rows(
  ahei_mped_9904 %>% mutate(cycle_group = "MPED_1999_2004"),
  ahei_fped_0518 %>% mutate(cycle_group = "FPED_2005_2018")
)

readr::write_csv(ahei_all, file.path(dir$output, "ahei_1999_2018_combined.csv"))

# Print quick distribution of total score (alcohol-included and excluded)
if (all(c("ahei_all","ahei_noetoh") %in% names(ahei_all))) {
  cat("\nAHEI_ALL summary:\n"); print(summary(ahei_all$ahei_all))
  cat("\nAHEI_NOETOH summary:\n"); print(summary(ahei_all$ahei_noetoh))
}

cat("\nWrote:\n  -", file.path(dir$output, "ahei_1999_2018_combined.csv"), "\n")




