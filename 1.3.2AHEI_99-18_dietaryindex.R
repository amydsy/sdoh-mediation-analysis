
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


# NEW SECTION III ================================================================
# AHEI (FPED era) 2005–2018 via dietaryindex::AHEI_NHANES_FPED
# One-shot script: installs packages (if needed), sets folders,
# runs all cycles, writes outputs, and summary tables.
# ------------------------------------------------
# Prereqs (local files expected):
#   NHANES IND-level IFF XPTs in  dir$nhanes  (e.g., DR1IFF_D.XPT, DR2IFF_D.XPT, …)
#   FPED IND-level SAS7BDAT in     dir$fped    (e.g., fped_dr1iff_0506.sas7bdat, …)
# Optional:
#   If you have 1999–2004 results already saved as
#     /output/ahei_1999_2004_mped.csv
#   they’ll be auto-loaded and combined; otherwise FPED-only outputs are produced.
# ================================================================

# 0) Setup --------------------------------------------------------
# (edit this if your root path is different)
root_path <- "/Users/dengshuyue/Desktop/SDOH/analysis"
setwd(root_path)

need <- c("dplyr","readr","haven","janitor","purrr","stringr","tidyr","rlang")
to_install <- need[!(need %in% rownames(installed.packages()))]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
lapply(need, require, character.only = TRUE)

if (!requireNamespace("dietaryindex", quietly = TRUE)) install.packages("dietaryindex", repos = "https://cloud.r-project.org")
library(dietaryindex)

dir <- list(
  root   = getwd(),
  data   = file.path(getwd(), "data"),
  nhanes = file.path(getwd(), "data", "nhanes_deit"),
  fped   = file.path(getwd(), "data", "fped"),
  output = file.path(getwd(), "output"),
  code   = file.path(getwd(), "code")
)
invisible(lapply(dir[c("data","nhanes","fped","output","code")],
                 dir.create, showWarnings = FALSE, recursive = TRUE))

# 1) Helpers ------------------------------------------------------
# Map 2-yr label -> NHANES suffix letter (for DR1TOT_*.XPT)
suffix_from_label <- function(label) {
  switch(label,
         "2005-2006" = "D", "2007-2008" = "E", "2009-2010" = "F",
         "2011-2012" = "G", "2013-2014" = "H", "2015-2016" = "I",
         "2017-2018" = "J", NA_character_
  )
}

# Adults (>=20y) & non-pregnant
drop_children_and_preg <- function(df, demo_path) {
  demo <- haven::read_xpt(demo_path) %>% janitor::clean_names()
  has_preg <- "ridexprg" %in% names(demo)
  demo2 <- demo %>% transmute(seqn, ridageyr, riagendr,
                              ridexprg = if (has_preg) ridexprg else NA_real_)
  df %>% janitor::clean_names() %>%
    dplyr::left_join(demo2, by = "seqn") %>%
    dplyr::filter(ridageyr >= 20, !(riagendr == 2 & ridexprg == 1)) %>%
    dplyr::select(-ridexprg)
}

# Attach WTDRD1 from DR1TOT_<suffix>.XPT (case-insensitive)
attach_wtdrd1_if_needed <- function(df, nhanes_dir, suffix_letter) {
  if ("wtdrd1" %in% names(df)) return(df)
  if (is.na(suffix_letter)) {
    warning("No suffix letter; cannot attach WTDRD1.")
    return(df)
  }
  patt <- paste0("(?i)^DR1TOT_", suffix_letter, "\\.XPT$")
  hit  <- list.files(nhanes_dir, pattern = patt, full.names = TRUE)
  if (length(hit) == 0) {
    warning("Missing DR1TOT_", suffix_letter, ".XPT — weighted summaries may be NA.")
    return(df)
  }
  wts <- haven::read_xpt(hit[1]) %>% janitor::clean_names() %>% dplyr::select(seqn, wtdrd1)
  out <- dplyr::left_join(df, wts, by = "seqn")
  if (!"wtdrd1" %in% names(out)) warning("WTDRD1 still missing after join.")
  out
}

# Weighted mean
w_mean <- function(x, w) {
  x <- as.numeric(x); w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# Filename tag: "2005-2006" -> "0506"
cycle_tag <- function(label) {
  m <- stringr::str_match(label, "^(\\d{4})-(\\d{4})$")
  if (is.na(m[1,1])) return(gsub("-", "", label))
  paste0(substr(m[1,2], 3, 4), substr(m[1,3], 3, 4))
}

# Paths & file guards
norm_if_exists <- function(p) if (file.exists(p)) normalizePath(p) else p
is_good_file   <- function(p, min_bytes = 20000L) {
  is.character(p) && length(p) == 1 && !is.na(p) && file.exists(p) && file.size(p) >= min_bytes
}

# 2) FPED config (numeric filenames) ------------------------------
fped_cfg <- list(
  `2005-2006` = list(
    demo = file.path(dir$nhanes, "DEMO_D.XPT"),
    fped1= file.path(dir$fped,  "fped_dr1iff_0506.sas7bdat"),
    dr1i = file.path(dir$nhanes, "DR1IFF_D.XPT"),
    fped2= file.path(dir$fped,  "fped_dr2iff_0506.sas7bdat"),
    dr2i = file.path(dir$nhanes, "DR2IFF_D.XPT")
  ),
  `2007-2008` = list(
    demo = file.path(dir$nhanes, "DEMO_E.XPT"),
    fped1= file.path(dir$fped,  "fped_dr1iff_0708.sas7bdat"),
    dr1i = file.path(dir$nhanes, "DR1IFF_E.XPT"),
    fped2= file.path(dir$fped,  "fped_dr2iff_0708.sas7bdat"),
    dr2i = file.path(dir$nhanes, "DR2IFF_E.XPT")
  ),
  `2009-2010` = list(
    demo = file.path(dir$nhanes, "DEMO_F.XPT"),
    fped1= file.path(dir$fped,  "fped_dr1iff_0910.sas7bdat"),
    dr1i = file.path(dir$nhanes, "DR1IFF_F.XPT"),
    fped2= file.path(dir$fped,  "fped_dr2iff_0910.sas7bdat"),
    dr2i = file.path(dir$nhanes, "DR2IFF_F.XPT")
  ),
  `2011-2012` = list(
    demo = file.path(dir$nhanes, "DEMO_G.XPT"),
    fped1= file.path(dir$fped,  "fped_dr1iff_1112.sas7bdat"),
    dr1i = file.path(dir$nhanes, "DR1IFF_G.XPT"),
    fped2= file.path(dir$fped,  "fped_dr2iff_1112.sas7bdat"),
    dr2i = file.path(dir$nhanes, "DR2IFF_G.XPT")
  ),
  `2013-2014` = list(
    demo = file.path(dir$nhanes, "DEMO_H.XPT"),
    fped1= file.path(dir$fped,  "fped_dr1iff_1314.sas7bdat"),
    dr1i = file.path(dir$nhanes, "DR1IFF_H.XPT"),
    fped2= file.path(dir$fped,  "fped_dr2iff_1314.sas7bdat"),
    dr2i = file.path(dir$nhanes, "DR2IFF_H.XPT")
  ),
  `2015-2016` = list(
    demo = file.path(dir$nhanes, "DEMO_I.XPT"),
    fped1= file.path(dir$fped,  "fped_dr1iff_1516.sas7bdat"),
    dr1i = file.path(dir$nhanes, "DR1IFF_I.XPT"),
    fped2= file.path(dir$fped,  "fped_dr2iff_1516.sas7bdat"),
    dr2i = file.path(dir$nhanes, "DR2IFF_I.XPT")
  ),
  `2017-2018` = list(
    demo = file.path(dir$nhanes, "DEMO_J.XPT"),
    fped1= file.path(dir$fped,  "fped_dr1iff_1718.sas7bdat"),
    dr1i = file.path(dir$nhanes, "DR1IFF_J.XPT"),
    fped2= file.path(dir$fped,  "fped_dr2iff_1718.sas7bdat"),
    dr2i = file.path(dir$nhanes, "DR2IFF_J.XPT")
  )
)

# 3) Per-cycle runner (robust; Day-2 fallback) -------------------
run_fped_cycle <- function(cfg, label) {
  message("FPED ", label, " …")
  # normalize
  demo  <- norm_if_exists(cfg$demo)
  dr1i  <- norm_if_exists(cfg$dr1i)
  dr2i  <- norm_if_exists(cfg$dr2i)
  fped1 <- norm_if_exists(cfg$fped1)
  fped2 <- norm_if_exists(cfg$fped2)
  
  # require Day-1; Day-2 optional
  stopifnot(is_good_file(fped1), is_good_file(dr1i))
  has_day2 <- is_good_file(fped2) && is_good_file(dr2i)
  
  # try Day-2; if it errors, fall back to Day-1 only
  di <- tryCatch({
    if (has_day2) {
      dietaryindex::AHEI_NHANES_FPED(
        FPED_IND_PATH      = fped1,
        NUTRIENT_IND_PATH  = dr1i,
        FPED_IND_PATH2     = fped2,
        NUTRIENT_IND_PATH2 = dr2i
      )
    } else {
      message("  • Day-2 not available/valid → using Day-1 only.")
      dietaryindex::AHEI_NHANES_FPED(
        FPED_IND_PATH      = fped1,
        NUTRIENT_IND_PATH  = dr1i
      )
    }
  }, error = function(e) {
    message("  • Day-2 read failed (", conditionMessage(e), ") → using Day-1 only.")
    dietaryindex::AHEI_NHANES_FPED(
      FPED_IND_PATH      = fped1,
      NUTRIENT_IND_PATH  = dr1i
    )
  }) %>% janitor::clean_names()
  
  # adult + non-pregnant
  di_adult <- drop_children_and_preg(di, demo)
  
  # ensure WTDRD1 using cycle label
  suf <- suffix_from_label(label)
  di_adult <- attach_wtdrd1_if_needed(di_adult, dir$nhanes, suf)
  
  # standardize totals + annotate Day-2
  di_adult <- di_adult %>%
    dplyr::mutate(cycle = label, used_day2 = has_day2) %>%
    dplyr::rename_with(~"AHEI_ALL",    dplyr::any_of(c("ahei_all","AHEI_ALL"))) %>%
    dplyr::rename_with(~"AHEI_NOETOH", dplyr::any_of(c("ahei_noetoh","AHEI_NOETOH")))
  
  # write per-cycle
  out <- file.path(dir$output, paste0("ahei_", cycle_tag(label), "_fped.csv"))
  readr::write_csv(di_adult, out)
  message("Saved: ", out)
  
  di_adult
}

# 4) Batch 2005–2018 ---------------------------------------------
ahei_fped_0518 <- purrr::imap(fped_cfg, run_fped_cycle) %>% dplyr::bind_rows()

# 5) Load 1999–2004 from saved CSV and combine -------------------
# Your saved file name:
#   ahei_9904_wjfrt_ssbfix.csv
# We'll look under `paths$output` if it exists, else `dir$output`.

output_dir <- if (exists("paths") && "output" %in% names(paths)) paths$output else dir$output
mped_csv   <- file.path(output_dir, "ahei_9904_wjfrt_ssbfix.csv")

if (!file.exists(mped_csv)) {
  stop("Could not find 1999–2004 file at: ", mped_csv,
       "\nUpdate the path/filename or place the file there and rerun.")
}

ahei_mped_9904 <- readr::read_csv(mped_csv, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  # standardize total columns just like FPED era
  dplyr::rename_with(~"AHEI_ALL",    dplyr::any_of(c("ahei_all","AHEI_ALL"))) %>%
  dplyr::rename_with(~"AHEI_NOETOH", dplyr::any_of(c("ahei_noetoh","AHEI_NOETOH"))) %>%
  # ensure we have a used_day2 flag (MPED era: not applicable)
  dplyr::mutate(used_day2 = NA)

# Bind MPED + FPED
ahei_all <- dplyr::bind_rows(
  ahei_mped_9904 %>% dplyr::mutate(cycle_group = "MPED_1999_2004"),
  ahei_fped_0518 %>% dplyr::mutate(cycle_group = "FPED_2005_2018")
)

# Write combined file
# readr::write_csv(ahei_all, file.path(output_dir, "ahei_1999_2018_combined.csv"))
# message("Wrote combined: ", file.path(output_dir, "ahei_1999_2018_combined.csv"))


summary(ahei_all$ahei_ssb_frtj)

summary(ahei_all)



# 6) Counts / weight sanity (combined) ---------------------------
counts_by_cycle_all <- ahei_all %>%
  dplyr::group_by(cycle, cycle_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    w_na = sum(is.na(wtdrd1)),
    w_na_pct = round(100 * mean(is.na(wtdrd1)), 2),
    used_day2 = dplyr::first(used_day2),
    .groups = "drop"
  ) %>% dplyr::arrange(cycle_group, cycle)

readr::write_csv(counts_by_cycle_all, file.path(output_dir, "ahei_counts_by_cycle_combined.csv"))
print(counts_by_cycle_all, n = Inf)






# ---- Attach WTDRD1 for MPED era (1999–2004) --------------------

# Map cycle -> possible NHANES suffix(es) used in older cycles
suffix_candidates_mped <- list(
  "1999-2000" = c("", "_A"),  # 99–00 often has no suffix; sometimes "_A"
  "2001-2002" = c("_B"),
  "2003-2004" = c("_C")
)

# Try to find a Day-1 totals file by cycle (DR1TOT_* or DRXTOT_*)
find_dr_day1_file <- function(nhanes_dir, cycle_label) {
  suffs <- suffix_candidates_mped[[cycle_label]]
  if (is.null(suffs)) suffs <- c("")  # fallback
  
  bases <- c("DR1TOT", "DRXTOT")
  for (base in bases) {
    for (suf in suffs) {
      patt <- paste0("(?i)^", base, suf, "\\.XPT$")
      hit <- list.files(nhanes_dir, pattern = patt, full.names = TRUE)
      if (length(hit) > 0) return(hit[1])
    }
  }
  
  # As a last resort: any DR1TOT/DRXTOT present at all
  hit <- list.files(nhanes_dir, pattern = "(?i)^(DR1TOT|DRXTOT).*\\.XPT$", full.names = TRUE)
  if (length(hit) > 0) return(hit[1])
  
  NA_character_
}

# Read weights from a DR1TOT/DRXTOT file (robust to column naming)
read_wtdrd1 <- function(path_xpt) {
  if (is.na(path_xpt) || !file.exists(path_xpt)) return(NULL)
  df <- haven::read_xpt(path_xpt) %>% janitor::clean_names()
  # prefer true day-1 dietary weight; fall back to interview 2-yr weight if truly absent (warn)
  weight_candidates <- c("wtdrd1", "wtdrd1d", "wtdrd1_", "wtint2yr")
  wcol <- intersect(weight_candidates, names(df))
  if (!length(wcol)) {
    warning("No WTDRD1-like column in: ", basename(path_xpt))
    return(NULL)
  }
  dplyr::select(df, seqn, wtdrd1 = !!rlang::sym(wcol[1]))
}

# Attach weights per cycle for MPED era
attach_mped_wts <- function(ahei_mped_9904, nhanes_dir) {
  cycs <- sort(unique(ahei_mped_9904$cycle))
  out_list <- vector("list", length(cycs))
  
  for (i in seq_along(cycs)) {
    cyc <- cycs[i]
    piece <- dplyr::filter(ahei_mped_9904, cycle == cyc)
    
    p <- find_dr_day1_file(nhanes_dir, cyc)
    if (is.na(p)) {
      warning("No DR1TOT/DRXTOT file found for ", cyc, "; leaving WTDRD1 as NA.")
      out_list[[i]] <- piece
      next
    }
    w <- read_wtdrd1(p)
    if (is.null(w)) {
      warning("Could not read WTDRD1 from ", basename(p), " for ", cyc, ".")
      out_list[[i]] <- piece
      next
    }
    piece2 <- dplyr::left_join(piece, w, by = "seqn")
    out_list[[i]] <- piece2
  }
  dplyr::bind_rows(out_list)
}

# ---- Run the MPED weight attach, then recompute combined -------

# a) Load your saved 99–04 CSV (already done in your flow)
#    ahei_mped_9904 <- readr::read_csv(file.path(output_dir, "ahei_9904_wjfrt_ssbfix.csv"), show_col_types = FALSE) %>%
#      janitor::clean_names() %>%
#      dplyr::rename_with(~"AHEI_ALL",    dplyr::any_of(c("ahei_all","AHEI_ALL"))) %>%
#      dplyr::rename_with(~"AHEI_NOETOH", dplyr::any_of(c("ahei_noetoh","AHEI_NOETOH"))) %>%
#      dplyr::mutate(used_day2 = NA)

# b) Attach weights for MPED era
ahei_mped_9904 <- attach_mped_wts(ahei_mped_9904, dir$nhanes)

# c) Combine MPED + FPED and write combined
ahei_all <- dplyr::bind_rows(
  ahei_mped_9904 %>% dplyr::mutate(cycle_group = "MPED_1999_2004"),
  ahei_fped_0518 %>% dplyr::mutate(cycle_group = "FPED_2005_2018")
)
readr::write_csv(ahei_all, file.path(output_dir, "ahei_1999_2018_combined.csv"))
message("Wrote combined: ", file.path(output_dir, "ahei_1999_2018_combined.csv"))

# d) Rebuild counts table (combined)
counts_by_cycle_all <- ahei_all %>%
  dplyr::group_by(cycle, cycle_group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    w_na = sum(is.na(wtdrd1)),
    w_na_pct = round(100 * mean(is.na(wtdrd1)), 2),
    used_day2 = dplyr::first(used_day2),
    .groups = "drop"
  ) %>% dplyr::arrange(cycle_group, cycle)
readr::write_csv(counts_by_cycle_all, file.path(output_dir, "ahei_counts_by_cycle_combined.csv"))
print(counts_by_cycle_all, n = Inf)

# e) Recompute combined summaries now that 99–04 has weights
#    (AHEI_ALL & AHEI_NOETOH by cycle + component means)
#    -- reuse your existing summary blocks below this point --



# 7) Summary tables (combined; AHEI_ALL as last row) -------------
# Weighted mean helper already defined earlier: w_mean()

# (a) AHEI_ALL by cycle (combined)
sum_all_comb <- ahei_all %>%
  dplyr::group_by(cycle) %>%
  dplyr::summarise(AHEI_ALL_mean_w = w_mean(AHEI_ALL, wtdrd1),
                   n = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(cycle)

sum_all_last_comb <- tibble::tibble(
  cycle = "AHEI_ALL",
  AHEI_ALL_mean_w = w_mean(ahei_all$AHEI_ALL, ahei_all$wtdrd1),
  n = nrow(ahei_all)
)

sum_all_out_comb <- dplyr::bind_rows(sum_all_comb, sum_all_last_comb)
# readr::write_csv(sum_all_out_comb, file.path(output_dir, "ahei_all_by_cycle_combined.csv"))

# (b) AHEI_NOETOH by cycle (combined)
sum_noetoh_comb <- ahei_all %>%
  dplyr::group_by(cycle) %>%
  dplyr::summarise(AHEI_NOETOH_mean_w = w_mean(AHEI_NOETOH, wtdrd1),
                   n = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(cycle)

sum_noetoh_last_comb <- tibble::tibble(
  cycle = "AHEI_ALL",
  AHEI_NOETOH_mean_w = w_mean(ahei_all$AHEI_NOETOH, ahei_all$wtdrd1),
  n = nrow(ahei_all)
)

sum_noetoh_out_comb <- dplyr::bind_rows(sum_noetoh_comb, sum_noetoh_last_comb)
# readr::write_csv(sum_noetoh_out_comb, file.path(output_dir, "ahei_noetoh_by_cycle_combined.csv"))


# (c) Component means by cycle (combined; weighted) --------------

# 1) Canonical -> explicit alias list
component_aliases <- list(
  VEG          = c("VEG","ahei_veg"),
  FRUIT        = c("FRUIT","ahei_frt"),
  WHOLEGRAIN   = c("WHOLEGRAIN","ahei_wgrain"),
  NUTSLEG      = c("NUTSLEG","ahei_nutsleg"),
  N3           = c("N3","ahei_n3fat"),
  PUFA         = c("PUFA","ahei_pufa"),
  TFA          = c("TFA","ahei_tfa"),           # may exist only in MPED era
  SSB          = c("SSB","ahei_ssb","ahei_ssb_frtj"),
  REDPROCMEAT  = c("REDPROCMEAT","ahei_redproc"),
  ALCOHOL      = c("ALCOHOL","ahei_alcohol")
)

# 2) Fallback regex patterns (if explicit aliases aren’t present)
pattern_map <- list(
  VEG          = "(^|_)(veg|vegetable|vegetables)$",
  FRUIT        = "(^|_)(frt|fruit|fruits)$",
  WHOLEGRAIN   = "(^|_)(w(hole)?_?grain|wgrain|wholegrains?)$",
  NUTSLEG      = "(^|_)(nuts?(_?leg(umes)?)?|nutsleg)$",
  N3           = "(^|_)(n3(_?fat)?|epa_dha|omega3|om?ega?3)$",
  PUFA         = "(^|_)pufa$",
  TFA          = "(^|_)tfa$",
  SSB          = "(^|_)(ssb(_?frtj)?|sugar(_|)?sweet(en(ed)?)?_?bev(erage)?s?)$",
  REDPROCMEAT  = "(^|_)(red(_?proc)?_?meat|redproc)$",
  ALCOHOL      = "(^|_)alcohol$"
)

exclude_cols <- c("AHEI_ALL","AHEI_NOETOH","cycle","wtdrd1","cycle_group","seqn","used_day2")
all_names <- names(ahei_all)
candidates <- setdiff(all_names, exclude_cols)

# 3) Build alias->canonical map from explicit aliases that actually exist
explicit_list <- purrr::imap(component_aliases, function(aliases, canon) {
  hits <- intersect(aliases, candidates)
  if (length(hits)) setNames(rep(canon, length(hits)), hits) else NULL
}) |> purrr::compact()

alias_to_canon_explicit <- if (length(explicit_list)) do.call(c, explicit_list) else c()

# 4) For components still missing, find one column via regex
missing_components <- setdiff(names(component_aliases), unname(unique(alias_to_canon_explicit)))
regex_list <- lapply(missing_components, function(canon) {
  pat <- pattern_map[[canon]]
  if (is.null(pat)) return(NULL)
  pool <- setdiff(candidates, names(alias_to_canon_explicit))
  idx <- which(grepl(pat, pool, ignore.case = TRUE, perl = TRUE))
  if (length(idx)) setNames(c(canon), pool[idx[1]]) else NULL   # named vector: alias -> canon
}) |> purrr::compact()

alias_to_canon_regex <- if (length(regex_list)) do.call(c, regex_list) else c()

# 5) Merge maps (names = aliases, values = canonical), no outer name pollution
alias_to_canon <- c(alias_to_canon_explicit, alias_to_canon_regex)

# Safety/diagnostics
if (!length(alias_to_canon)) {
  cat("\nNo component columns detected. Available columns:\n",
      paste0(" - ", sort(all_names)), sep = "\n")
  stop("Could not detect component columns in `ahei_all`. Adjust aliases/patterns.")
}

# Optional: inspect what we mapped
# print(alias_to_canon)

# 6) Save mapping for transparency
alias_map_df <- tibble::tibble(
  alias = names(alias_to_canon),
  component = unname(alias_to_canon)
)
# readr::write_csv(alias_map_df, file.path(output_dir, "ahei_component_alias_mapping_combined.csv"))

# 7) Pivot long and compute weighted means (robust to prefixed keys) ----

# If alias_to_canon exists from a prior run, clean its names (strip "CANON." prefix)
if (exists("alias_to_canon")) {
  alias_keys <- names(alias_to_canon)
  if (!is.null(alias_keys)) {
    # remove any "<prefix>." before the real column name
    alias_keys_clean <- sub("^[^.]+\\.", "", alias_keys)
    names(alias_to_canon) <- alias_keys_clean
    
    # de-duplicate keys after cleaning (keep first)
    alias_to_canon <- alias_to_canon[!duplicated(names(alias_to_canon))]
  }
} else {
  stop("alias_to_canon not found. Re-run the alias/pattern mapping block first.")
}

# Keep only aliases that actually exist as columns in ahei_all
cols <- intersect(names(ahei_all), names(alias_to_canon))
if (length(cols) == 0L) {
  cat("\nNo usable component columns after cleaning.\nColumns in ahei_all include:\n",
      paste0(" - ", sort(names(ahei_all))), sep = "\n")
  stop("No matching component columns found in `ahei_all`.")
}

# Reduce the map to those columns
alias_to_canon <- alias_to_canon[cols]

# Do the pivot + weighted summary
comp_df_combined <- ahei_all %>%
  dplyr::select(cycle, wtdrd1, dplyr::all_of(cols)) %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(cols),
    names_to  = "alias",
    values_to = "score"
  ) %>%
  dplyr::mutate(
    component = unname(alias_to_canon[alias]),
    score = as.numeric(score)
  ) %>%
  dplyr::group_by(cycle, component) %>%
  dplyr::summarise(
    mean_w = w_mean(score, wtdrd1),
    n      = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(cycle, component)

readr::write_csv(
  comp_df_combined,
  file.path(output_dir, "ahei_component_means_by_cycle_combined.csv")
)

# (optional) also write the cleaned alias -> canonical mapping actually used
alias_map_df <- tibble::tibble(
  alias = names(alias_to_canon),
  component = unname(alias_to_canon)
)
# readr::write_csv(alias_map_df, file.path(output_dir, "ahei_component_alias_mapping_combined.csv"))

message("Wrote: ",
        file.path(output_dir, "ahei_component_means_by_cycle_combined.csv"),
        " and alias map: ",
        file.path(output_dir, "ahei_component_alias_mapping_combined.csv"))




# push to git
# setwd("/Users/dengshuyue/Desktop/SDOH/analysis/code")




