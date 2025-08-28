
# ================================================================
# AHEI (1999–2004) using dietaryindex::AHEI_NHANES_MPED
# - Reads MPED per-100g equivalences from .sas7bdat or fixed-width .txt
# - Dummy WJFRT to satisfy internals
# - MPED-safe SSB handling (SSB_code = numeric(0))
# - Runs 99–00, 01–02 (Day 1 only), 03–04 (Day 1 + Day 2 if present)
# - Saves raw outputs + unweighted & survey-weighted summaries
# ================================================================

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


# 5) Survey-weighted means (and SE) by cycle---------

suppressPackageStartupMessages({ library(survey) })

# helper: get weight file & variable for a cycle using your 'cycles' map
get_weight_info <- function(cycle_name) {
  cf <- cycles[[cycle_name]]
  stopifnot(!is.null(cf))
  if (is.null(cf$dr2tot)) {
    list(weight_file = file.path(paths$nhanes, cf$dr1tot), weight_var = "WTDRD1")
  } else {
    list(weight_file = file.path(paths$nhanes, cf$dr2tot), weight_var = "WTDR2D")
  }
}

weighted_by_cycle <- imap_dfr(results, ~{
  cyc <- .y; df <- .x
  cf <- cycles[[cyc]]
  info <- get_weight_info(cyc)
  
  nutr <- haven::read_xpt(info$weight_file) %>%
    transmute(SEQN, WT = .data[[info$weight_var]])
  demo <- haven::read_xpt(file.path(paths$nhanes, cf$demo)) %>%
    select(SEQN, SDMVPSU, SDMVSTRA)
  
  analytic <- df %>%
    select(SEQN, any_of(ahei_cols)) %>%
    inner_join(nutr, by = "SEQN") %>%
    inner_join(demo, by = "SEQN") %>%
    filter(!is.na(WT) & WT > 0)
  
  des <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT,
    data = analytic, nest = TRUE
  )
  
  purrr::map_dfr(intersect(ahei_cols, names(df)), function(v) {
    est <- try(svymean(as.formula(paste0("~", v)), design = des, na.rm = TRUE), silent = TRUE)
    if (inherits(est, "try-error")) {
      tibble(cycle = cyc, component = v, mean = NA_real_, SE = NA_real_, n = sum(!is.na(analytic[[v]])))
    } else {
      tibble(cycle = cyc, component = v,
             mean = as.numeric(coef(est)), SE = as.numeric(SE(est)),
             n = sum(!is.na(analytic[[v]])))
    }
  })
}) %>%
  arrange(cycle, match(component, ahei_cols))

weighted_by_cycle %>% print(n = 50)
# Save if you want:
# readr::write_csv(weighted_by_cycle, file.path(paths$output, "ahei_weighted_by_cycle.csv"))





















data("NHANES_20032004")

AHEI_NHANES_MPED(MPED_PER_100_GRAM_PATH = NHANES_20032004$MPED_PER_100_GRAM, WJFRT = NHANES_20032004$WJFRT, NUTRIENT_PATH = NHANES_20032004$NUTRIENT, NUTRIENT_IND_PATH = NHANES_20032004$NUTRIENT_IND, DEMO_PATH = NHANES_20032004$DEMO, NUTRIENT_PATH2 = NHANES_20032004$NUTRIENT2, NUTRIENT_IND_PATH2 = NHANES_20032004$NUTRIENT_IND2)




# 1) Load pkg + example bundle
suppressPackageStartupMessages({
  library(dietaryindex); library(dplyr); library(readr); library(survey)
})
data("NHANES_20032004")  # gives: MPED_PER_100_GRAM, WJFRT, NUTRIENT, NUTRIENT_IND, DEMO, NUTRIENT2, NUTRIENT_IND2

# 2) Run AHEI for 2003–2004 (Day1+Day2 from the bundle)
#    IMPORTANT: pass empty SSB_code to avoid the package’s 2017–18 default SSB list
ahei_0304 <- AHEI_NHANES_MPED(
  MPED_PER_100_GRAM_PATH = NHANES_20032004$MPED_PER_100_GRAM,
  WJFRT                  = NHANES_20032004$WJFRT,
  NUTRIENT_PATH          = NHANES_20032004$NUTRIENT,
  NUTRIENT_IND_PATH      = NHANES_20032004$NUTRIENT_IND,
  DEMO_PATH              = NHANES_20032004$DEMO,
  NUTRIENT_PATH2         = NHANES_20032004$NUTRIENT2,
  NUTRIENT_IND_PATH2     = NHANES_20032004$NUTRIENT_IND2,
  SSB_code               = numeric(0)
)

# 3) Save raw results
out_dir <- "output"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write_csv(ahei_0304, file.path(out_dir, "ahei_mped_2003_2004_from_bundle.csv"))

# 4) Quick survey-weighted means (uses WTDR2D for 2003–2004)
#    The bundle’s NUTRIENT2 has WTDR2D; DEMO has SDMVPSU/SDMVSTRA.
nutr_wt <- NHANES_20032004$NUTRIENT2 %>% dplyr::transmute(SEQN, WT = WTDR2D)
demo    <- NHANES_20032004$DEMO      %>% dplyr::select(SEQN, SDMVPSU, SDMVSTRA)

analytic <- ahei_0304 %>%
  inner_join(nutr_wt, by = "SEQN") %>%
  inner_join(demo,    by = "SEQN") %>%
  filter(!is.na(WT) & WT > 0)

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT, data = analytic, nest = TRUE)

ahei_cols <- c("AHEI_ALL","AHEI_NOETOH","AHEI_VEG","AHEI_FRT","AHEI_WGRAIN",
               "AHEI_NUTSLEG","AHEI_N3FAT","AHEI_PUFA","AHEI_SSB_FRTJ",
               "AHEI_REDPROC_MEAT","AHEI_SODIUM","AHEI_ALCOHOL")
ahei_cols <- intersect(ahei_cols, names(analytic))

weighted_means <- do.call(rbind, lapply(ahei_cols, function(v){
  est <- svymean(as.formula(paste0("~", v)), des, na.rm = TRUE)
  data.frame(component = v, mean = as.numeric(coef(est)), SE = as.numeric(SE(est)))
}))
readr::write_csv(weighted_means, file.path(out_dir, "ahei_mped_2003_2004_weighted_means_from_bundle.csv"))
print(weighted_means)



# FPED era (2005–2018) ----------------------------------------------------
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




