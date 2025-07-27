
# 1 Setup: Packages and Directories-------

# Set working directory to analysis folder
setwd("/Users/dengshuyue/Desktop/SDOH/analysis")

# Define directory structure
dir <- list()
dir$root    <- getwd()
dir$data    <- file.path(dir$root, "data")
dir$output  <- file.path(dir$root, "output")
dir$code    <- file.path(dir$root, "code")


# Load Required Packages------
# List of required packages
want <- c("dplyr", "survey", "foreign", "Hmisc", "data.table", "tidyr", 
          "tibble", "readr", "flextable", "officer", "usethis", "gert")

# Install any missing packages
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)

# Load all required packages
lapply(want, function(pkg) require(pkg, character.only = TRUE))
rm(want, need)


# 2 Generate Demographic Summary Table Using Survey Design-------

# 2.1. Load NHANES dataset -----
df <- fread(file.path(dir$data, "SODH_diet_mort_depr.csv"))
df <- df %>% filter(!is.na(wt10) & !is.na(sdmvstra) & !is.na(sdmvpsu))

# 2.2. Reassign BMI category -------
# revised the coding when merging/data cleaning so no need this step 

# df$bmic <- with(df, ifelse(bmi > 0 & bmi < 18.5, 0,
#                           ifelse(bmi >= 18.5 & bmi < 25, 1,
#                                  ifelse(bmi >= 25 & bmi < 30, 2,
#                                         ifelse(bmi >= 30, 3, 4)))))

table(df$pir)
table(df$SNAP)
table(df$bmic)

# 2.3. Define labeled variable names ------
variable_labels <- c(
  SEX = "Sex", RACE = "Race", EDU = "Education", pir = "Family income to poverty ratio", FS = "Food Insecurity",
  SNAP3 = "SNAP", SMK = "Smoking status", ALCG2 = "Drinking status", bmic = "BMI",
  DIABE = "Diabetes", CVD = "CVD", dm_rx = "DiabetesRx", chol_rx = "Cholestory",
  angina = "Angina", cancer = "Cancer", lung_disease = "Lung-disease", MORTSTAT = "Death",
  RIDAGEYR = "Age, years", met_hr = "Physical activity, median (SE)", hba1c = "HbA1c",
  sbp = "Systolic Blood Pressure", dbp = "Diastolic Blood Pressure",
  hdl = "High-Density Lipoprotein", ldl = "Low-Density Lipoprotein", tg = "Triglycerides",
  HEI2015_TOTAL_SCORE = "HEI2015", probable_depression = "Depression"
)

map_variable_labels_once <- function(var_vector) {
  sapply(var_vector, function(v) if (v %in% names(variable_labels)) variable_labels[[v]] else v, USE.NAMES = FALSE)
}

# 2.4. Rename df columns just once ------
df_labeled <- df
names(df_labeled) <- map_variable_labels_once(names(df_labeled))

###### Create formula-safe function
as_var_formula <- function(v) as.formula(paste0("~`", v, "`"))

# 2.5. Survey design -----
nhanes_design <- svydesign(
  id = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wt10,
  data = df_labeled,
  nest = TRUE
)

# 2.6. Categorical variables ------
cat_vars <- c("Sex", "Race", "Education", "Family income to poverty ratio", "SNAP", "Smoking status", "Drinking status", "BMI", "Food Insecurity")

cat_results <- lapply(cat_vars, function(v) {
  unweighted_tab <- table(df_labeled[[v]])
  tab <- svytable(as_var_formula(v), design = nhanes_design)
  pct <- prop.table(tab) * 100
  means <- svymean(as_var_formula(v), design = nhanes_design, na.rm = TRUE)
  se <- SE(means) * 100
  
  data.frame(
    Variable = v,
    Category = names(tab),
    Count = as.vector(unweighted_tab[names(tab)]),
    Mean_or_Percent = round(as.vector(pct), 1),
    SE = round(as.vector(se)[names(tab)], 2),
    Type = "Categorical"
  )
}) %>% bind_rows()

# 2.7. Binary variables ------
binary_vars <- c("Diabetes", "CVD", "DiabetesRx", "Cholestory", "Angina", "Cancer", "Lung-disease", "Death", "Depression")

binary_results <- lapply(binary_vars, function(v) {
  count <- sum(df_labeled[[v]] == 1, na.rm = TRUE)
  mean_obj <- svymean(as_var_formula(v), nhanes_design, na.rm = TRUE)
  
  data.frame(
    Variable = v,
    Category = "1",
    Count = count,
    Mean_or_Percent = round(100 * coef(mean_obj)[[1]], 1),
    SE = round(100 * SE(mean_obj)[[1]], 2),
    Type = "Binary"
  )
}) %>% bind_rows()

# 2.8. Continuous variables ------
cont_vars <- c("Age, years", "Physical activity, median (SE)", "HbA1c", "Systolic Blood Pressure", 
               "Diastolic Blood Pressure", "High-Density Lipoprotein", 
               "Low-Density Lipoprotein", "Triglycerides", "HEI2015")

#cont_results <- lapply(cont_vars, function(v) {
#  count <- sum(!is.na(df_labeled[[v]]))
#  stat <- svymean(as_var_formula(v), nhanes_design, na.rm = TRUE)
  
#  data.frame(
#    Variable = v,
#    Count = count,
#    Mean_or_Percent = round(coef(stat)[[1]], 2),
#    SE = round(SE(stat)[[1]], 2),
#    Type = "Continuous"
#  )
#}) %>% bind_rows()


cont_results <- lapply(cont_vars, function(v) {
  count <- sum(!is.na(df_labeled[[v]]))
  
  if (v == "Physical activity, median (SE)") {
    med <- svyquantile(as.formula(paste0("~`", v, "`")), nhanes_design, quantiles = 0.5, na.rm = TRUE)
    
    return(data.frame(
      Variable = v,
      Count = count,
      Mean_or_Percent = round(med[[1]], 2),
      SE = NA,
      Type = "Continuous (Median)"
    ))
    
  } else {
    stat <- svymean(as.formula(paste0("~`", v, "`")), nhanes_design, na.rm = TRUE)
    
    return(data.frame(
      Variable = v,
      Count = count,
      Mean_or_Percent = round(coef(stat)[[1]], 2),
      SE = round(SE(stat)[[1]], 2),
      Type = "Continuous"
    ))
  }
}) %>% bind_rows()


cont_results$Mean_or_Percent[cont_results$Variable == "Physical activity, median (SE)"] <-
  cont_results$Mean_or_Percent.quantile[cont_results$Variable == "Physical activity, median (SE)"]

cont_results$SE[cont_results$Variable == "Physical activity, median (SE)"] <-
  cont_results$Mean_or_Percent.se[cont_results$Variable == "Physical activity, median (SE)"]



# 2.9. Category labels ------
category_labels <- list(
  Sex = c("1" = "Male", "2" = "Female"),
  Race = c("1" = "Non-Hispanic White", "2" = "Non-Hispanic Black", "3" = "Hispanic", "4" = "Other"),
  Education = c("1" = "Less than high school", "2" = "High school or equivalent", "3" = "Some college", "4" = "College or above"),
  `Family income to poverty ratio` = c("1" = "< 1.3", "2" = "1.3~2.99", "3" = ">=3", "5" = "Missing"),
  `Food Insecurity` = c( "0" = "Insecure", "1" = "Secure"), 
  SNAP = c("0" = "Not participant", "1" = "Participant", "2" = "income eligible non-participant"),
  `Smoking status` = c("1" = "Nonsmokers", "2" = "Former smokers", "3" = "<15 cigarettes/day", "4" = "15-24.9 cigarettes/day", "5" = "≥ 25 cigarettes/day"),
  `Drinking status` = c("1" = "Nondrinkers", "2" = "Moderate drinker", "3" = "Heavy drinker", "4" = "Missing"),
  BMI = c("0" = "<18.5 kg/m2", "1" = "18-24.9 kg/m2", "2" = "25-29.9 kg/m2", "3" = "≥30 kg/m2")
)

cat_results <- cat_results %>%
  rowwise() %>%
  mutate(
    Category = if (Variable %in% names(category_labels)) {
      label_map <- category_labels[[Variable]]
      ifelse(Category %in% names(label_map), label_map[Category], Category)
    } else {
      Category
    }
  ) %>%
  ungroup()

# 2.10. Combine all ------
full_summary <- bind_rows(cat_results, binary_results, cont_results)

##### Formatting----------
bmi_levels_ordered <- c("<18.5 kg/m2", "18-24.9 kg/m2", "25-29.9 kg/m2", "≥30 kg/m2")
mean_sd_vars <- c("Age, years", "Physical activity, median (SE)", "HbA1c", 
                  "Systolic Blood Pressure", "Diastolic Blood Pressure",
                  "High-Density Lipoprotein", "Low-Density Lipoprotein",
                  "Triglycerides", "HEI2015")

full_summary <- full_summary %>%
  mutate(
    Category = if_else(Variable == "BMI", as.character(factor(Category, levels = bmi_levels_ordered)), Category),
    `Primary population` = ifelse(Variable %in% mean_sd_vars,
                                  paste0(Mean_or_Percent, " (", SE, ")"),
                                  paste0(Count, " (", Mean_or_Percent, "%)"))
  )

##### Add groupings ------
groupings <- list(
  "Sociodemographics" = c("Sex", "Race", "Education", "Family income to poverty ratio", "SNAP", "Food Insecurity"),
  "Health Behaviors" = c("Smoking status", "Drinking status", "Physical activity, median (SE)"),
  "Clinical Characteristics" = c("BMI", "HbA1c", "Diabetes", "DiabetesRx", "CVD", "Angina", "Cancer", "Cholestory", "Lung-disease", "Depression", "Death"),
  "Dietary & Physiologic Measures" = c("Age, years", "Systolic Blood Pressure", "Diastolic Blood Pressure", 
                                       "High-Density Lipoprotein", "Low-Density Lipoprotein", "Triglycerides", "HEI2015")
)

group_df <- enframe(groupings, name = "Group", value = "Variable") %>% unnest(cols = c(Variable))

full_summary <- full_summary %>%
  mutate(Variable = as.character(Variable)) %>%  # force as character, so those variable lose previous order
  left_join(group_df, by = "Variable") %>%
  mutate(
    Group = factor(Group, levels = names(groupings)),
    Variable = factor(Variable, levels = unlist(groupings))
  ) %>%
  arrange(Group, Variable) %>%  # force sorting
  mutate(
    Category = ifelse(as.character(Variable) %in% c(mean_sd_vars, binary_vars), "", Category),
    Variable = ifelse(duplicated(Variable), "", as.character(Variable)),
    Group = ifelse(duplicated(Group), "", as.character(Group))
  ) %>%
  select(Group, Variable, Category, `Primary population`)


# 2.11. Export and Display-------
write_csv(full_summary, "/Users/dengshuyue/Desktop/SDOH/analysis/output/demo_summary_table.csv")

demo_flex <- flextable(full_summary)

# Optional: Auto fit columns
# demo_flex <- autofit(demo_flex)

# Create a new Word document and add the flextable
doc <- read_docx() %>%
  body_add_par("Table 1. Descriptive characteristics of the primary population", style = "heading 1") %>%
  body_add_flextable(demo_flex)

# Save the document
print(doc, target = "/Users/dengshuyue/Desktop/SDOH/analysis/output/demo_summary_table.docx")


# end



