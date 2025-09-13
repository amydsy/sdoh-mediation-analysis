
# 
install.packages("duckdb")   # fast, no compilers needed
library(DBI)
library(duckdb)

p <- "/Users/dengshuyue/Desktop/SDOH/analysis/output/cov_core_1999_2023.parquet"

con <- dbConnect(duckdb(), dbdir=":memory:")



# Create a convenient view so you don't repeat read_parquet(...)
dbExecute(con, sprintf("CREATE VIEW core AS SELECT * FROM read_parquet('%s')", p))

# 1) Quick schema (column names + types)-----
dbGetQuery(con, "
  SELECT column_name, data_type
  FROM information_schema.columns
  WHERE table_name = 'core'
  ORDER BY ordinal_position
")

# 1.1) Peek a few rows-------
head_core <- dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') LIMIT 10", p))
print(head_core)


# 2) Row count and distinct cycles
dbGetQuery(con, "SELECT COUNT(*) AS n FROM core")
dbGetQuery(con, "SELECT SDDSRVYR, COUNT(*) AS rows FROM core GROUP BY SDDSRVYR ORDER BY SDDSRVYR")

# 3) PA coverage by cycle (non-missing LTPA)
dbGetQuery(con, "
  SELECT
    SDDSRVYR,
    COUNT(*)                      AS rows,
    COUNT(LTPA)                   AS rows_with_pa,
    ROUND(100.0 * COUNT(LTPA)/COUNT(*), 1) AS pa_cov_pct
  FROM core
  GROUP BY SDDSRVYR
  ORDER BY SDDSRVYR
")

# 4) Peek one cycle only (e.g., cycle 66)
dbGetQuery(con, "
  SELECT SEQN, SDDSRVYR, WTMEC2YR, LTPA, METSCORE, BMI, BMI_CLAS
  FROM core
  WHERE SDDSRVYR = 66
  LIMIT 10
")

# 5) Pull a small subset into an R data.frame
core_small <- dbGetQuery(con, "
  SELECT SEQN, SDDSRVYR, SMK_STATUS, ALCOHOL_CAT, LTPA, METSCORE, BMI
  FROM core
  WHERE SDDSRVYR IN (12,66)
  LIMIT 1000
")


# 6) (Optional) Export to CSV if you want a lightweight copy
dbExecute(con, sprintf("
  COPY (SELECT * FROM core WHERE SDDSRVYR IN (12,66))
  TO '%s' (FORMAT CSV, HEADER, DELIMITER ',')
", "~/localdata/core_2011_2023_subset.csv"))

# Done
dbDisconnect(con, shutdown = TRUE)



