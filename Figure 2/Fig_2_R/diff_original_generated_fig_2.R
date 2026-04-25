r_gdp <- read.csv("est_1_R.csv")
r_pce <- read.csv("est_2_R.csv")
r_ffr <- read.csv("est_3_R.csv")

orig_gdp <- read.csv("../Fig_2/est_1.csv",
                     header = TRUE, skip = 1)
orig_pce <- read.csv("../Fig_2/est_2.csv",
                     header = TRUE, skip = 1)
orig_ffr <- read.csv("../Fig_2/est_3.csv",
                     header = TRUE, skip = 1)

cat("=== SUMMARY OF DIFFERENCES ===\n")
for (col in 1:3) {
  diff_gdp <- abs(r_gdp[1:20, col] - as.numeric(orig_gdp[1:20, col]))
  diff_pce <- abs(r_pce[1:20, col] - as.numeric(orig_pce[1:20, col]))
  diff_ffr <- abs(r_ffr[1:20, col] - as.numeric(orig_ffr[1:20, col]))
  cat(sprintf("Column %d - Max diff: GDP=%.4f, PCE=%.4f, FFR=%.4f\n",
              col,
              max(diff_gdp, na.rm=TRUE),
              max(diff_pce, na.rm=TRUE),
              max(diff_ffr, na.rm=TRUE)))
}

# Check column names of your results
cat("=== R Results ===\n")
print(colnames(r_gdp))
print(head(r_gdp, 3))

cat("\n=== Original EViews Results ===\n")
print(colnames(orig_gdp))
print(head(orig_gdp, 3))