# =============================================================================
# Compute_Ct.R
# Diebold-Yilmaz total connectedness index Ct (housing-only, 51-state VAR).
#
# Method:
#   - Rolling window (10y = 120 months), step 1 month.
#   - Within each window: equation-by-equation Elastic-Net VAR(1) on the 51
#     state-level housing returns; generalised FEVD with 10-step horizon;
#     Diebold-Yilmaz (2012) spillover table; scalar total connectedness.
#
# All numerical helpers live in connectedness_helpers.R (sourced below).
#
# INPUT:  EI_DATA.csv  (51 state housing returns in columns 4:54)
# OUTPUT: Ct_monthly.csv  (one Ct value per rolling window, dated at window
#         centre -> series spans 1981Q1 to 2016Q4)
# =============================================================================

rm(list = ls())

# Load all helper functions (LVAR, Bcoef, Acoef, Phi, res_fn,
# fevd_generalised, sptable2012). library(glmnet) is loaded inside.
source("connectedness_helpers.R")


# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
raw <- read.csv("EI_DATA.csv", header = TRUE)

# Columns 4-54: 51 state housing return series (AK..WY)
sym        <- raw[, c(4:60)]
data_House <- as.matrix(sym[, c(1:51)])

# 'data' is referenced globally by sptable2012() to label rows/columns.
# Required for the helper to attach state names; do not rename.
dat  <- data_House
data <- sym[, c(1:51)]

cat("Data loaded:", nrow(dat), "monthly observations,", ncol(dat), "states\n")
cat("Sample: Year", raw$Year[1], "Month", raw$Month[1],
    "to Year", raw$Year[nrow(raw)], "Month", raw$Month[nrow(raw)], "\n\n")


# ---------------------------------------------------------------------------
# 2. Rolling window parameters
# ---------------------------------------------------------------------------
window  <- 120   # 10-year rolling window (months)
p       <- 1     # VAR lag order
n.ahead <- 10    # GFEVD forecast horizon

last <- nrow(dat) - window - p + 1
cat("Rolling windows:", last, "\n")
cat("Ct dated at centre of each window; series spans 1981Q1 to 2016Q4\n\n")


# ---------------------------------------------------------------------------
# 3. Rolling Elastic-Net VAR + Diebold-Yilmaz spillover loop
# ---------------------------------------------------------------------------
matspill <- matrix(0, last, 1)

cat("Starting computation (est. 1-3 hours)...\n\n")
start_time <- proc.time()

for (i in seq_len(last)) {

  data1       <- dat[i:(window - 1 + p + i), ]
  Rolling_VAR <- LVAR(data1, p = p)
  rolling_fe  <- fevd_generalised(Rolling_VAR, n.ahead = n.ahead) * 100
  th0         <- sptable2012(rolling_fe)

  # Total connectedness Ct = bottom-right cell (Net row, From column)
  matspill[i] <- th0[(Rolling_VAR$K + 2), (Rolling_VAR$K + 1)]

  if (i %% 10 == 0 || i == 1 || i == last) {
    elapsed <- (proc.time() - start_time)[3]
    eta     <- if (i > 1) elapsed / i * (last - i) else NA
    cat(sprintf("Window %d / %d  |  Ct = %.3f  |  Elapsed: %.1f min  |  ETA: %.1f min\n",
                i, last, matspill[i], elapsed / 60,
                ifelse(is.na(eta), NA, eta / 60)))
  }
}

cat("\nDone! Total:", round((proc.time() - start_time)[3] / 60, 1), "min\n")


# ---------------------------------------------------------------------------
# 4. Save results
# ---------------------------------------------------------------------------
write.csv(matspill, file = "Ct_monthly.csv", row.names = TRUE)
cat("Saved: Ct_monthly.csv\n")
