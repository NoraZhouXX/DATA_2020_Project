# =============================================================================
# Compute_Ct_macro.R
# Diebold-Yilmaz total connectedness index WITH macroeconomic controls.
#
# Method:
#   - Rolling window (10y = 120 months), step 1 month.
#   - Within each window: Elastic-Net VAR(1) on the FULL system (51 housing
#     returns + macro variables in EI_DATA.csv columns 55+).
#   - Compute generalised FEVD on the full system (K_total x K_total).
#   - Extract the K_house x K_house housing sub-block, re-normalise rows to
#     sum to 100, and compute Ct from the off-diagonal elements only.
#
# Economic interpretation:
#   The macro variables enter the VAR as controls (they help estimate cleaner
#   between-state spillover coefficients), but the connectedness index itself
#   is computed only on the housing block, so the index measures "state-to-
#   state housing spillover net of common macro drivers."
#
# All numerical helpers live in connectedness_helpers.R (sourced below).
#
# INPUT:  EI_DATA.csv  (51 housing + macro variables in columns 4:end)
# OUTPUT: Ct_macro_monthly.csv
# =============================================================================

rm(list = ls())

# Load all helper functions. sptable2012() is NOT used in this script (we
# extract a sub-block manually), but the rest of the pipeline is identical.
source("connectedness_helpers.R")


# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
raw <- read.csv("../data/EI_DATA.csv", header = TRUE)
sym <- raw[, c(4:ncol(raw))]   # all variable columns (housing + macro)

K_house <- 51                  # number of housing series (always first 51)
dat_all <- as.matrix(sym)      # full system (housing + macro)
K_total <- ncol(dat_all)

cat("Data loaded:", nrow(dat_all), "monthly observations |",
    K_house, "housing +", K_total - K_house, "macro variables\n\n")


# ---------------------------------------------------------------------------
# 2. Rolling window parameters (identical to housing-only benchmark)
# ---------------------------------------------------------------------------
window  <- 120
p       <- 1
n.ahead <- 10

last <- nrow(dat_all) - window - p + 1
cat("Rolling windows:", last, "\n\n")


# ---------------------------------------------------------------------------
# 3. Rolling Elastic-Net VAR + housing sub-block connectedness loop
# ---------------------------------------------------------------------------
matspill <- matrix(0, last, 1)

cat("Starting computation with macro controls...\n\n")
start_time <- proc.time()

for (i in seq_len(last)) {

  data1       <- dat_all[i:(window - 1 + p + i), ]
  Rolling_VAR <- LVAR(data1, p = p)

  # Full-system GFEVD (K_total x K_total, row-normalised, scaled by 100)
  fevd_full <- fevd_generalised(Rolling_VAR, n.ahead = n.ahead) * 100

  # Extract housing sub-block (first K_house rows/cols) and re-normalise
  sp_sub <- fevd_full[1:K_house, 1:K_house]
  sp_sub <- sp_sub / rowSums(sp_sub) * 100

  # Total connectedness on housing sub-block (DY 2012 definition)
  matspill[i] <- (sum(sp_sub) - sum(diag(sp_sub))) / sum(sp_sub) * 100

  if (i %% 10 == 0 || i == 1 || i == last) {
    elapsed <- (proc.time() - start_time)[3]
    eta     <- if (i > 1) elapsed / i * (last - i) else NA
    cat(sprintf("Window %d / %d  |  Ct_macro = %.3f  |  Elapsed: %.1f min  |  ETA: %.1f min\n",
                i, last, matspill[i], elapsed / 60,
                ifelse(is.na(eta), NA, eta / 60)))
  }
}

cat("\nDone! Total:", round((proc.time() - start_time)[3] / 60, 1), "min\n")


# ---------------------------------------------------------------------------
# 4. Save results
# ---------------------------------------------------------------------------
write.csv(matspill, file = "../data/Ct_macro_monthly.csv", row.names = TRUE)
cat("Saved: Ct_macro_monthly.csv\n")
