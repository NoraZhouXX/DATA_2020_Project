# =============================================================================
# extension/compute_directional.R
#
# L2 directional connectedness extension to Lee & Ma (2025).
# ---------------------------------------------------------------------------
# Diebold & Yilmaz (2014, Journal of Econometrics) define connectedness at
# three hierarchical levels:
#   L1 - pairwise directional        (51 x 51 = 2,550 measures per window)
#   L2 - total directional           (51 To + 51 From + 51 Net per window)
#   L3 - system-wide total           (1 scalar per window: this is what the
#                                     paper calls C_t = Ct_monthly.csv)
#
# Lee & Ma (2025) only use L3 in their main analysis. This extension computes
# and saves L2 measures so we can ask:
#   - Which states are persistent "net transmitters" of housing-market shocks?
#   - Which states are persistent "net receivers"?
#   - Do these roles change between high- and low-connectedness regimes?
#
# Cost: same 20 min as Compute_Ct.R (we re-run the rolling VAR loop;
# the underlying spillover table sptable2012() already computes To/From/Net
# at every window, but Lee & Ma's code only saves the bottom-right scalar).
# Consequence: outputs for L2 are exactly consistent with the published
# Ct_monthly.csv values (no re-estimation drift), and we additionally write
# matspill alongside as a sanity check.
#
# OUTPUTS (all written to ./outputs/):
#   directional_To_monthly.csv     wide: date + 51 state columns of To
#   directional_From_monthly.csv   wide: date + 51 state columns of From
#   directional_Net_monthly.csv    wide: date + 51 state columns of Net
#   directional_long.csv           long: date, state, To, From, Net
#                                  (22,032 rows = 432 windows x 51 states)
#   state_ranking_by_Net.csv       51 rows: per-state time-averaged
#                                  mean_To, mean_From, mean_Net, sd_Net,
#                                  pct_pos_Net, ranked by mean_Net
#
# USAGE: from RStudio, set working directory to this folder (extension/),
# then Source. Run time ~20 min; progress prints every 10 windows.
# =============================================================================

# ---------------------------------------------------------------------------
# 0. Sanity: confirm we are running from extension/ (relative paths assume so)
# ---------------------------------------------------------------------------
if (!file.exists("../scripts/connectedness_helpers.R")) {
  stop("Cannot find ../connectedness_helpers.R. ",
       "Set RStudio working directory to the extension/ folder before sourcing.")
}
if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
  cat("Created outputs/ directory.\n")
}

# ---------------------------------------------------------------------------
# 1. Load shared helpers + data
# ---------------------------------------------------------------------------
source("../scripts/connectedness_helpers.R")

raw        <- read.csv("../data/EI_DATA.csv", header = TRUE)
sym        <- raw[, c(4:60)]
data_House <- as.matrix(sym[, c(1:51)])

# 'data' is referenced globally by sptable2012() to label rows/columns.
# Required for the helper to attach state names; do not rename.
dat  <- data_House
data <- sym[, c(1:51)]

state_names <- colnames(data_House)   # AK, AL, ..., WY

cat("Data loaded:", nrow(dat), "monthly observations,", ncol(dat), "states\n")
cat("Sample: Year", raw$Year[1], "Month", raw$Month[1],
    "to Year", raw$Year[nrow(raw)], "Month", raw$Month[nrow(raw)], "\n\n")

# ---------------------------------------------------------------------------
# 2. Rolling-window parameters (must match Compute_Ct.R for L3 consistency)
# ---------------------------------------------------------------------------
window  <- 120   # 10-year rolling window (months)
p       <- 1     # VAR lag order
n.ahead <- 10    # GFEVD forecast horizon

last <- nrow(dat) - window - p + 1
cat("Rolling windows:", last, "\n")
cat("L2 vectors dated at centre of each window; series spans 1981m1 to 2016m12\n\n")

# Date vector matching the Lee & Ma (2025) convention used in Fig_1.R / Fig_2.R:
# 432 monthly observations starting January 1981.
date_vec <- seq(as.Date("1981-01-01"), by = "month", length.out = last)

# ---------------------------------------------------------------------------
# 3. Storage
#    - matspill : L3 scalar (= Ct, sanity-check vs Ct_monthly.csv)
#    - mat_to   : L2 To   per state per window (last x 51)
#    - mat_from : L2 From per state per window
#    - mat_net  : L2 Net  per state per window
# ---------------------------------------------------------------------------
matspill <- matrix(0, last, 1)
mat_to   <- matrix(0, last, length(state_names))
mat_from <- matrix(0, last, length(state_names))
mat_net  <- matrix(0, last, length(state_names))
colnames(mat_to)   <- state_names
colnames(mat_from) <- state_names
colnames(mat_net)  <- state_names

# ---------------------------------------------------------------------------
# 4. Rolling Elastic-Net VAR + Diebold-Yilmaz spillover loop
#    Identical computation to Compute_Ct.R; we just keep more of the output
#    of sptable2012() at each window.
# ---------------------------------------------------------------------------
cat("Starting computation (est. 1-3 hours)...\n\n")
start_time <- proc.time()

for (i in seq_len(last)) {

  data1       <- dat[i:(window - 1 + p + i), ]
  Rolling_VAR <- LVAR(data1, p = p)
  rolling_fe  <- fevd_generalised(Rolling_VAR, n.ahead = n.ahead) * 100
  th0         <- sptable2012(rolling_fe)

  K <- Rolling_VAR$K   # = 51

  # L3 scalar (bottom-right cell)
  matspill[i] <- th0[K + 2, K + 1]

  # L2 vectors
  #   th0[K + 1, 1:K] -> "To"   row (col-sums of off-diagonal of GFEVD)
  #   th0[1:K,   K+1] -> "From" col (row-sums of off-diagonal of GFEVD)
  #   th0[K + 2, 1:K] -> "Net"  row (= To - From per state)
  mat_to[i,   ] <- th0[K + 1, 1:K]
  mat_from[i, ] <- th0[1:K,   K + 1]
  mat_net[i,  ] <- th0[K + 2, 1:K]

  if (i %% 10 == 0 || i == 1 || i == last) {
    elapsed <- (proc.time() - start_time)[3]
    eta     <- if (i > 1) elapsed / i * (last - i) else NA
    cat(sprintf(
      "Window %d / %d  |  Ct = %.3f  |  Elapsed: %.1f min  |  ETA: %.1f min\n",
      i, last, matspill[i], elapsed / 60,
      ifelse(is.na(eta), NA, eta / 60)
    ))
  }
}

cat("\nDone! Total:",
    round((proc.time() - start_time)[3] / 60, 1), "min\n\n")

# ---------------------------------------------------------------------------
# 5. Sanity check: matspill should match existing Ct_monthly.csv exactly
# ---------------------------------------------------------------------------
existing_csv <- "../data/Ct_monthly.csv"
if (file.exists(existing_csv)) {
  existing_ct <- read.csv(existing_csv)[, 2]
  if (length(existing_ct) == last) {
    diffs <- abs(as.numeric(matspill) - existing_ct)
    cat(sprintf(
      "Sanity check vs Ct_monthly.csv:\n  max abs diff = %.3e  |  %s\n\n",
      max(diffs),
      if (max(diffs) < 1e-9) "PASS (bit-for-bit match)"
      else if (max(diffs) < 1e-3) "ok (within rounding tolerance)"
      else "FAIL: investigate"
    ))
  }
}

# ---------------------------------------------------------------------------
# 6. Save outputs (5 CSVs, all to outputs/)
# ---------------------------------------------------------------------------

# (a) Three wide tables: one per measure, rows = month, cols = state
write.csv(
  data.frame(date = date_vec, mat_to),
  file = "outputs/directional_To_monthly.csv", row.names = FALSE
)
write.csv(
  data.frame(date = date_vec, mat_from),
  file = "outputs/directional_From_monthly.csv", row.names = FALSE
)
write.csv(
  data.frame(date = date_vec, mat_net),
  file = "outputs/directional_Net_monthly.csv", row.names = FALSE
)

# (b) One long table: 22,032 rows, ready for ggplot / dplyr group-by
n_states <- length(state_names)
long_df <- data.frame(
  date  = rep(date_vec, n_states),
  state = rep(state_names, each = last),
  To    = as.vector(mat_to),
  From  = as.vector(mat_from),
  Net   = as.vector(mat_net),
  stringsAsFactors = FALSE
)
write.csv(long_df, file = "outputs/directional_long.csv", row.names = FALSE)

# (c) Ranking table: time-averaged statistics per state, sorted by mean_Net
ranking <- data.frame(
  state = state_names,
  mean_To     = colMeans(mat_to),
  mean_From   = colMeans(mat_from),
  mean_Net    = colMeans(mat_net),
  sd_Net      = apply(mat_net, 2, sd),
  pct_pos_Net = colMeans(mat_net > 0) * 100,   # % of months with Net > 0
  stringsAsFactors = FALSE
)
ranking <- ranking[order(-ranking$mean_Net), ]
rownames(ranking) <- NULL
write.csv(ranking, file = "outputs/state_ranking_by_Net.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 7. Summary report
# ---------------------------------------------------------------------------
cat("=== Files written to outputs/ ===\n")
cat("  directional_To_monthly.csv     ",
    nrow(mat_to),   "x", ncol(mat_to) + 1, " (date + 51 states)\n")
cat("  directional_From_monthly.csv   ",
    nrow(mat_from), "x", ncol(mat_from) + 1, "\n")
cat("  directional_Net_monthly.csv    ",
    nrow(mat_net),  "x", ncol(mat_net) + 1, "\n")
cat("  directional_long.csv           ",
    nrow(long_df),  "x", ncol(long_df), "\n")
cat("  state_ranking_by_Net.csv        51 x 6\n\n")

cat("=== Top 5 net transmitters (highest mean_Net) ===\n")
print(head(ranking, 5), row.names = FALSE)
cat("\n=== Top 5 net receivers (lowest mean_Net) ===\n")
print(tail(ranking, 5), row.names = FALSE)

# Mathematical sanity: at every t, sum(To) should equal sum(From) (network
# property), so sum(Net) should be ~0.
net_sum_per_t <- rowSums(mat_net)
cat(sprintf("\nNetwork sanity check: max |sum_t(Net)| = %.3e (should be ~0)\n",
            max(abs(net_sum_per_t))))

cat("\nDONE. Hand off the 5 CSVs to Cris (Person B) for visualization.\n")
