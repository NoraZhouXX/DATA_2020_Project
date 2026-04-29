# =============================================================================
# Table_1.R  (lab version)
# Reproduce Table 1: Cumulative LP IRFs at h = 4, 8, 12, 16, 20
# Reference: Lee (2025), table_1.prg
#
# Lab methods used:
#   - lm()                  (Lab 4): OLS for each cumulative LP horizon
#   - sandwich::NeweyWest() (Lab 4/5): Newey-West HAC standard errors
#   - Dummy interactions     (Lab 9 DiD style): High vs Low connectedness
#   - Wald test via chi-sq   (Lab 4): test beta_H = beta_L
#
# INPUTS:  All_data.xls, Ct_dummies.csv
# OUTPUTS: Table_1.csv
# =============================================================================

rm(list = ls())

# install.packages(c("readxl","sandwich"))   # run once if needed
library(readxl)
library(sandwich)    # NeweyWest()

# ---------------------------------------------------------------------------
# 1. Load macro data and connectedness dummies
# ---------------------------------------------------------------------------
macro   <- read_excel("../data/All_data.xls", sheet = "Macro_data")
log_GDP <- log(as.numeric(macro$RGDP)) * 100
log_CPI <- log(as.numeric(macro$PCE))  * 100
FFR     <- as.numeric(macro$FFR)
n       <- nrow(macro)
yr      <- seq(1975.0, by = 0.25, length.out = n)

dum     <- read.csv("../data/Ct_dummies.csv", header = TRUE)
idx_dum <- match(round(dum$year, 4), round(yr, 4))
p1 <- rep(NA, n);  p1[idx_dum] <- dum$p1
n1 <- rep(NA, n);  n1[idx_dum] <- dum$n1

smpl   <- yr >= 1981.0 & yr <= 2007.75
T_smpl <- sum(smpl)
bw     <- floor(0.75 * T_smpl^(1/3))

cat("Macro data:", n, "quarters |",
    "Sample:", T_smpl, "quarters | NW bandwidth:", bw, "\n\n")

# ---------------------------------------------------------------------------
# 2. Build regression data frame (lags + interaction columns)
#    Same Lab 9 (DiD) interaction structure as Fig_2.R
# ---------------------------------------------------------------------------
lag_vec <- function(x, l) c(rep(NA, l), x[seq_len(n - l)])

make_reg_df <- function(cumY) {
  df <- data.frame(
    Y      = cumY,
    p1     = p1,
    n1     = n1,
    trend  = seq_len(n),
    trend2 = seq_len(n)^2
  )
  # Interaction: shock × regime (Lab 9 DiD style)
  df$p1_FFR <- p1 * FFR
  df$n1_FFR <- n1 * FFR
  # Control lags × regime
  for (l in 0:4) {
    df[[paste0("gdp_l", l)]]    <- lag_vec(log_GDP, l)
    df[[paste0("cpi_l", l)]]    <- lag_vec(log_CPI, l)
    df[[paste0("p1_gdp_l", l)]] <- p1 * df[[paste0("gdp_l", l)]]
    df[[paste0("n1_gdp_l", l)]] <- n1 * df[[paste0("gdp_l", l)]]
    df[[paste0("p1_cpi_l", l)]] <- p1 * df[[paste0("cpi_l", l)]]
    df[[paste0("n1_cpi_l", l)]] <- n1 * df[[paste0("cpi_l", l)]]
  }
  for (l in 1:4) {
    df[[paste0("ffr_l", l)]]    <- lag_vec(FFR, l)
    df[[paste0("p1_ffr_l", l)]] <- p1 * df[[paste0("ffr_l", l)]]
    df[[paste0("n1_ffr_l", l)]] <- n1 * df[[paste0("ffr_l", l)]]
  }
  df[smpl, ]
}

# ---------------------------------------------------------------------------
# 3. Run state-dependent cumulative LP at horizon h
#    Uses lm() (Lab 4) + NeweyWest() (sandwich)
# ---------------------------------------------------------------------------
run_lp_h <- function(cumY) {

  df <- make_reg_df(cumY)
  df <- df[complete.cases(df), ]

  # State-dependent specification: no intercept (p1 + n1 = 1 acts as intercept)
  # This is the Lab 9 DiD equivalent: interacting all controls with regime dummies
  sd_rhs <- c("0", "p1", "n1", "p1_FFR", "n1_FFR",
               paste0("p1_gdp_l", 0:4), paste0("n1_gdp_l", 0:4),
               paste0("p1_cpi_l", 0:4), paste0("n1_cpi_l", 0:4),
               paste0("p1_ffr_l", 1:4), paste0("n1_ffr_l", 1:4),
               "trend", "trend2")

  # lm() — Lab 4 method
  m <- lm(as.formula(paste("Y ~", paste(sd_rhs, collapse = " + "))), data = df)

  # Newey-West HAC SE — sandwich package (same spirit as heteroskedasticity-
  #   robust SE from Lab 4, extended for serial correlation)
  V  <- NeweyWest(m, lag = bw, prewhite = FALSE)
  se <- sqrt(diag(V))

  b_H <- -coef(m)["p1_FFR"]
  b_L <- -coef(m)["n1_FFR"]
  s_H <- se["p1_FFR"]
  s_L <- se["n1_FFR"]

  # Wald test: H0: beta_H = beta_L  (chi-sq(1), Lab 4 F-test equivalent)
  diff_var <- V["p1_FFR","p1_FFR"] + V["n1_FFR","n1_FFR"] -
              2 * V["p1_FFR","n1_FFR"]
  W        <- (coef(m)["p1_FFR"] - coef(m)["n1_FFR"])^2 / diff_var
  pval     <- 1 - pchisq(W, df = 1)

  list(coef_H    = as.numeric(b_H),
       coef_L    = as.numeric(b_L),
       se_H      = as.numeric(s_H),
       se_L      = as.numeric(s_L),
       pval_wald = as.numeric(pval))
}

# ---------------------------------------------------------------------------
# 4. Loop h = 1 to 20 with cumulative dependent variable
#    (cumulative LP: Y_t + Y_{t+1} + ... + Y_{t+h})
# ---------------------------------------------------------------------------
H_max   <- 20
cum_GDP <- log_GDP
cum_CPI <- log_CPI
cum_FFR <- FFR
results <- vector("list", H_max)

cat("Running cumulative local projections (h = 1 to", H_max, ")...\n")
for (h in 1:H_max) {
  # Accumulate h-period-ahead values
  cum_GDP <- cum_GDP + c(log_GDP[(h+1):n], rep(NA, h))
  cum_CPI <- cum_CPI + c(log_CPI[(h+1):n], rep(NA, h))
  cum_FFR <- cum_FFR + c(FFR[(h+1):n],     rep(NA, h))

  results[[h]] <- list(
    GDP = run_lp_h(cum_GDP),
    CPI = run_lp_h(cum_CPI),
    FFR = run_lp_h(cum_FFR)
  )

  cat(sprintf("h = %2d  |  GDP: High=%7.4f  Low=%7.4f  (Wald p=%.4f)\n",
              h, results[[h]]$GDP$coef_H,
                 results[[h]]$GDP$coef_L,
                 results[[h]]$GDP$pval_wald))
}

# ---------------------------------------------------------------------------
# 5. Build Table 1 at selected horizons
# ---------------------------------------------------------------------------
horizons <- c(4, 8, 12, 16, 20)

tab <- data.frame(
  Variable = c("GDP - High Connectedness", "GDP - Low Connectedness", "GDP - Wald p-value",
               "CPI - High Connectedness", "CPI - Low Connectedness", "CPI - Wald p-value",
               "FFR - High Connectedness", "FFR - Low Connectedness", "FFR - Wald p-value"),
  stringsAsFactors = FALSE
)

for (h in horizons) {
  r <- results[[h]]
  tab[[paste0("h=", h)]] <- round(c(
    r$GDP$coef_H, r$GDP$coef_L, r$GDP$pval_wald,
    r$CPI$coef_H, r$CPI$coef_L, r$CPI$pval_wald,
    r$FFR$coef_H, r$FFR$coef_L, r$FFR$pval_wald
  ), 4)
}

cat("\n========== TABLE 1 ==========\n")
print(tab, row.names = FALSE)

write.csv(tab, file = "Table_1.csv", row.names = FALSE)
cat("\nSaved: Table_1.csv\n")
