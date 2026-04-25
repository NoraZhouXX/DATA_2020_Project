# =============================================================================
# LP_Ct.R
# Local Projection (LP) function + Table 1 replication
#
# Reference: Lee (2025), "Housing market connectedness and transmission of
#            monetary policy", Economic Inquiry
# Adapted from ICPSR_216361-V3.1/Table_1/table_1.prg
#
# INPUTS:
#   All_data.xls   (Macro_data sheet: RGDP, PCE, FFR)
#   Ct_dummies.csv (from Fig_1.R: p1, n1 quarterly dummies)
#
# OUTPUT:
#   Table_1.csv
#
# No external packages required.
# =============================================================================

rm(list = ls())

library(readxl)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# newey_west_vcov: Newey-West HAC variance-covariance matrix
# Args:
#   X   : T x k regressor matrix
#   e   : T x 1 residual vector
#   bw  : bandwidth (number of lags); if NULL, uses floor(4*(T/100)^(2/9))
# Returns: k x k HAC covariance matrix of OLS coefficients
# -----------------------------------------------------------------------------
newey_west_vcov <- function(X, e, bw = NULL) {
  T <- nrow(X)
  k <- ncol(X)
  if (is.null(bw)) bw <- floor(4 * (T / 100)^(2/9))

  XtX_inv <- solve(t(X) %*% X)

  # Score matrix: S_t = x_t * e_t
  S <- X * as.numeric(e)   # T x k

  # Newey-West long-run variance: Omega = Gamma_0 + sum_l w_l*(Gamma_l + Gamma_l')
  Omega <- t(S) %*% S / T  # Gamma_0

  for (l in 1:bw) {
    w       <- 1 - l / (bw + 1)   # Bartlett kernel weight
    Gamma_l <- t(S[(l+1):T, , drop=FALSE]) %*% S[1:(T-l), , drop=FALSE] / T
    Omega   <- Omega + w * (Gamma_l + t(Gamma_l))
  }

  # V = (X'X)^{-1} * T * Omega * (X'X)^{-1}
  V <- XtX_inv %*% (T * Omega) %*% XtX_inv
  return(V)
}


# -----------------------------------------------------------------------------
# lp_run: Run one local projection regression at a given horizon h
#
# Model (following Lee 2025 / table_1.prg):
#   cumY_h(t) ~ p1(t)*shock(t) + n1(t)*shock(t)
#             + p1(t)*[gdp(t:t-4)] + p1(t)*[ffr(t-1:t-4)] + p1(t)*[cpi(t:t-4)]
#             + n1(t)*[gdp(t:t-4)] + n1(t)*[ffr(t-1:t-4)] + n1(t)*[cpi(t:t-4)]
#             + trend + trend^2
#
# Note: coef(1) = beta_H (high connectedness), coef(2) = beta_L (low)
#       Both multiplied by -1 so a positive FFR shock → positive reported response
#
# Args:
#   cumY  : T-vector, cumulative dependent variable at horizon h
#   shock : T-vector, monetary policy shock (FFR)
#   p1    : T-vector, high-connectedness dummy (lagged 1 quarter)
#   n1    : T-vector, low-connectedness dummy (lagged 1 quarter)
#   gdp   : T-vector, log GDP
#   cpi   : T-vector, log CPI
#   ffr   : T-vector, FFR level
#   smpl  : logical index vector selecting the estimation sample
#   bw    : Newey-West bandwidth (NULL = automatic)
#
# Returns: list with coef_H, coef_L, se_H, se_L, pval_wald
# -----------------------------------------------------------------------------
lp_run <- function(cumY, shock, p1, n1, gdp, cpi, ffr, smpl, bw = NULL) {

  T_full <- length(cumY)

  # Build regressors (full-length vectors, then subset)
  build_lags <- function(x, lags) {
    # lags: integer vector of lag orders (e.g., 0:4 or 1:4)
    sapply(lags, function(l) {
      c(rep(NA, l), x[1:(T_full - l)])
    })
  }

  gdp_lags <- build_lags(gdp, 0:4)   # 5 columns: t, t-1, ..., t-4
  cpi_lags <- build_lags(cpi, 0:4)   # 5 columns
  ffr_lags <- build_lags(ffr, 1:4)   # 4 columns: t-1, ..., t-4
  trend    <- 1:T_full

  # Interaction with state dummies
  X_raw <- cbind(
    p1 * shock,          # coef 1: beta_H (high connectedness × shock)
    n1 * shock,          # coef 2: beta_L (low connectedness × shock)
    p1 * gdp_lags,       # 5 cols
    p1 * ffr_lags,       # 4 cols
    p1 * cpi_lags,       # 5 cols
    n1 * gdp_lags,       # 5 cols
    n1 * ffr_lags,       # 4 cols
    n1 * cpi_lags,       # 5 cols
    trend,               # 1 col
    trend^2              # 1 col
  )                      # Total: 2 + 5+4+5 + 5+4+5 + 2 = 32 columns

  # Apply sample restriction and drop NAs
  X_s <- X_raw[smpl, , drop = FALSE]
  Y_s <- cumY[smpl]
  ok  <- complete.cases(X_s, Y_s)
  X_s <- X_s[ok, , drop = FALSE]
  Y_s <- Y_s[ok]

  # OLS
  beta <- solve(t(X_s) %*% X_s) %*% t(X_s) %*% Y_s
  e    <- Y_s - X_s %*% beta

  # HAC standard errors
  V    <- newey_west_vcov(X_s, e, bw = bw)
  se   <- sqrt(diag(V))

  # Wald test: H0: beta_H = beta_L  (i.e., coef[1] = coef[2])
  # W = (beta_H - beta_L)^2 / Var(beta_H - beta_L)
  # Var(beta_H - beta_L) = V[1,1] + V[2,2] - 2*V[1,2]
  diff_var <- V[1,1] + V[2,2] - 2 * V[1,2]
  W        <- (beta[1] - beta[2])^2 / diff_var
  pval     <- 1 - pchisq(W, df = 1)

  list(
    coef_H    = as.numeric(beta[1]) * -1,   # sign flip (see paper)
    coef_L    = as.numeric(beta[2]) * -1,
    se_H      = se[1],
    se_L      = se[2],
    pval_wald = as.numeric(pval)
  )
}


# =============================================================================
# LOAD DATA
# =============================================================================

# --- Macro data (quarterly, 1975Q1 - 2020Q4) ---------------------------------
macro <- read_excel("All_data.xls", sheet = "Macro_data")

log_GDP <- log(as.numeric(macro$RGDP)) * 100
log_CPI <- log(as.numeric(macro$PCE))  * 100
FFR     <- as.numeric(macro$FFR)

# Quarterly date vector (decimal year, 1975Q1 = 1975.0)
n_macro <- nrow(macro)
yr_macro <- seq(1975.0, by = 0.25, length.out = n_macro)

# --- Connectedness dummies (quarterly, 1981Q1 - 2016Q4) ----------------------
dum  <- read.csv("Ct_dummies.csv", header = TRUE)
yr_d <- dum$year   # 1981.0, 1981.25, ..., 2016.75

# Align dummies with macro data by matching year values
idx_dum_in_macro <- match(round(yr_d, 4), round(yr_macro, 4))

p1_full <- rep(NA, n_macro)
n1_full <- rep(NA, n_macro)
p1_full[idx_dum_in_macro] <- dum$p1
n1_full[idx_dum_in_macro] <- dum$n1

shock_full <- FFR   # monetary policy shock = FFR level

# --- Estimation sample: 1981Q1 to 2007Q4 ------------------------------------
smpl <- yr_macro >= 1981.0 & yr_macro <= 2007.75

cat("Macro data:", n_macro, "quarters (", yr_macro[1], "-", tail(yr_macro,1), ")\n")
cat("Estimation sample:", sum(smpl), "quarters (1981Q1 - 2007Q4)\n\n")


# =============================================================================
# RUN LOCAL PROJECTIONS: h = 1 to 20
# =============================================================================

H_max  <- 20
results <- vector("list", H_max)

# Cumulative dependent variables (updated each iteration, like pc_gdp in EViews)
cum_GDP <- log_GDP
cum_CPI <- log_CPI
cum_FFR <- FFR

cat("Running local projections (h = 1 to", H_max, ")...\n")

for (h in 1:H_max) {

  # Cumulate by adding h-period ahead value
  cum_GDP <- cum_GDP + c(log_GDP[(h+1):n_macro], rep(NA, h))
  cum_CPI <- cum_CPI + c(log_CPI[(h+1):n_macro], rep(NA, h))
  cum_FFR <- cum_FFR + c(FFR[(h+1):n_macro],     rep(NA, h))

  # Newey-West bandwidth: fixed EViews default = floor(0.75 * T^(1/3))
  T_smpl <- sum(smpl)
  bw_h   <- floor(0.75 * T_smpl^(1/3))

  res_GDP <- lp_run(cum_GDP, shock_full, p1_full, n1_full,
                    log_GDP, log_CPI, FFR, smpl, bw = bw_h)
  res_CPI <- lp_run(cum_CPI, shock_full, p1_full, n1_full,
                    log_GDP, log_CPI, FFR, smpl, bw = bw_h)
  res_FFR <- lp_run(cum_FFR, shock_full, p1_full, n1_full,
                    log_GDP, log_CPI, FFR, smpl, bw = bw_h)

  results[[h]] <- list(GDP = res_GDP, CPI = res_CPI, FFR = res_FFR)
  cat(sprintf("h = %2d  |  GDP: High=%.3f Low=%.3f (p=%.4f)\n",
              h, res_GDP$coef_H, res_GDP$coef_L, res_GDP$pval_wald))
}


# =============================================================================
# BUILD TABLE 1 (horizons h = 4, 8, 12, 16, 20)
# =============================================================================

horizons <- c(4, 8, 12, 16, 20)

tab <- data.frame(
  Variable  = c("GDP - High Connectedness", "GDP - Low Connectedness", "GDP - P-Value",
                "CPI - High Connectedness", "CPI - Low Connectedness", "CPI - P-Value",
                "FFR - High Connectedness", "FFR - Low Connectedness", "FFR - P-Value"),
  stringsAsFactors = FALSE
)

for (h in horizons) {
  col <- sprintf("h=%d", h)
  r   <- results[[h]]
  tab[[col]] <- c(
    round(r$GDP$coef_H,    6), round(r$GDP$coef_L,    6), round(r$GDP$pval_wald, 4),
    round(r$CPI$coef_H,    6), round(r$CPI$coef_L,    6), round(r$CPI$pval_wald, 4),
    round(r$FFR$coef_H,    6), round(r$FFR$coef_L,    6), round(r$FFR$pval_wald, 4)
  )
}

cat("\n========== TABLE 1 ==========\n")
print(tab, row.names = FALSE)

write.csv(tab, file = "Table_1.csv", row.names = FALSE)
cat("\nSaved: Table_1.csv\n")
