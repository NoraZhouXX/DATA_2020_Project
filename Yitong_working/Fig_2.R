# =============================================================================
# Fig_2.R
# Reproduce Figure 2: Baseline IRFs (non-cumulative local projections)
#
# For each horizon h = 0 to 20, regresses log_GDP(t+h), log_CPI(t+h),
# FFR(t+h) on the monetary shock, with and without state-dependent dummies.
# Plots 3 rows (GDP / CPI / FFR) x 3 cols (Three Models / Linear / State-Dep)
#
# Reference: Lee (2025), adapted from ICPSR_216361-V3.1/Fig_2/Fig_2.prg + .m
#
# INPUTS:  All_data.xls, Ct_dummies.csv
# OUTPUT:  Fig_2.png, Fig_2_GDP.png, Fig_2_CPI.png, Fig_2_FFR.png
#
# No external packages required.
# =============================================================================

rm(list = ls())
library(readxl)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

newey_west_vcov <- function(X, e, bw = NULL) {
  T <- nrow(X)
  if (is.null(bw)) bw <- floor(0.75 * T^(1/3))
  XtX_inv <- solve(t(X) %*% X)
  S <- X * as.numeric(e)
  Omega <- t(S) %*% S / T
  for (l in 1:bw) {
    w       <- 1 - l / (bw + 1)
    Gamma_l <- t(S[(l+1):T, , drop=FALSE]) %*% S[1:(T-l), , drop=FALSE] / T
    Omega   <- Omega + w * (Gamma_l + t(Gamma_l))
  }
  XtX_inv %*% (T * Omega) %*% XtX_inv
}

# Non-cumulative LP: regress Y(t+h) on shock and controls
# Returns: list(coef, se) for specified coefficient indices
lp_point <- function(Y_lead, shock, p1, n1, gdp, cpi, ffr,
                     smpl, type = c("linear", "state"), bw = NULL) {
  type <- match.arg(type)
  T_full <- length(Y_lead)

  build_lags <- function(x, lags) {
    sapply(lags, function(l) c(rep(NA, l), x[1:(T_full - l)]))
  }

  gdp_lags <- build_lags(gdp, 0:4)
  cpi_lags <- build_lags(cpi, 0:4)
  ffr_lags <- build_lags(ffr, 1:4)
  trend    <- 1:T_full

  if (type == "linear") {
    # Linear model: intercept + shock + ffr lags + gdp lags + cpi lags + trend
    X_raw <- cbind(1, shock,
                   ffr_lags, gdp_lags, cpi_lags,
                   trend, trend^2)
    shock_idx <- 2   # coefficient on shock (FFR)
  } else {
    # State-dependent model (no intercept; p1+n1 = 1 acts as intercept)
    X_raw <- cbind(
      p1 * shock, n1 * shock,
      p1 * gdp_lags, p1 * ffr_lags, p1 * cpi_lags,
      n1 * gdp_lags, n1 * ffr_lags, n1 * cpi_lags,
      trend, trend^2
    )
    shock_idx <- 1:2   # coef 1 = High, coef 2 = Low
  }

  X_s <- X_raw[smpl, , drop = FALSE]
  Y_s <- Y_lead[smpl]
  ok  <- complete.cases(X_s, Y_s)
  X_s <- X_s[ok, , drop = FALSE]
  Y_s <- Y_s[ok]

  beta <- solve(t(X_s) %*% X_s) %*% t(X_s) %*% Y_s
  e    <- Y_s - X_s %*% beta
  V    <- newey_west_vcov(X_s, e, bw = bw)
  se   <- sqrt(diag(V))

  list(coef = as.numeric(beta[shock_idx]) * -1,
       se   = se[shock_idx])
}

# =============================================================================
# LOAD DATA
# =============================================================================
macro   <- read_excel("All_data.xls", sheet = "Macro_data")
log_GDP <- log(as.numeric(macro$RGDP)) * 100
log_CPI <- log(as.numeric(macro$PCE))  * 100
FFR     <- as.numeric(macro$FFR)
n_macro <- nrow(macro)
yr_macro <- seq(1975.0, by = 0.25, length.out = n_macro)

dum  <- read.csv("Ct_dummies.csv", header = TRUE)
yr_d <- dum$year
idx_dum <- match(round(yr_d, 4), round(yr_macro, 4))
p1_full <- rep(NA, n_macro); p1_full[idx_dum] <- dum$p1
n1_full <- rep(NA, n_macro); n1_full[idx_dum] <- dum$n1

smpl <- yr_macro >= 1981.0 & yr_macro <= 2007.75
T_smpl <- sum(smpl)
bw_fix <- floor(0.75 * T_smpl^(1/3))

cat("Sample:", T_smpl, "quarters\n")
cat("NW bandwidth:", bw_fix, "\n\n")

# =============================================================================
# RUN LP: h = 0 to 20
# =============================================================================
H_max <- 20
h_seq <- 0:H_max
nh    <- length(h_seq)

# Storage: [point, lower_CI, upper_CI] for each variable x model
res <- list(
  GDP = list(lin = matrix(NA, nh, 3), high = matrix(NA, nh, 3),
             low  = matrix(NA, nh, 3)),
  CPI = list(lin = matrix(NA, nh, 3), high = matrix(NA, nh, 3),
             low  = matrix(NA, nh, 3)),
  FFR = list(lin = matrix(NA, nh, 3), high = matrix(NA, nh, 3),
             low  = matrix(NA, nh, 3))
)

vars <- list(GDP = log_GDP, CPI = log_CPI, FFR = FFR)

cat("Running LP for h = 0 to", H_max, "...\n")
for (hi in seq_along(h_seq)) {
  h <- h_seq[hi]

  for (vname in c("GDP", "CPI", "FFR")) {
    y    <- vars[[vname]]
    # Lead: y(t+h) at position t
    Y_lead <- c(y[(h+1):n_macro], rep(NA, h))

    # Linear model
    r_lin <- lp_point(Y_lead, FFR, p1_full, n1_full,
                      log_GDP, log_CPI, FFR, smpl,
                      type = "linear", bw = bw_fix)
    res[[vname]]$lin[hi, ] <- c(r_lin$coef - r_lin$se,
                                 r_lin$coef,
                                 r_lin$coef + r_lin$se)

    # State-dependent model
    r_sd <- lp_point(Y_lead, FFR, p1_full, n1_full,
                     log_GDP, log_CPI, FFR, smpl,
                     type = "state", bw = bw_fix)
    # High: coef[1], Low: coef[2]
    res[[vname]]$high[hi, ] <- c(r_sd$coef[1] - r_sd$se[1],
                                  r_sd$coef[1],
                                  r_sd$coef[1] + r_sd$se[1])
    res[[vname]]$low[hi, ]  <- c(r_sd$coef[2] - r_sd$se[2],
                                  r_sd$coef[2],
                                  r_sd$coef[2] + r_sd$se[2])
  }
  cat(sprintf("h=%2d  GDP: Lin=%.3f H=%.3f L=%.3f\n",
              h, res$GDP$lin[hi,2], res$GDP$high[hi,2], res$GDP$low[hi,2]))
}

# =============================================================================
# PLOT FUNCTION (one row = one variable)
# =============================================================================
plot_row <- function(lin, high, low, h_seq,
                     ylim1, ylim2, ylim3, title_var) {

  zero <- rep(0, length(h_seq))

  # --- Panel 1: Three Models (point estimates only) ---
  plot(h_seq, lin[,2], type = "o", pch = 16, col = "black", lwd = 1.5,
       xlim = c(-0.1, 20), ylim = ylim1,
       xlab = "Horizon (Quarter)", ylab = "",
       xaxt = "n", bty = "l", main = "Three Models")
  axis(1, at = c(0, 5, 10, 15, 20))
  lines(h_seq, high[,2], type = "o", pch = 2,  col = "blue", lwd = 1.5)
  lines(h_seq, low[,2],  type = "o", pch = 0,  col = "red",  lwd = 1.5)
  abline(h = 0, col = "black", lwd = 0.5)
  grid()
  if (title_var == "GDP")
    legend("topleft", legend = c("Linear","High","Low"),
           col = c("black","blue","red"), pch = c(16,2,0),
           lwd = 1.5, bty = "n", cex = 0.85)

  # --- Panel 2: Linear model with ±1 SE band ---
  plot(h_seq, lin[,2], type = "n",
       xlim = c(-0.1, 20), ylim = ylim2,
       xlab = "Horizon (Quarter)", ylab = "",
       xaxt = "n", bty = "l", main = "Linear Model")
  axis(1, at = c(0, 5, 10, 15, 20))
  polygon(c(h_seq, rev(h_seq)),
          c(lin[,3], rev(lin[,1])),
          col = rgb(0.85,0.85,0.85), border = NA)
  lines(h_seq, lin[,2], type = "o", pch = 16, col = "black", lwd = 1.5)
  lines(h_seq, lin[,1], lty = 2, col = "black", lwd = 1)
  lines(h_seq, lin[,3], lty = 2, col = "black", lwd = 1)
  abline(h = 0, col = "black", lwd = 0.5)
  grid()

  # --- Panel 3: State-Dependent (High=blue, Low=red) ---
  plot(h_seq, high[,2], type = "n",
       xlim = c(-0.1, 20), ylim = ylim3,
       xlab = "Horizon (Quarter)", ylab = "",
       xaxt = "n", bty = "l", main = "State-Dependent Models")
  axis(1, at = c(0, 5, 10, 15, 20))
  # High CI band (light blue)
  polygon(c(h_seq, rev(h_seq)),
          c(high[,3], rev(high[,1])),
          col = rgb(0.85,0.85,1.0), border = NA)
  lines(h_seq, high[,2], type = "o", pch = 2, col = "blue", lwd = 1.5)
  lines(h_seq, high[,1], lty = 2, col = "blue", lwd = 1)
  lines(h_seq, high[,3], lty = 2, col = "blue", lwd = 1)
  # Low (red dotted)
  lines(h_seq, low[,2], type = "o", pch = 0, col = "red", lwd = 1.5)
  lines(h_seq, low[,1], lty = 3, col = "red", lwd = 2)
  lines(h_seq, low[,3], lty = 3, col = "red", lwd = 2)
  abline(h = 0, col = "black", lwd = 0.5)
  grid()
}

# =============================================================================
# SAVE FIGURES
# =============================================================================

# --- Full 3x3 figure (one PNG per row, matching paper layout) ---

# GDP row
png("Fig_2_GDP.png", width = 12, height = 4, units = "in", res = 200)
par(mfrow = c(1,3), mar = c(4,3,3,1))
plot_row(res$GDP$lin, res$GDP$high, res$GDP$low, h_seq,
         ylim1 = c(-0.5, 1.5), ylim2 = range(res$GDP$lin) + c(-0.2, 0.2),
         ylim3 = range(c(res$GDP$high, res$GDP$low)) + c(-0.2, 0.2),
         title_var = "GDP")
mtext("GDP", side = 3, line = -1.5, outer = TRUE, cex = 1.2, font = 2)
dev.off()

# CPI row
png("Fig_2_CPI.png", width = 12, height = 4, units = "in", res = 200)
par(mfrow = c(1,3), mar = c(4,3,3,1))
plot_row(res$CPI$lin, res$CPI$high, res$CPI$low, h_seq,
         ylim1 = range(c(res$CPI$lin[,2], res$CPI$high[,2], res$CPI$low[,2])) + c(-0.1,0.1),
         ylim2 = c(-0.5, 1.0),
         ylim3 = c(-0.5, 2.0),
         title_var = "CPI")
mtext("PCE Deflator", side = 3, line = -1.5, outer = TRUE, cex = 1.2, font = 2)
dev.off()

# FFR row
png("Fig_2_FFR.png", width = 12, height = 4, units = "in", res = 200)
par(mfrow = c(1,3), mar = c(4,3,3,1))
plot_row(res$FFR$lin, res$FFR$high, res$FFR$low, h_seq,
         ylim1 = range(c(res$FFR$lin[,2], res$FFR$high[,2], res$FFR$low[,2])) + c(-0.2,0.2),
         ylim2 = range(res$FFR$lin) + c(-0.2, 0.2),
         ylim3 = c(-2, 1),
         title_var = "FFR")
mtext("FFR", side = 3, line = -1.5, outer = TRUE, cex = 1.2, font = 2)
dev.off()

# Combined 3x3 figure
png("Fig_2.png", width = 12, height = 10, units = "in", res = 200)
par(mfrow = c(3,3), mar = c(4,3,3,1), oma = c(0,2,0,0))
# Row 1: GDP
plot_row(res$GDP$lin, res$GDP$high, res$GDP$low, h_seq,
         ylim1 = c(-0.5, 1.5),
         ylim2 = range(res$GDP$lin) + c(-0.2, 0.2),
         ylim3 = range(c(res$GDP$high, res$GDP$low)) + c(-0.2, 0.2),
         title_var = "GDP")
mtext("GDP", side = 2, line = 0.5, at = par("usr")[4]*0.5, cex = 0.9, font = 2, las = 0)

# Row 2: CPI
plot_row(res$CPI$lin, res$CPI$high, res$CPI$low, h_seq,
         ylim1 = range(c(res$CPI$lin[,2], res$CPI$high[,2], res$CPI$low[,2])) + c(-0.1,0.1),
         ylim2 = c(-0.5, 1.0),
         ylim3 = c(-0.5, 2.0),
         title_var = "CPI")

# Row 3: FFR
plot_row(res$FFR$lin, res$FFR$high, res$FFR$low, h_seq,
         ylim1 = range(c(res$FFR$lin[,2], res$FFR$high[,2], res$FFR$low[,2])) + c(-0.2,0.2),
         ylim2 = range(res$FFR$lin) + c(-0.2, 0.2),
         ylim3 = c(-2, 1),
         title_var = "FFR")
dev.off()

cat("Saved: Fig_2.png, Fig_2_GDP.png, Fig_2_CPI.png, Fig_2_FFR.png\n")
