# ============================================================
# Fig_2.prg - R Translation (v2)
# Replication of Figure 2: State-Dependent Local Projections
# Variables: GDP, PCE Deflator, FFR
# ============================================================

# Install packages if needed
# install.packages(c("readxl", "mFilter", "sandwich", "lmtest"))

library(readxl)
library(mFilter)
library(sandwich)
library(lmtest)

# ============================================================
# STEP 1 - Load Data
# ============================================================

# Load macro data (adjust path if needed) 
# The macroeconomic variables
macro <- read_excel("~/Documents/Brown/Classes/DATA2020/Final Project/DATA_2020_Project/Figure 2/Fig_2_R/All_data.xls", sheet = "Macro_data") 

# Load connectedness index (monthly)
# The connectedness index
data2 <- read.csv("data2.csv")

# ============================================================
# STEP 2 - Convert connectedness from monthly to quarterly
# ============================================================

ci_monthly <- data2$V1  # 432 monthly values

T_q <- length(ci_monthly) / 3  # = 144 quarters
ci_quarterly <- numeric(T_q)
for (j in 1:T_q) {
  ci_quarterly[j] <- mean(ci_monthly[(3*(j-1)+1):(3*j)])
}

# Assign dates 1981Q1 to 2016Q4 (144 quarters, following data1.xlsx)

ci_dates <- seq(from = as.Date("1980-10-01"),
                by   = "quarter",
                length.out = T_q)

ci_df <- data.frame(date_q   = ci_dates,
                    c_index  = ci_quarterly)

# ============================================================
# STEP 3 - Prepare macro data
# ============================================================

# Parse dates
macro$date_q <- as.Date(paste0(
  substr(macro$date, 1, 4), "-",
  (as.numeric(substr(macro$date, 6, 6)) - 1) * 3 + 1, "-01"
))

# Create log variables (x100 to get percentages)
# Taking the log of GDP and PCE means the regression coefficients can be interpreted as percentage changes. 
macro$log_GDP <- log(macro$RGDP) * 100
macro$log_cpi <- log(macro$PCE)  * 100
macro$shock   <- macro$FFR

# ============================================================
# STEP 4 - Apply HP filter to FULL connectedness series first
# ============================================================

# Apply HP filter on raw values across full 144 quarters
# Original series = Trend + Cycle
# If C is above its trend → high connectedness period
# If C is below its trend → low connectedness period
log_ci <- log(ci_quarterly / 100) * 100
hp_result <- hpfilter(log_ci, freq = 1600)

ci_df$c_trend <- hp_result$trend
ci_df$c_cycle <- hp_result$cycle
ci_df$c_cycle <- -ci_df$c_cycle

# ============================================================
# STEP 5 - Create dummy variables on full series BEFORE merging
# ============================================================

# High connectedness: cycle >= 0
# Low connectedness:  cycle < 0
ci_df$pshock <- ifelse(ci_df$c_cycle >= 0, 1, 0)
ci_df$nshock <- ifelse(ci_df$c_cycle <  0, 1, 0)

# Lag by 1 period before merging (D_H,t-1 and D_L,t-1)
# We cannot use the current period's state in the regression — 
# that would create a simultaneity problem (using today's information to explain today's outcome)
n_ci    <- nrow(ci_df)
ci_df$p1 <- c(NA, ci_df$pshock[-n_ci])
ci_df$n1 <- c(NA, ci_df$nshock[-n_ci])


# ============================================================
# STEP 6 - Merge macro and connectedness data
# ============================================================

df <- merge(macro, ci_df, by = "date_q", all.x = TRUE)
df <- df[order(df$date_q), ]
rownames(df) <- NULL


# ---------- Step 7 (NEW): build lags on the full df ----------
make_lags <- function(x, max_lag) {
  n <- length(x)
  out <- matrix(NA, nrow = n, ncol = max_lag)
  for (i in 1:max_lag) out[(i+1):n, i] <- x[1:(n-i)]
  out
}

gdp_lags <- make_lags(df$log_GDP, 4); colnames(gdp_lags) <- paste0("gdp_lag", 1:4)
cpi_lags <- make_lags(df$log_cpi, 4); colnames(cpi_lags) <- paste0("cpi_lag", 1:4)
ffr_lags <- make_lags(df$shock,   4); colnames(ffr_lags) <- paste0("ffr_lag", 1:4)

df <- cbind(df, gdp_lags, cpi_lags, ffr_lags)

# Time trend on the FULL df, but indexed so trend = 1 at 1981Q1
# (matches EViews @trend behavior up to an affine shift, which
#  is absorbed by trend^2 + intercept and does not affect the
#  shock coefficient).
start_date <- as.Date("1981-01-01")
end_date   <- as.Date("2007-10-01")

df$trend  <- seq_len(nrow(df)) - which(df$date_q == start_date) + 1
df$trend2 <- df$trend^2

# Mark which rows belong to the 1981Q1-2007Q4 estimation window
df$in_sample <- df$date_q >= start_date & df$date_q <= end_date
N_full <- nrow(df)

# ---------- Step 8 (NEW): regression loop ----------
H    <- 21
conf <- 1

est_1 <- matrix(NA, H, 12); est_2 <- matrix(NA, H, 12); est_3 <- matrix(NA, H, 12)
cn <- c("point_linear","point_high","point_low",
        "linear1","linear2","linear3",
        "H1","H2","H3","L1","L2","L3")
colnames(est_1) <- colnames(est_2) <- colnames(est_3) <- cn

ctrl_lin <- "shock +
             ffr_lag1 + ffr_lag2 + ffr_lag3 + ffr_lag4 +
             log_GDP + gdp_lag1 + gdp_lag2 + gdp_lag3 + gdp_lag4 +
             log_cpi + cpi_lag1 + cpi_lag2 + cpi_lag3 + cpi_lag4 +
             trend + trend2"

sd_formula <- "
    I(p1 * shock) + I(n1 * shock) +
    p1 +
    I(p1 * log_GDP)  + I(p1 * gdp_lag1) + I(p1 * gdp_lag2) +
    I(p1 * gdp_lag3) + I(p1 * gdp_lag4) +
    I(p1 * ffr_lag1) + I(p1 * ffr_lag2) +
    I(p1 * ffr_lag3) + I(p1 * ffr_lag4) +
    I(p1 * log_cpi)  + I(p1 * cpi_lag1) + I(p1 * cpi_lag2) +
    I(p1 * cpi_lag3) + I(p1 * cpi_lag4) +
    n1 +
    I(n1 * log_GDP)  + I(n1 * gdp_lag1) + I(n1 * gdp_lag2) +
    I(n1 * gdp_lag3) + I(n1 * gdp_lag4) +
    I(n1 * ffr_lag1) + I(n1 * ffr_lag2) +
    I(n1 * ffr_lag3) + I(n1 * ffr_lag4) +
    I(n1 * log_cpi)  + I(n1 * cpi_lag1) + I(n1 * cpi_lag2) +
    I(n1 * cpi_lag3) + I(n1 * cpi_lag4) +
    trend + trend2 - 1"

for (h in 0:20) {
  
  # Leads built on the FULL df (so they can extend past 2007Q4)
  if (h > 0) {
    gdp_lead <- c(df$log_GDP[(h+1):N_full], rep(NA, h))
    cpi_lead <- c(df$log_cpi[(h+1):N_full], rep(NA, h))
    ffr_lead <- c(df$shock  [(h+1):N_full], rep(NA, h))
  } else {
    gdp_lead <- df$log_GDP
    cpi_lead <- df$log_cpi
    ffr_lead <- df$shock
  }
  
  # Now restrict to in-sample rows AND drop only rows whose
  # required regressors / dep var are missing.
  use <- df$in_sample & complete.cases(
    cbind(gdp_lead, cpi_lead, ffr_lead,
          df[, c("shock","p1","n1","trend","trend2",
                 "log_GDP","gdp_lag1","gdp_lag2","gdp_lag3","gdp_lag4",
                 "log_cpi","cpi_lag1","cpi_lag2","cpi_lag3","cpi_lag4",
                 "ffr_lag1","ffr_lag2","ffr_lag3","ffr_lag4")])
  )
  
  temp <- df[use, ]
  gdp_y <- gdp_lead[use]; cpi_y <- cpi_lead[use]; ffr_y <- ffr_lead[use]
  
  if (h <= 2) cat("h =", h, " | nobs =", nrow(temp), "\n")
  
  lpm1 <- lm(reformulate(ctrl_lin, "gdp_y"), data = temp)
  lpm2 <- lm(reformulate(ctrl_lin, "cpi_y"), data = temp)
  lpm3 <- lm(reformulate(ctrl_lin, "ffr_y"), data = temp)
  se1 <- sqrt(diag(NeweyWest(lpm1, prewhite = FALSE)))
  se2 <- sqrt(diag(NeweyWest(lpm2, prewhite = FALSE)))
  se3 <- sqrt(diag(NeweyWest(lpm3, prewhite = FALSE)))
  
  adlpm1 <- lm(as.formula(paste("gdp_y ~", sd_formula)), data = temp)
  adlpm2 <- lm(as.formula(paste("cpi_y ~", sd_formula)), data = temp)
  adlpm3 <- lm(as.formula(paste("ffr_y ~", sd_formula)), data = temp)
  sd1 <- sqrt(diag(NeweyWest(adlpm1, prewhite = FALSE)))
  sd2 <- sqrt(diag(NeweyWest(adlpm2, prewhite = FALSE)))
  sd3 <- sqrt(diag(NeweyWest(adlpm3, prewhite = FALSE)))
  
  fill <- function(est, lpm, adlpm, se, sd) {
    c(-coef(lpm)["shock"],
      -coef(adlpm)["I(p1 * shock)"],
      -coef(adlpm)["I(n1 * shock)"],
      -coef(lpm)["shock"]   - conf*se["shock"],
      -coef(lpm)["shock"],
      -coef(lpm)["shock"]   + conf*se["shock"],
      -coef(adlpm)["I(p1 * shock)"] - conf*sd["I(p1 * shock)"],
      -coef(adlpm)["I(p1 * shock)"],
      -coef(adlpm)["I(p1 * shock)"] + conf*sd["I(p1 * shock)"],
      -coef(adlpm)["I(n1 * shock)"] - conf*sd["I(n1 * shock)"],
      -coef(adlpm)["I(n1 * shock)"],
      -coef(adlpm)["I(n1 * shock)"] + conf*sd["I(n1 * shock)"])
  }
  est_1[h+1, ] <- fill(est_1, lpm1, adlpm1, se1, sd1)
  est_2[h+1, ] <- fill(est_2, lpm2, adlpm2, se2, sd2)
  est_3[h+1, ] <- fill(est_3, lpm3, adlpm3, se3, sd3)
}

write.csv(est_1, "est_1_R.csv", row.names = FALSE)
write.csv(est_2, "est_2_R.csv", row.names = FALSE)
write.csv(est_3, "est_3_R.csv", row.names = FALSE)


cat("\nDone! Results saved to est_1_R.csv, est_2_R.csv, est_3_R.csv\n")

