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
macro <- read_excel("All_data.xls", sheet = "Macro_data") 

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
ci_dates <- seq(from       = as.Date("1981-01-01"),
                by         = "quarter",
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
hp_result <- hpfilter(ci_quarterly, freq = 1600)

ci_df$c_trend <- hp_result$trend
ci_df$c_cycle <- hp_result$cycle

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

ci_df$p1[1] <- ci_df$pshock[1]
ci_df$n1[1] <- ci_df$nshock[1]

# ============================================================
# STEP 6 - Merge macro and connectedness data
# ============================================================

df <- merge(macro, ci_df, by = "date_q", all.x = TRUE)
df <- df[order(df$date_q), ]
rownames(df) <- NULL

# ============================================================
# STEP 7 - Restrict to baseline sample 1981Q1 - 2007Q4
# ============================================================

start_date <- as.Date("1981-01-01")
end_date   <- as.Date("2007-10-01")

sample_df <- df[df$date_q >= start_date & df$date_q <= end_date, ]
sample_df <- sample_df[complete.cases(sample_df[, c("log_GDP", "log_cpi",
                                                    "shock", "p1", "n1")]), ]
rownames(sample_df) <- NULL

# Quick check
cat("Rows in sample_df:", nrow(sample_df), "\n")
cat("First date:", as.character(sample_df$date_q[1]), "\n")
cat("Last date:",  as.character(sample_df$date_q[nrow(sample_df)]), "\n")
cat("High connectedness (p1=1):", sum(sample_df$p1 == 1), "\n")
cat("Low connectedness  (n1=1):", sum(sample_df$n1 == 1), "\n")

# ============================================================
# STEP 8 - Create lagged variables
# ============================================================

make_lags <- function(x, max_lag) {
  n      <- length(x)
  result <- matrix(NA, nrow = n, ncol = max_lag)
  for (i in 1:max_lag) {
    result[(i+1):n, i] <- x[1:(n-i)]
  }
  return(result)
}

# Lags of GDP, CPI, FFR (up to 4 lags)
gdp_lags <- make_lags(sample_df$log_GDP, 4)
cpi_lags <- make_lags(sample_df$log_cpi, 4)
ffr_lags <- make_lags(sample_df$shock,   4)

colnames(gdp_lags) <- paste0("gdp_lag", 1:4)
colnames(cpi_lags) <- paste0("cpi_lag", 1:4)
colnames(ffr_lags) <- paste0("ffr_lag", 1:4)

sample_df <- cbind(sample_df, gdp_lags, cpi_lags, ffr_lags)

# Time trend
sample_df$trend  <- 1:nrow(sample_df) # captures linear growth over time
sample_df$trend2 <- sample_df$trend^2 # captures any acceleration or deceleration in that growth (quadratic)

# ============================================================
# STEP 9 - Local Projections Loop (h = 0 to 20)
# ============================================================

H    <- 21   # horizons 0 to 20
conf <- 1    # 1 standard error bands (68% confidence)

# Storage matrices
est_1 <- matrix(NA, nrow = H, ncol = 12)  # GDP
est_2 <- matrix(NA, nrow = H, ncol = 12)  # PCE
est_3 <- matrix(NA, nrow = H, ncol = 12)  # FFR

colnames(est_1) <- colnames(est_2) <- colnames(est_3) <-
  c("point_linear", "point_high", "point_low",
    "linear1",      "linear2",    "linear3",
    "H1",           "H2",         "H3",
    "L1",           "L2",         "L3")

N <- nrow(sample_df)

for (h in 0:20) {
  
  # ---- Create lead variables (dependent at horizon h) ----
  # we directly regress the future value of GDP on today's shock.
  if (h > 0) {
    gdp_lead <- c(sample_df$log_GDP[(h+1):N], rep(NA, h))
    cpi_lead <- c(sample_df$log_cpi[(h+1):N], rep(NA, h))
    ffr_lead <- c(sample_df$shock[(h+1):N],   rep(NA, h))
  } else {
    gdp_lead <- sample_df$log_GDP
    cpi_lead <- sample_df$log_cpi
    ffr_lead <- sample_df$shock
  }
  
  # ---- Keep only complete rows for this horizon ----
  # On average, what happens to GDP h quarters after a 1% increase in FFR?
  valid <- complete.cases(data.frame(
    gdp_lead, cpi_lead, ffr_lead,
    sample_df[, c("shock", "p1", "n1", "trend", "trend2",
                  "log_GDP", "gdp_lag1", "gdp_lag2", "gdp_lag3", "gdp_lag4",
                  "log_cpi", "cpi_lag1", "cpi_lag2", "cpi_lag3", "cpi_lag4",
                  "ffr_lag1", "ffr_lag2", "ffr_lag3", "ffr_lag4")]
  ))
  
  temp_df  <- sample_df[valid, ]
  gdp_y    <- gdp_lead[valid]
  cpi_y    <- cpi_lead[valid]
  ffr_y    <- ffr_lead[valid]
  
  # ---- LINEAR MODELS ----
  controls <- "shock +
               ffr_lag1 + ffr_lag2 + ffr_lag3 + ffr_lag4 +
               log_GDP + gdp_lag1 + gdp_lag2 + gdp_lag3 + gdp_lag4 +
               log_cpi + cpi_lag1 + cpi_lag2 + cpi_lag3 + cpi_lag4 +
               trend + trend2"
  
  lpm1 <- lm(as.formula(paste("gdp_y ~", controls)), data = temp_df)
  lpm2 <- lm(as.formula(paste("cpi_y ~", controls)), data = temp_df)
  lpm3 <- lm(as.formula(paste("ffr_y ~", controls)), data = temp_df)
  
  # HAC standard errors (Newey-West, no prewhitening)
  # Newey-West corrects by accounting for:
  #   > Heteroskedasticity — variance of errors changes over time
  #   > Autocorrelation — errors are correlated across time periods
  se1 <- sqrt(diag(NeweyWest(lpm1, prewhite = FALSE)))
  se2 <- sqrt(diag(NeweyWest(lpm2, prewhite = FALSE)))
  se3 <- sqrt(diag(NeweyWest(lpm3, prewhite = FALSE)))
  
  # ---- STATE-DEPENDENT MODELS (Formula 2) ----
  # βH = effect of FFR shock in high connectedness state
  # βL = effect of FFR shock in low connectedness state
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
  
  adlpm1 <- lm(as.formula(paste("gdp_y ~", sd_formula)), data = temp_df)
  adlpm2 <- lm(as.formula(paste("cpi_y ~", sd_formula)), data = temp_df)
  adlpm3 <- lm(as.formula(paste("ffr_y ~", sd_formula)), data = temp_df)
  
  # HAC standard errors (Newey-West, no prewhitening)
  sd_se1 <- sqrt(diag(NeweyWest(adlpm1, prewhite = FALSE)))
  sd_se2 <- sqrt(diag(NeweyWest(adlpm2, prewhite = FALSE)))
  sd_se3 <- sqrt(diag(NeweyWest(adlpm3, prewhite = FALSE)))
  
  # ---- STORE RESULTS (x -1 for expansionary shock) ----
  row <- h + 1
  
  # GDP (est_1)
  est_1[row,  1] <- -coef(lpm1)["shock"]
  est_1[row,  2] <- -coef(adlpm1)["I(p1 * shock)"]
  est_1[row,  3] <- -coef(adlpm1)["I(n1 * shock)"]
  est_1[row,  4] <- -coef(lpm1)["shock"] - conf * se1["shock"]
  est_1[row,  5] <- -coef(lpm1)["shock"]
  est_1[row,  6] <- -coef(lpm1)["shock"] + conf * se1["shock"]
  est_1[row,  7] <- -coef(adlpm1)["I(p1 * shock)"] - conf * sd_se1["I(p1 * shock)"]
  est_1[row,  8] <- -coef(adlpm1)["I(p1 * shock)"]
  est_1[row,  9] <- -coef(adlpm1)["I(p1 * shock)"] + conf * sd_se1["I(p1 * shock)"]
  est_1[row, 10] <- -coef(adlpm1)["I(n1 * shock)"] - conf * sd_se1["I(n1 * shock)"]
  est_1[row, 11] <- -coef(adlpm1)["I(n1 * shock)"]
  est_1[row, 12] <- -coef(adlpm1)["I(n1 * shock)"] + conf * sd_se1["I(n1 * shock)"]
  
  # PCE (est_2)
  est_2[row,  1] <- -coef(lpm2)["shock"]
  est_2[row,  2] <- -coef(adlpm2)["I(p1 * shock)"]
  est_2[row,  3] <- -coef(adlpm2)["I(n1 * shock)"]
  est_2[row,  4] <- -coef(lpm2)["shock"] - conf * se2["shock"]
  est_2[row,  5] <- -coef(lpm2)["shock"]
  est_2[row,  6] <- -coef(lpm2)["shock"] + conf * se2["shock"]
  est_2[row,  7] <- -coef(adlpm2)["I(p1 * shock)"] - conf * sd_se2["I(p1 * shock)"]
  est_2[row,  8] <- -coef(adlpm2)["I(p1 * shock)"]
  est_2[row,  9] <- -coef(adlpm2)["I(p1 * shock)"] + conf * sd_se2["I(p1 * shock)"]
  est_2[row, 10] <- -coef(adlpm2)["I(n1 * shock)"] - conf * sd_se2["I(n1 * shock)"]
  est_2[row, 11] <- -coef(adlpm2)["I(n1 * shock)"]
  est_2[row, 12] <- -coef(adlpm2)["I(n1 * shock)"] + conf * sd_se2["I(n1 * shock)"]
  
  # FFR (est_3)
  est_3[row,  1] <- -coef(lpm3)["shock"]
  est_3[row,  2] <- -coef(adlpm3)["I(p1 * shock)"]
  est_3[row,  3] <- -coef(adlpm3)["I(n1 * shock)"]
  est_3[row,  4] <- -coef(lpm3)["shock"] - conf * se3["shock"]
  est_3[row,  5] <- -coef(lpm3)["shock"]
  est_3[row,  6] <- -coef(lpm3)["shock"] + conf * se3["shock"]
  est_3[row,  7] <- -coef(adlpm3)["I(p1 * shock)"] - conf * sd_se3["I(p1 * shock)"]
  est_3[row,  8] <- -coef(adlpm3)["I(p1 * shock)"]
  est_3[row,  9] <- -coef(adlpm3)["I(p1 * shock)"] + conf * sd_se3["I(p1 * shock)"]
  est_3[row, 10] <- -coef(adlpm3)["I(n1 * shock)"] - conf * sd_se3["I(n1 * shock)"]
  est_3[row, 11] <- -coef(adlpm3)["I(n1 * shock)"]
  est_3[row, 12] <- -coef(adlpm3)["I(n1 * shock)"] + conf * sd_se3["I(n1 * shock)"]
  
  cat("Horizon", h, "done\n")
}

# ============================================================
# STEP 10 - Save results
# ============================================================

write.csv(est_1, "est_1_R.csv", row.names = FALSE)
write.csv(est_2, "est_2_R.csv", row.names = FALSE)
write.csv(est_3, "est_3_R.csv", row.names = FALSE)

cat("\nDone! Results saved to est_1_R.csv, est_2_R.csv, est_3_R.csv\n")
