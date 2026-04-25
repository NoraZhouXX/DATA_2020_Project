# Reproduce Figure 1: Dynamic Housing Market Connectedness Index

# Reads Ct_monthly.csv (output of Compute_Ct.R), converts to quarterly,
# applies HP filter (lambda=1600), and plots level + trend with NBER recession

rm(list = ls())

# 1. Load and aggregate monthly Ct to quarterly
ct_raw <- read.csv("Ct_monthly.csv", header = TRUE)[, 2]  # 432 monthly values

T_q  <- length(ct_raw) / 3
ct_q <- sapply(1:T_q, function(j) mean(ct_raw[(3*(j-1)+1):(3*j)]))

# Decimal year: 1981Q1 = 1981.0, step = 0.25
yr <- seq(1981.0, by = 0.25, length.out = T_q)

cat("Quarterly Ct:", T_q, "observations,", yr[1], "to", tail(yr,1), "\n")
cat("Range:", round(range(ct_q), 3), "\n\n")

# 2. HP filter (lambda = 1600, no external packages needed)
# Analytical solution: trend = (I + lambda * K'K)^{-1} y
hp_filter <- function(y, lambda = 1600) {
  n <- length(y)
  I <- diag(n)
  D <- diff(I, differences = 2)
  trend <- solve(I + lambda * t(D) %*% D, y)
  list(trend = as.numeric(trend), cycle = as.numeric(y - trend))
}

hp       <- hp_filter(ct_q, lambda = 1600)
ct_trend <- hp$trend
ct_cycle <- hp$cycle

# 3. Dummy variables (saved for later regression use)
# pshock = 1: above-trend connectedness (high); nshock = 1: below-trend (low)
# p1, n1: lagged one quarter
pshock <- as.integer(ct_cycle >= 0)
nshock <- as.integer(ct_cycle <  0)
p1     <- c(NA, pshock[-length(pshock)])
n1     <- c(NA, nshock[-length(nshock)])

write.csv(data.frame(year=yr, pshock=pshock, nshock=nshock, p1=p1, n1=n1),
          file = "Ct_dummies.csv", row.names = FALSE)
cat("Saved: Ct_dummies.csv\n")

# 4. Plot Figure 1 (base R, no packages required)

# Restrict to paper's sample: 1981Q1 - 2007Q4
x_min <- 1981.0
x_max <- 2007.75
idx   <- yr >= x_min & yr <= x_max
yr_p  <- yr[idx];  ct_p <- ct_q[idx];  tr_p <- ct_trend[idx]

# NBER recession shading bands (start, end in decimal years)
rec <- data.frame(
  start = c(1981.50, 1990.75, 2001.25, 2007.75),
  end   = c(1982.75, 1991.00, 2001.75, 2007.75)  # clipped to x_max
)

# Save to PDF and PNG
for (device in c("pdf", "png")) {
  if (device == "pdf") pdf("Fig_1.pdf", width = 7, height = 4)
  if (device == "png") png("Fig_1.png", width = 7, height = 4,
                           units = "in", res = 300)

  par(mar = c(4, 4, 1, 1))

  # Empty plot
  plot(yr_p, ct_p, type = "n",
       xlim = c(x_min, x_max), ylim = c(88, 96),
       xlab = "Year", ylab = "Connectedness Index (%)",
       xaxt = "n", bty = "l")
  axis(1, at = seq(1982, 2006, by = 4))
  grid(nx = NA, ny = NULL, lty = 1, col = "grey90")

  # Recession shading
  for (i in seq_len(nrow(rec))) {
    rect(rec$start[i], -1e9, rec$end[i], 1e9,
         col = rgb(0.8, 0.8, 0.8, 0.6), border = NA)
  }

  # Connectedness level (black solid)
  lines(yr_p, ct_p, col = "black", lwd = 1.5)

  # HP trend (blue dashed)
  lines(yr_p, tr_p, col = "blue", lwd = 2, lty = 2)

  # Legend
  legend("topleft", legend = c("Connectedness (Level)", "Connectedness (Trend)"),
         col = c("black", "blue"), lwd = c(1.5, 2), lty = c(1, 2),
         bty = "n", cex = 0.9)

  dev.off()
}

cat("Saved: Fig_1.pdf and Fig_1.png\n")

