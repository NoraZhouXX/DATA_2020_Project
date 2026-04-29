# =============================================================================
# Fig_1.R  (lab version)
# Reproduce Figure 1: Dynamic Housing Market Connectedness Index
#
# Lab methods used:
#   - ggplot2  (Lab 2 / 4 / 9): time-series plot with recession shading
#   - sapply / data.frame (Lab 3): aggregate monthly → quarterly
#
# INPUTS:  Ct_monthly.csv  (output of Compute_Ct.R)
# OUTPUTS: Ct_dummies.csv, Fig_1.png
# =============================================================================

rm(list = ls())

# install.packages(c("ggplot2"))   # run once if needed
library(ggplot2)

# ---------------------------------------------------------------------------
# 1. Load monthly Ct and aggregate to quarterly averages
# ---------------------------------------------------------------------------
ct_raw <- read.csv("../data/Ct_monthly.csv", header = TRUE)[, 2]   # 432 monthly values
T_q    <- length(ct_raw) / 3
ct_q   <- sapply(1:T_q, function(j) mean(ct_raw[(3*(j-1)+1):(3*j)]))
yr     <- seq(1981.0, by = 0.25, length.out = T_q)

cat("Quarterly Ct:", T_q, "obs |", yr[1], "to", tail(yr, 1), "\n")
cat("Range:", round(range(ct_q), 3), "\n\n")

# ---------------------------------------------------------------------------
# 2. HP filter (lambda = 1600) — analytical solution (no external package)
#    trend = (I + lambda * D'D)^{-1} y
# ---------------------------------------------------------------------------
hp_filter <- function(y, lambda = 1600) {
  n <- length(y)
  D <- diff(diag(n), differences = 2)
  trend <- solve(diag(n) + lambda * t(D) %*% D, y)
  list(trend = as.numeric(trend), cycle = as.numeric(y - trend))
}

hp       <- hp_filter(ct_q, lambda = 1600)
ct_trend <- hp$trend
ct_cycle <- hp$cycle

# ---------------------------------------------------------------------------
# 3. Regime dummies → save Ct_dummies.csv (used by Table_1.R / Fig_2.R)
#    pshock / nshock: above / below HP trend
#    p1 / n1: lagged one quarter (instrument enters at t-1)
# ---------------------------------------------------------------------------
pshock <- as.integer(ct_cycle >= 0)
nshock <- as.integer(ct_cycle <  0)
p1     <- c(NA, pshock[-length(pshock)])
n1     <- c(NA, nshock[-length(nshock)])

write.csv(data.frame(year   = yr,
                     pshock = pshock,
                     nshock = nshock,
                     p1     = p1,
                     n1     = n1),
          file = "../data/Ct_dummies.csv", row.names = FALSE)
cat("Saved: Ct_dummies.csv\n")

# ---------------------------------------------------------------------------
# 4. ggplot2 figure  (Lab 2 / 4 / 9 style)
# ---------------------------------------------------------------------------

# Restrict to paper sample: 1981Q1–2007Q4
idx    <- yr >= 1981.0 & yr <= 2007.75
df_plt <- data.frame(
  year   = rep(yr[idx], 2),
  Ct     = c(ct_q[idx], ct_trend[idx]),
  Series = rep(c("Connectedness (Level)", "Connectedness (Trend)"),
               each = sum(idx))
)
df_plt$Series <- factor(df_plt$Series,
                        levels = c("Connectedness (Level)",
                                   "Connectedness (Trend)"))

# NBER recession shading bands
rec <- data.frame(
  start = c(1981.50, 1990.75, 2001.25),
  end   = c(1982.75, 1991.00, 2001.75)
)

p <- ggplot(df_plt, aes(x = year, y = Ct,
                         colour   = Series,
                         linetype = Series)) +
  # Recession shading (drawn first so lines appear on top)
  geom_rect(data = rec, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey80", alpha = 0.6) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values   = c("black", "blue")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_continuous(breaks = seq(1982, 2006, by = 4),
                     limits = c(1981, 2008)) +
  labs(x        = "Year",
       y        = "Connectedness Index (%)",
       colour   = NULL,
       linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position   = "top",
        legend.key.width  = unit(1.5, "cm"),
        panel.grid.minor  = element_blank())

ggsave("../output/Fig_1.png", p, width = 7, height = 4, dpi = 300)
cat("Saved: Fig_1.png\n")
