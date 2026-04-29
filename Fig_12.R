# =============================================================================
# Fig_12.R
# Connectedness index: Benchmark vs. Macro-controls
#
# Two lines: Connectedness (Benchmark) [black solid]  from Ct_monthly.csv
#            Connectedness (Macro)     [blue dashed]  from Ct_macro_monthly.csv
# Gray shaded regions: NBER recession periods (1981–2007)
#
# Lab methods: ggplot2 (Lab 2/4/9)
# INPUTS:  Ct_monthly.csv, Ct_macro_monthly.csv
# OUTPUT:  Fig_12.png
# =============================================================================
rm(list = ls())
library(ggplot2)

# ---------------------------------------------------------------------------
# 1. Load Ct series
# ---------------------------------------------------------------------------
ct_bench <- read.csv("Ct_monthly.csv", header = TRUE)[, 2]
ct_macro <- read.csv("Ct_macro_monthly.csv", header = TRUE)[, 2]

stopifnot(length(ct_bench) == length(ct_macro))
T_total <- length(ct_bench)

# Dates: rolling window of 120 months, p=1; first Ct is dated at end of
# first window → month 121 = January 1981 (data starts Jan 1971)
# Adjust start month to match your EI_DATA.csv start date if needed.
start_year  <- 1981
start_month <- 1
dates <- seq(as.Date(paste(start_year, start_month, "01", sep = "-")),
             by = "month", length.out = T_total)

df <- data.frame(
  date      = dates,
  Benchmark = ct_bench,
  Macro     = ct_macro
)

cat("Correlation (Benchmark, Macro):", round(cor(ct_bench, ct_macro), 3), "\n")

# ---------------------------------------------------------------------------
# 2. NBER recession bands (monthly; within 1981-2007 sample)
# ---------------------------------------------------------------------------
recessions <- data.frame(
  start = as.Date(c("1981-07-01", "1990-07-01", "2001-03-01")),
  end   = as.Date(c("1982-11-30", "1991-03-31", "2001-11-30"))
)

# ---------------------------------------------------------------------------
# 3. Plot
# ---------------------------------------------------------------------------
fig12 <- ggplot(df, aes(x = date)) +
  # NBER recession shading (drawn first, behind lines)
  geom_rect(data = recessions,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "grey70", alpha = 0.5) +
  # Benchmark: black solid
  geom_line(aes(y = Benchmark, colour = "Connectedness (Benchmark)",
                linetype = "Connectedness (Benchmark)"), linewidth = 0.8) +
  # Macro: blue dashed
  geom_line(aes(y = Macro, colour = "Connectedness (Macro)",
                linetype = "Connectedness (Macro)"), linewidth = 0.8) +
  scale_colour_manual(
    values = c("Connectedness (Benchmark)" = "black",
               "Connectedness (Macro)"     = "#2166ac")
  ) +
  scale_linetype_manual(
    values = c("Connectedness (Benchmark)" = "solid",
               "Connectedness (Macro)"     = "dashed")
  ) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y",
               limits = c(as.Date("1981-01-01"), as.Date("2007-12-31"))) +
  labs(x = "Year", y = NULL, colour = NULL, linetype = NULL) +
  theme_bw(base_size = 10) +
  theme(
    legend.position   = c(0.25, 0.88),
    legend.background = element_rect(fill = "white", colour = "black",
                                     linewidth = 0.3),
    legend.key.width  = unit(1.8, "cm"),
    panel.grid.minor  = element_blank()
  )

ggsave("Fig_12.png", fig12, width = 10, height = 5, dpi = 200)
cat("Saved: Fig_12.png\n")
