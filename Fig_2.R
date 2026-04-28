# =============================================================================
# Fig_2.R  (lab version)
# Reproduce Figure 2: Baseline IRFs — Local Projections (Lee 2025)
#
# Lab methods used:
#   - lm()               (Lab 4): OLS for each LP horizon
#   - sandwich::NeweyWest() (Lab 4/5): HAC standard errors
#   - ggplot2 + patchwork   (Lab 2/4/9): 3x3 IRF panel figure
#   - Dummy interactions    (Lab 9 DiD style): state-dependent LP
#
# INPUTS:  All_data.xls, Ct_dummies.csv
# OUTPUTS: Fig_2.png
# =============================================================================

rm(list = ls())

# install.packages(c("readxl","sandwich","ggplot2","patchwork"))  # run once
library(readxl)
library(sandwich)    # NeweyWest() — HAC covariance matrix
library(ggplot2)     # Lab 2/4/9
library(patchwork)   # combine ggplot panels

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
macro   <- read_excel("All_data.xls", sheet = "Macro_data")
log_GDP <- log(as.numeric(macro$RGDP)) * 100
log_CPI <- log(as.numeric(macro$PCE))  * 100
FFR     <- as.numeric(macro$FFR)
n       <- nrow(macro)
yr      <- seq(1975.0, by = 0.25, length.out = n)

# Align connectedness dummies with macro data
dum     <- read.csv("Ct_dummies.csv", header = TRUE)
idx_dum <- match(round(dum$year, 4), round(yr, 4))
p1 <- rep(NA, n);  p1[idx_dum] <- dum$p1
n1 <- rep(NA, n);  n1[idx_dum] <- dum$n1

smpl   <- yr >= 1981.0 & yr <= 2007.75
T_smpl <- sum(smpl)
bw_fix <- floor(0.75 * T_smpl^(1/3))
cat("Sample:", T_smpl, "quarters | NW bandwidth:", bw_fix, "\n\n")

# ---------------------------------------------------------------------------
# 2. Build regression data frame (lags + interaction columns)
#    Dummy interactions follow Lab 9 (DiD) style:
#    High regime = p1 * FFR,  Low regime = n1 * FFR
# ---------------------------------------------------------------------------
lag_vec <- function(x, l) c(rep(NA, l), x[seq_len(n - l)])

make_reg_df <- function(Y_lead) {
  df <- data.frame(
    Y      = Y_lead,
    FFR    = FFR,
    p1     = p1,
    n1     = n1,
    trend  = seq_len(n),
    trend2 = seq_len(n)^2
  )
  # Lag columns for GDP, CPI, FFR
  for (l in 0:4) {
    df[[paste0("gdp_l", l)]] <- lag_vec(log_GDP, l)
    df[[paste0("cpi_l", l)]] <- lag_vec(log_CPI, l)
  }
  for (l in 1:4) {
    df[[paste0("ffr_l", l)]] <- lag_vec(FFR, l)
  }
  # State-dependent interaction columns (Lab 9 style)
  df$p1_FFR <- df$p1 * df$FFR
  df$n1_FFR <- df$n1 * df$FFR
  for (l in 0:4) {
    df[[paste0("p1_gdp_l", l)]] <- df$p1 * df[[paste0("gdp_l", l)]]
    df[[paste0("n1_gdp_l", l)]] <- df$n1 * df[[paste0("gdp_l", l)]]
    df[[paste0("p1_cpi_l", l)]] <- df$p1 * df[[paste0("cpi_l", l)]]
    df[[paste0("n1_cpi_l", l)]] <- df$n1 * df[[paste0("cpi_l", l)]]
  }
  for (l in 1:4) {
    df[[paste0("p1_ffr_l", l)]] <- df$p1 * df[[paste0("ffr_l", l)]]
    df[[paste0("n1_ffr_l", l)]] <- df$n1 * df[[paste0("ffr_l", l)]]
  }
  df[smpl, ]
}

# ---------------------------------------------------------------------------
# 3. LP at horizon h using lm() + NeweyWest()  (Lab 4 / sandwich)
# ---------------------------------------------------------------------------
run_lp_h <- function(Y_lead, bw) {

  df <- make_reg_df(Y_lead)
  df <- df[complete.cases(df), ]

  # --- Linear model (Lab 4: lm) ---
  lin_rhs <- c("FFR",
               paste0("gdp_l", 0:4), paste0("cpi_l", 0:4),
               paste0("ffr_l", 1:4), "trend", "trend2")
  m_lin <- lm(as.formula(paste("Y ~", paste(lin_rhs, collapse = " + "))),
              data = df)
  V_lin <- NeweyWest(m_lin, lag = bw, prewhite = FALSE)
  b_lin <- -coef(m_lin)["FFR"]
  s_lin <- sqrt(V_lin["FFR", "FFR"])

  # --- State-dependent model (Lab 9 DiD style: no intercept, dummies interact
  #     all controls — High connectedness = p1, Low = n1) ---
  sd_rhs <- c("0", "p1_FFR", "n1_FFR",
               paste0("p1_gdp_l", 0:4), paste0("n1_gdp_l", 0:4),
               paste0("p1_cpi_l", 0:4), paste0("n1_cpi_l", 0:4),
               paste0("p1_ffr_l", 1:4), paste0("n1_ffr_l", 1:4),
               "trend", "trend2")
  m_sd <- lm(as.formula(paste("Y ~", paste(sd_rhs, collapse = " + "))),
             data = df)
  V_sd <- NeweyWest(m_sd, lag = bw, prewhite = FALSE)
  b_H  <- -coef(m_sd)["p1_FFR"]
  b_L  <- -coef(m_sd)["n1_FFR"]
  s_H  <- sqrt(V_sd["p1_FFR", "p1_FFR"])
  s_L  <- sqrt(V_sd["n1_FFR", "n1_FFR"])

  list(b_lin = b_lin, s_lin = s_lin,
       b_H   = b_H,   s_H   = s_H,
       b_L   = b_L,   s_L   = s_L)
}

# ---------------------------------------------------------------------------
# 4. Run LP for all horizons h = 0 to 20 and 3 outcome variables
# ---------------------------------------------------------------------------
H_max   <- 20
h_seq   <- 0:H_max
vars_ls <- list(GDP = log_GDP, CPI = log_CPI, FFR = FFR)
results <- list()

cat("Running LP h = 0 to", H_max, "...\n")
for (hi in seq_along(h_seq)) {
  h      <- h_seq[hi]
  res_h  <- list()
  for (vname in names(vars_ls)) {
    y_v              <- vars_ls[[vname]]
    Y_lead           <- c(y_v[(h+1):n], rep(NA, h))
    res_h[[vname]]   <- run_lp_h(Y_lead, bw_fix)
  }
  results[[hi]] <- res_h
  if (h %% 5 == 0)
    cat(sprintf("h=%2d  GDP: Lin=%.3f  High=%.3f  Low=%.3f\n",
                h, results[[hi]]$GDP$b_lin,
                   results[[hi]]$GDP$b_H,
                   results[[hi]]$GDP$b_L))
}

# ---------------------------------------------------------------------------
# 5. Compile tidy IRF data frame for ggplot2
# ---------------------------------------------------------------------------
irf_df <- do.call(rbind, lapply(seq_along(h_seq), function(hi) {
  h   <- h_seq[hi]
  res <- results[[hi]]
  do.call(rbind, lapply(names(vars_ls), function(vname) {
    r <- res[[vname]]
    data.frame(h     = h,
               var   = vname,
               model = c("Linear", "High", "Low"),
               est   = c(r$b_lin, r$b_H, r$b_L),
               se    = c(r$s_lin, r$s_H, r$s_L),
               stringsAsFactors = FALSE)
  }))
}))
irf_df$lower <- irf_df$est - irf_df$se
irf_df$upper <- irf_df$est + irf_df$se
irf_df$var   <- factor(irf_df$var,   levels = c("GDP","CPI","FFR"),
                       labels = c("GDP","PCE Deflator","FFR"))
irf_df$model <- factor(irf_df$model, levels = c("Linear","High","Low"))

# ---------------------------------------------------------------------------
# 6. ggplot2 plot functions  (Lab 2/4/9)
# ---------------------------------------------------------------------------

# Colour / fill / shape scheme
col_map  <- c(Linear = "black", High = "#2166ac", Low = "#d6604d")
fill_map <- c(Linear = "grey80", High = "#d1e5f0", Low = "#fddbc7")
shp_map  <- c(Linear = 16, High = 17, Low = 15)

# Panel A: Three models overlaid (point estimates, no CI)
plot_three <- function(var_label) {
  sub <- irf_df[irf_df$var == var_label, ]
  ggplot(sub, aes(x = h, y = est, colour = model, shape = model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = col_map) +
    scale_shape_manual(values  = shp_map) +
    scale_x_continuous(breaks  = c(0,5,10,15,20)) +
    labs(x = "Horizon (quarters)", y = var_label,
         title = "Three Models",
         colour = NULL, shape = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position   = "top",
          panel.grid.minor  = element_blank())
}

# Panel B: Linear model with ±1 SE ribbon
plot_linear <- function(var_label) {
  sub <- irf_df[irf_df$var == var_label & irf_df$model == "Linear", ]
  ggplot(sub, aes(x = h)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.6) +
    geom_line(aes(y = est), colour = "black", linewidth = 0.8) +
    geom_point(aes(y = est), colour = "black", size = 1.8) +
    scale_x_continuous(breaks = c(0,5,10,15,20)) +
    labs(x = "Horizon (quarters)", y = var_label,
         title = "Linear Model") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
}

# Panel C: State-dependent (High = blue, Low = red) with ±1 SE ribbons
plot_state <- function(var_label) {
  sub_H <- irf_df[irf_df$var == var_label & irf_df$model == "High", ]
  sub_L <- irf_df[irf_df$var == var_label & irf_df$model == "Low",  ]
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_ribbon(data = sub_H, aes(x = h, ymin = lower, ymax = upper),
                fill = "#d1e5f0", alpha = 0.6) +
    geom_ribbon(data = sub_L, aes(x = h, ymin = lower, ymax = upper),
                fill = "#fddbc7", alpha = 0.6) +
    geom_line(data  = sub_H, aes(x = h, y = est),
              colour = "#2166ac", linewidth = 0.8) +
    geom_point(data = sub_H, aes(x = h, y = est),
               colour = "#2166ac", shape = 17, size = 1.8) +
    geom_line(data  = sub_L, aes(x = h, y = est),
              colour = "#d6604d", linewidth = 0.8, linetype = "dashed") +
    geom_point(data = sub_L, aes(x = h, y = est),
               colour = "#d6604d", shape = 15, size = 1.8) +
    scale_x_continuous(breaks = c(0,5,10,15,20)) +
    labs(x = "Horizon (quarters)", y = var_label,
         title = "State-Dependent Models") +
    annotate("text", x = 18, y = Inf, vjust = 1.5, size = 2.8,
             label = "▲ High   ■ Low", colour = "grey30") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
}

# ---------------------------------------------------------------------------
# 7. Assemble 3×3 grid with patchwork and save
# ---------------------------------------------------------------------------
var_labs <- c("GDP", "PCE Deflator", "FFR")

panels <- lapply(var_labs, function(vl)
  plot_three(vl) | plot_linear(vl) | plot_state(vl))

fig2 <- panels[[1]] / panels[[2]] / panels[[3]]

ggsave("Fig_2.png", fig2, width = 12, height = 10, dpi = 200)
cat("Saved: Fig_2.png\n")
