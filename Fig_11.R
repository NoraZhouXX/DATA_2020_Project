# =============================================================================
# Fig_11.R
# IRFs including ZLB (extended sample 1981Q1-2016Q4), Shadow Rate shock
#
# Layout: 2 rows × 3 cols
#   Row (a) "Headline Variables":           GDP / PCE / Shadow Rate
#   Row (b) "Other Expenditure Variables":  Consumption / Investment / Residential Investment
# Each panel: Linear (black solid) + High (blue dashed + CI) + Low (red dotted + CI)
#
# Lab methods: lm() (Lab 4), sandwich::NeweyWest(), ggplot2 (Lab 2/4/9)
# INPUTS:  All_data.xls, Ct_dummies.csv
# OUTPUT:  Fig_11.png
# =============================================================================
rm(list = ls())
library(readxl); library(sandwich); library(ggplot2); library(patchwork)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
macro    <- read_excel("All_data.xls", sheet = "Macro_data")
n        <- nrow(macro)
yr       <- seq(1975.0, by = 0.25, length.out = n)
log_GDP  <- log(as.numeric(macro$RGDP))        * 100
log_CPI  <- log(as.numeric(macro$PCE))         * 100
SHADOW   <- as.numeric(macro$SHADOW_RATE)      # Wu-Xia shadow FFR
log_CON  <- log(as.numeric(macro$CONSUMPTION)) * 100
log_INV  <- log(as.numeric(macro$INVESTMENT))  * 100
log_RIN  <- log(as.numeric(macro$R_INVESTMET)) * 100  # residential investment

dum     <- read.csv("Ct_dummies.csv")
idx_dum <- match(round(dum$year, 4), round(yr, 4))
p1 <- rep(NA, n); p1[idx_dum] <- dum$p1
n1 <- rep(NA, n); n1[idx_dum] <- dum$n1

# Extended sample: 1981Q1 – 2016Q4
smpl   <- yr >= 1981.0 & yr <= 2016.75
T_smpl <- sum(smpl)
bw_fix <- floor(0.75 * T_smpl^(1/3))
cat("Extended sample:", T_smpl, "quarters | NW bandwidth:", bw_fix, "\n")

# ---------------------------------------------------------------------------
# 2. LP helpers
# ---------------------------------------------------------------------------
lag_vec <- function(x, l) c(rep(NA, l), x[seq_len(n - l)])

run_lp_h <- function(Y_lead, bw) {
  df <- data.frame(Y = Y_lead, SHOCK = SHADOW, p1 = p1, n1 = n1,
                   trend = seq_len(n), trend2 = seq_len(n)^2)
  df$p1_SHOCK <- p1 * SHADOW;  df$n1_SHOCK <- n1 * SHADOW
  for (l in 0:4) {
    df[[paste0("gdp_l", l)]]    <- lag_vec(log_GDP, l)
    df[[paste0("cpi_l", l)]]    <- lag_vec(log_CPI, l)
    df[[paste0("p1_gdp_l", l)]] <- p1 * df[[paste0("gdp_l", l)]]
    df[[paste0("n1_gdp_l", l)]] <- n1 * df[[paste0("gdp_l", l)]]
    df[[paste0("p1_cpi_l", l)]] <- p1 * df[[paste0("cpi_l", l)]]
    df[[paste0("n1_cpi_l", l)]] <- n1 * df[[paste0("cpi_l", l)]]
  }
  # Shadow rate lags as controls (shadow rate replaces FFR in extended sample)
  for (l in 1:4) {
    df[[paste0("ffr_l", l)]]    <- lag_vec(SHADOW, l)
    df[[paste0("p1_ffr_l", l)]] <- p1 * df[[paste0("ffr_l", l)]]
    df[[paste0("n1_ffr_l", l)]] <- n1 * df[[paste0("ffr_l", l)]]
  }
  df <- df[smpl, ];  df <- df[complete.cases(df), ]

  lin_rhs <- c("SHOCK", paste0("gdp_l", 0:4), paste0("cpi_l", 0:4),
               paste0("ffr_l", 1:4), "trend", "trend2")
  m_lin <- lm(as.formula(paste("Y ~", paste(lin_rhs, collapse = " + "))), data = df)
  V_lin <- NeweyWest(m_lin, lag = bw, prewhite = FALSE)

  sd_rhs <- c("0", "p1", "n1", "p1_SHOCK", "n1_SHOCK",
              paste0("p1_gdp_l", 0:4), paste0("n1_gdp_l", 0:4),
              paste0("p1_cpi_l", 0:4), paste0("n1_cpi_l", 0:4),
              paste0("p1_ffr_l", 1:4), paste0("n1_ffr_l", 1:4),
              "trend", "trend2")
  m_sd <- lm(as.formula(paste("Y ~", paste(sd_rhs, collapse = " + "))), data = df)
  V_sd <- NeweyWest(m_sd, lag = bw, prewhite = FALSE)

  list(b_lin = -coef(m_lin)["SHOCK"],   s_lin = sqrt(V_lin["SHOCK","SHOCK"]),
       b_H   = -coef(m_sd)["p1_SHOCK"], s_H   = sqrt(V_sd["p1_SHOCK","p1_SHOCK"]),
       b_L   = -coef(m_sd)["n1_SHOCK"], s_L   = sqrt(V_sd["n1_SHOCK","n1_SHOCK"]))
}

# ---------------------------------------------------------------------------
# 3. Run LP for all variables
# ---------------------------------------------------------------------------
H_max   <- 20;  h_seq <- 0:H_max
vars_ls <- list(GDP = log_GDP, PCE = log_CPI, `Shadow Rate` = SHADOW,
                Consumption = log_CON, Investment = log_INV,
                `Residential Investment` = log_RIN)

cat("Running LP (Shadow Rate shock, extended sample) h = 0 to", H_max, "...\n")
results <- lapply(seq_along(h_seq), function(hi) {
  h <- h_seq[hi]
  lapply(names(vars_ls), function(vname) {
    y_v <- vars_ls[[vname]]
    run_lp_h(c(y_v[(h+1):n], rep(NA, h)), bw_fix)
  }) |> setNames(names(vars_ls))
})

# ---------------------------------------------------------------------------
# 4. Compile IRF data frame
# ---------------------------------------------------------------------------
irf_df <- do.call(rbind, lapply(seq_along(h_seq), function(hi) {
  h <- h_seq[hi];  res <- results[[hi]]
  do.call(rbind, lapply(names(vars_ls), function(vname) {
    r <- res[[vname]]
    data.frame(
      h = h, var = vname,
      model = c("Linear","High","Low"),
      est   = c(r$b_lin, r$b_H,          r$b_L),
      lower = c(r$b_lin - r$s_lin, r$b_H - r$s_H, r$b_L - r$s_L),
      upper = c(r$b_lin + r$s_lin, r$b_H + r$s_H, r$b_L + r$s_L),
      stringsAsFactors = FALSE
    )
  }))
}))

var_order <- c("GDP","PCE","Shadow Rate","Consumption","Investment","Residential Investment")
irf_df$var   <- factor(irf_df$var,   levels = var_order)
irf_df$model <- factor(irf_df$model, levels = c("Linear","High","Low"))

# ---------------------------------------------------------------------------
# 5. Panel builder (same style as Fig_10)
#    Linear = black solid; High = blue dashed + blue CI band; Low = red dotted + dotted CI lines
# ---------------------------------------------------------------------------
make_panel <- function(vname, title_label) {
  sub  <- irf_df[irf_df$var == vname, ]
  lin  <- sub[sub$model == "Linear", ]
  high <- sub[sub$model == "High",   ]
  low  <- sub[sub$model == "Low",    ]

  ggplot() +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    # High CI band (blue shading)
    geom_ribbon(data = high, aes(x = h, ymin = lower, ymax = upper),
                fill = "#2166ac", alpha = 0.2) +
    # Low CI boundary lines (red dotted)
    geom_line(data = low, aes(x = h, y = lower),
              colour = "#d6604d", linetype = "dotted", linewidth = 0.7) +
    geom_line(data = low, aes(x = h, y = upper),
              colour = "#d6604d", linetype = "dotted", linewidth = 0.7) +
    # Point estimates
    geom_line(data  = lin,  aes(x = h, y = est),
              colour = "black",   linetype = "solid",  linewidth = 1.0) +
    geom_point(data = lin,  aes(x = h, y = est),
               colour = "black",  shape = 16, size = 1.5) +
    geom_line(data  = high, aes(x = h, y = est),
              colour = "#2166ac", linetype = "dashed", linewidth = 0.8) +
    geom_point(data = high, aes(x = h, y = est),
               colour = "#2166ac", shape = 17, size = 1.5) +
    geom_line(data  = low,  aes(x = h, y = est),
              colour = "#d6604d", linetype = "dotted", linewidth = 0.8) +
    geom_point(data = low,  aes(x = h, y = est),
               colour = "#d6604d", shape = 15, size = 1.5) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20)) +
    labs(x = "Horizon (Quarter)", y = NULL, title = title_label) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
}

# ---------------------------------------------------------------------------
# 6. Assemble 2×3 grid with centered row titles
#    Row (a): GDP / PCE / Shadow Rates
#    Row (b): Consumption / Investment / Residential Investment
# ---------------------------------------------------------------------------

# Row (a) — Headline Variables
p_gdp  <- make_panel("GDP",          "GDP") +
  labs(tag = "(a)") +
  theme(plot.tag = element_text(face = "bold", size = 13))
p_pce  <- make_panel("PCE",          "PCE")
p_shd  <- make_panel("Shadow Rate",  "Shadow Rates")

row_a <- (p_gdp | p_pce | p_shd) +
  plot_annotation(
    title = "Headline Variables",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 11,
                                face = "bold", margin = margin(b = 4))
    )
  )

# Row (b) — Other Expenditure Variables
p_con  <- make_panel("Consumption",           "Consumption") +
  labs(tag = "(b)") +
  theme(plot.tag = element_text(face = "bold", size = 13))
p_inv  <- make_panel("Investment",            "Investment")
p_rin  <- make_panel("Residential Investment","Residential Investment")

row_b <- (p_con | p_inv | p_rin) +
  plot_annotation(
    title = "Other Expenditure Variables",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 11,
                                face = "bold", margin = margin(b = 4))
    )
  )

fig11 <- row_a / row_b
ggsave("Fig_11.png", fig11, width = 12, height = 8, dpi = 200)
cat("Saved: Fig_11.png\n")
