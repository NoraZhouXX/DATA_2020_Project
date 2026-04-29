# =============================================================================
# Fig_4.R
# Other macro variables: Employment + Real Wage + House Price  (FFR shock)
# Lab methods: lm() (Lab 4), sandwich::NeweyWest(), ggplot2 (Lab 2/4/9)
# INPUTS:  All_data.xls, Ct_dummies.csv
# OUTPUT:  Fig_4.png
# =============================================================================
rm(list = ls())
library(readxl); library(sandwich); library(ggplot2); library(patchwork)

# --- Load data ---------------------------------------------------------------
macro   <- read_excel("All_data.xls", sheet = "Macro_data")
n       <- nrow(macro)
yr      <- seq(1975.0, by = 0.25, length.out = n)
log_GDP <- log(as.numeric(macro$RGDP)) * 100
log_CPI <- log(as.numeric(macro$PCE))  * 100
FFR     <- as.numeric(macro$FFR)
log_EMP <- log(as.numeric(macro$EMPLOYMENT))  * 100
log_WAG <- log(as.numeric(macro$REAL_WAGE))   * 100
log_HP  <- log(as.numeric(macro$HOUSE_PRICE)) * 100

dum     <- read.csv("Ct_dummies.csv")
idx_dum <- match(round(dum$year, 4), round(yr, 4))
p1 <- rep(NA, n); p1[idx_dum] <- dum$p1
n1 <- rep(NA, n); n1[idx_dum] <- dum$n1

smpl   <- yr >= 1981.0 & yr <= 2007.75
T_smpl <- sum(smpl)
bw_fix <- floor(0.75 * T_smpl^(1/3))
cat("Sample:", T_smpl, "quarters | NW bandwidth:", bw_fix, "\n")

# --- Helper functions --------------------------------------------------------
lag_vec <- function(x, l) c(rep(NA, l), x[seq_len(n - l)])

make_reg_df <- function(Y_lead) {
  df <- data.frame(Y = Y_lead, SHOCK = FFR, p1 = p1, n1 = n1,
                   trend = seq_len(n), trend2 = seq_len(n)^2)
  df$p1_SHOCK <- p1 * FFR;  df$n1_SHOCK <- n1 * FFR
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

run_lp_h <- function(Y_lead, bw) {
  df <- make_reg_df(Y_lead);  df <- df[complete.cases(df), ]

  lin_rhs <- c("SHOCK", paste0("gdp_l", 0:4), paste0("cpi_l", 0:4),
               paste0("ffr_l", 1:4), "trend", "trend2")
  m_lin <- lm(as.formula(paste("Y ~", paste(lin_rhs, collapse = " + "))), data = df)
  V_lin <- NeweyWest(m_lin, lag = bw, prewhite = FALSE)

  sd_rhs <- c("0", "p1_SHOCK", "n1_SHOCK",
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

# --- Run LP ------------------------------------------------------------------
H_max   <- 20;  h_seq <- 0:H_max
vars_ls <- list(Employment = log_EMP, `Real Wage` = log_WAG, `House Price` = log_HP)
results <- list()

cat("Running LP h = 0 to", H_max, "...\n")
for (hi in seq_along(h_seq)) {
  h <- h_seq[hi];  res_h <- list()
  for (vname in names(vars_ls)) {
    y_v <- vars_ls[[vname]]
    res_h[[vname]] <- run_lp_h(c(y_v[(h+1):n], rep(NA, h)), bw_fix)
  }
  results[[hi]] <- res_h
}

# --- Compile IRF data frame --------------------------------------------------
irf_df <- do.call(rbind, lapply(seq_along(h_seq), function(hi) {
  h <- h_seq[hi];  res <- results[[hi]]
  do.call(rbind, lapply(names(vars_ls), function(vname) {
    r <- res[[vname]]
    data.frame(h = h, var = vname, model = c("Linear","High","Low"),
               est = c(r$b_lin, r$b_H, r$b_L),
               se  = c(r$s_lin, r$s_H, r$s_L), stringsAsFactors = FALSE)
  }))
}))
irf_df$lower <- irf_df$est - irf_df$se
irf_df$upper <- irf_df$est + irf_df$se
irf_df$var   <- factor(irf_df$var,   levels = c("Employment","Real Wage","House Price"))
irf_df$model <- factor(irf_df$model, levels = c("Linear","High","Low"))

# --- ggplot2 panel functions -------------------------------------------------
col_map <- c(Linear = "black", High = "#2166ac", Low = "#d6604d")
shp_map <- c(Linear = 16, High = 17, Low = 15)

plot_three <- function(vl) {
  sub <- irf_df[irf_df$var == vl, ]
  ggplot(sub, aes(x = h, y = est, colour = model, shape = model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
    scale_colour_manual(values = col_map) + scale_shape_manual(values = shp_map) +
    scale_x_continuous(breaks = c(0,5,10,15,20)) +
    labs(x = "Horizon (quarters)", y = vl, title = "Three Models",
         colour = NULL, shape = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}
plot_linear <- function(vl) {
  sub <- irf_df[irf_df$var == vl & irf_df$model == "Linear", ]
  ggplot(sub, aes(x = h)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.6) +
    geom_line(aes(y = est), colour = "black", linewidth = 0.8) +
    geom_point(aes(y = est), colour = "black", size = 1.8) +
    scale_x_continuous(breaks = c(0,5,10,15,20)) +
    labs(x = "Horizon (quarters)", y = vl, title = "Linear Model") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
}
plot_state <- function(vl) {
  sH <- irf_df[irf_df$var == vl & irf_df$model == "High", ]
  sL <- irf_df[irf_df$var == vl & irf_df$model == "Low",  ]
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_ribbon(data = sH, aes(x=h, ymin=lower, ymax=upper), fill="#d1e5f0", alpha=0.6) +
    geom_ribbon(data = sL, aes(x=h, ymin=lower, ymax=upper), fill="#fddbc7", alpha=0.6) +
    geom_line(data=sH, aes(x=h,y=est), colour="#2166ac", linewidth=0.8) +
    geom_point(data=sH, aes(x=h,y=est), colour="#2166ac", shape=17, size=1.8) +
    geom_line(data=sL, aes(x=h,y=est), colour="#d6604d", linewidth=0.8, linetype="dashed") +
    geom_point(data=sL, aes(x=h,y=est), colour="#d6604d", shape=15, size=1.8) +
    scale_x_continuous(breaks = c(0,5,10,15,20)) +
    labs(x = "Horizon (quarters)", y = vl, title = "State-Dependent Models") +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
}

# --- Assemble and save -------------------------------------------------------
var_labs <- c("Employment","Real Wage","House Price")
panels   <- lapply(var_labs, function(vl) plot_three(vl) | plot_linear(vl) | plot_state(vl))
fig4 <- panels[[1]] / panels[[2]] / panels[[3]]
ggsave("Fig_4.png", fig4, width = 12, height = 10, dpi = 200)
cat("Saved: Fig_4.png\n")
