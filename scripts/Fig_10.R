# =============================================================================
# Fig_10.R
# Robustness: HP filter smoothing κ = 3, 5, 10  (λ = κ × 160)
#   κ=10 → λ=1600 (baseline);  κ=5 → λ=800;  κ=3 → λ=480
#
# Layout: 3 rows (κ=3 / κ=5 / κ=10) × 3 columns (GDP / PCE / FFR)
# Each panel: Linear (black solid) + High (blue dashed + CI) + Low (red dotted + CI)
#
# Lab methods: lm() (Lab 4), sandwich::NeweyWest(), ggplot2 (Lab 2/4/9)
# INPUTS:  All_data.xls, Ct_monthly.csv
# OUTPUT:  Fig_10.png
# =============================================================================
rm(list = ls())
library(readxl); library(sandwich); library(ggplot2); library(patchwork)

# ---------------------------------------------------------------------------
# 1. Load macro data
# ---------------------------------------------------------------------------
macro   <- read_excel("../data/All_data.xls", sheet = "Macro_data")
n       <- nrow(macro)
yr      <- seq(1975.0, by = 0.25, length.out = n)
log_GDP <- log(as.numeric(macro$RGDP)) * 100
log_CPI <- log(as.numeric(macro$PCE))  * 100
FFR     <- as.numeric(macro$FFR)

smpl   <- yr >= 1981.0 & yr <= 2007.75
T_smpl <- sum(smpl)
bw_fix <- floor(0.75 * T_smpl^(1/3))

# ---------------------------------------------------------------------------
# 2. Load monthly Ct → quarterly, then HP filter for each κ
# ---------------------------------------------------------------------------
ct_raw <- read.csv("../data/Ct_monthly.csv", header = TRUE)[, 2]
T_q    <- length(ct_raw) / 3
ct_q   <- sapply(1:T_q, function(j) mean(ct_raw[(3*(j-1)+1):(3*j)]))
yr_q   <- seq(1981.0, by = 0.25, length.out = T_q)

hp_filter <- function(y, lambda) {
  nn <- length(y)
  D  <- diff(diag(nn), differences = 2)
  trend <- solve(diag(nn) + lambda * t(D) %*% D, y)
  list(trend = as.numeric(trend), cycle = as.numeric(y - trend))
}

# κ × 160 = λ  (κ=10 → λ=1600 = quarterly HP baseline)
kappa_lambda <- c(`3` = 480, `5` = 800, `10` = 1600)

get_dummies <- function(lambda) {
  hp     <- hp_filter(ct_q, lambda)
  p1_q   <- c(NA, as.integer(hp$cycle >= 0)[-T_q])
  n1_q   <- c(NA, as.integer(hp$cycle <  0)[-T_q])
  idx    <- match(round(yr_q, 4), round(yr, 4))
  p1     <- rep(NA, n);  p1[idx] <- p1_q
  n1     <- rep(NA, n);  n1[idx] <- n1_q
  list(p1 = p1, n1 = n1)
}

# ---------------------------------------------------------------------------
# 3. LP helpers
# ---------------------------------------------------------------------------
lag_vec <- function(x, l) c(rep(NA, l), x[seq_len(n - l)])

run_lp_h <- function(Y_lead, p1, n1, bw) {
  df <- data.frame(Y = Y_lead, SHOCK = FFR, p1 = p1, n1 = n1,
                   trend = seq_len(n), trend2 = seq_len(n)^2)
  df$p1_SHOCK <- p1 * FFR;  df$n1_SHOCK <- n1 * FFR
  for (l in 0:4) {
    df[[paste0("gdp_l",l)]]    <- lag_vec(log_GDP, l)
    df[[paste0("cpi_l",l)]]    <- lag_vec(log_CPI, l)
    df[[paste0("p1_gdp_l",l)]] <- p1 * df[[paste0("gdp_l",l)]]
    df[[paste0("n1_gdp_l",l)]] <- n1 * df[[paste0("gdp_l",l)]]
    df[[paste0("p1_cpi_l",l)]] <- p1 * df[[paste0("cpi_l",l)]]
    df[[paste0("n1_cpi_l",l)]] <- n1 * df[[paste0("cpi_l",l)]]
  }
  for (l in 1:4) {
    df[[paste0("ffr_l",l)]]    <- lag_vec(FFR, l)
    df[[paste0("p1_ffr_l",l)]] <- p1 * df[[paste0("ffr_l",l)]]
    df[[paste0("n1_ffr_l",l)]] <- n1 * df[[paste0("ffr_l",l)]]
  }
  df <- df[smpl, ];  df <- df[complete.cases(df), ]

  # Linear model (Lab 4: lm)
  lin_rhs <- c("SHOCK", paste0("gdp_l",0:4), paste0("cpi_l",0:4),
               paste0("ffr_l",1:4), "trend", "trend2")
  m_lin <- lm(as.formula(paste("Y ~", paste(lin_rhs, collapse=" + "))), data=df)
  V_lin <- NeweyWest(m_lin, lag=bw, prewhite=FALSE)

  # State-dependent model (Lab 9 DiD style)
  sd_rhs <- c("0", "p1", "n1", "p1_SHOCK","n1_SHOCK",
              paste0("p1_gdp_l",0:4), paste0("n1_gdp_l",0:4),
              paste0("p1_cpi_l",0:4), paste0("n1_cpi_l",0:4),
              paste0("p1_ffr_l",1:4), paste0("n1_ffr_l",1:4),
              "trend","trend2")
  m_sd <- lm(as.formula(paste("Y ~", paste(sd_rhs, collapse=" + "))), data=df)
  V_sd <- NeweyWest(m_sd, lag=bw, prewhite=FALSE)

  list(b_lin = -coef(m_lin)["SHOCK"],          s_lin = sqrt(V_lin["SHOCK","SHOCK"]),
       b_H   = -coef(m_sd)["p1_SHOCK"],        s_H   = sqrt(V_sd["p1_SHOCK","p1_SHOCK"]),
       b_L   = -coef(m_sd)["n1_SHOCK"],        s_L   = sqrt(V_sd["n1_SHOCK","n1_SHOCK"]))
}

# ---------------------------------------------------------------------------
# 4. Run LP for all κ values, all horizons, all variables
# ---------------------------------------------------------------------------
H_max   <- 20;  h_seq <- 0:H_max
vars_ls <- list(GDP = log_GDP, PCE = log_CPI, FFR = FFR)

cat("Running LP for κ = 3, 5, 10 ...\n")
all_res <- lapply(names(kappa_lambda), function(kname) {
  cat("  κ =", kname, "(λ =", kappa_lambda[kname], ")\n")
  dum <- get_dummies(kappa_lambda[kname])
  p1  <- dum$p1;  n1 <- dum$n1
  lapply(seq_along(h_seq), function(hi) {
    h <- h_seq[hi]
    lapply(names(vars_ls), function(vname) {
      y_v <- vars_ls[[vname]]
      run_lp_h(c(y_v[(h+1):n], rep(NA, h)), p1, n1, bw_fix)
    }) |> setNames(names(vars_ls))
  })
}) |> setNames(names(kappa_lambda))

# ---------------------------------------------------------------------------
# 5. Compile tidy IRF data frame
# ---------------------------------------------------------------------------
irf_df <- do.call(rbind, lapply(names(kappa_lambda), function(kname) {
  do.call(rbind, lapply(seq_along(h_seq), function(hi) {
    h   <- h_seq[hi];  res <- all_res[[kname]][[hi]]
    do.call(rbind, lapply(names(vars_ls), function(vname) {
      r <- res[[vname]]
      data.frame(
        h = h, var = vname, kappa = kname,
        model = c("Linear","High","Low"),
        est   = c(r$b_lin, r$b_H,          r$b_L),
        lower = c(r$b_lin - r$s_lin, r$b_H - r$s_H, r$b_L - r$s_L),
        upper = c(r$b_lin + r$s_lin, r$b_H + r$s_H, r$b_L + r$s_L),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

irf_df$var   <- factor(irf_df$var,   levels = c("GDP","PCE","FFR"))
irf_df$kappa <- factor(irf_df$kappa, levels = c("3","5","10"))
irf_df$model <- factor(irf_df$model, levels = c("Linear","High","Low"))

# ---------------------------------------------------------------------------
# 6. ggplot2: one panel per (κ, variable)
#    Linear = black solid (no CI);  High = blue dashed + CI;  Low = red dotted + CI
# ---------------------------------------------------------------------------
make_panel <- function(kname, vname) {
  sub  <- irf_df[irf_df$kappa == kname & irf_df$var == vname, ]
  lin  <- sub[sub$model == "Linear", ]
  high <- sub[sub$model == "High",   ]
  low  <- sub[sub$model == "Low",    ]

  ggplot() +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    # High CI band (blue shading)
    geom_ribbon(data = high, aes(x = h, ymin = lower, ymax = upper),
                fill = "#2166ac", alpha = 0.2) +
    # Low CI band (red, shown as dashed boundary lines only)
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
    scale_x_continuous(breaks = c(0,5,10,15,20)) +
    labs(x = "Horizon (Quarter)", y = NULL, title = vname) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
}

# ---------------------------------------------------------------------------
# 7. Assemble 3×3 grid: rows = κ, cols = GDP/PCE/FFR
#    Row subtitles: (a) κ = 3  /  (b) κ = 5  /  (c) κ = 10
# ---------------------------------------------------------------------------
kappa_labels <- c("3","5","10")
row_letters  <- c("(a)", "(b)", "(c)")
var_labels   <- c("GDP","PCE","FFR")

rows <- lapply(seq_along(kappa_labels), function(ki) {
  k      <- kappa_labels[ki]
  letter <- row_letters[ki]

  p_gdp <- make_panel(k, "GDP") +
    labs(tag = letter) +
    theme(plot.tag = element_text(face = "bold", size = 13))

  # κ label on the middle (PCE) panel as subtitle, centered
  p_pce <- make_panel(k, "PCE") +
    labs(subtitle = bquote(italic(kappa) == .(as.integer(k)))) +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 12,
                                       margin = margin(t = 0, b = 3)))

  p_ffr <- make_panel(k, "FFR")

  p_gdp | p_pce | p_ffr
})

fig10 <- rows[[1]] / rows[[2]] / rows[[3]]

ggsave("../output/Fig_10.png", fig10, width = 12, height = 10, dpi = 200)
cat("Saved: Fig_10.png\n")
