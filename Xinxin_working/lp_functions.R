# ==============================================================================
# lp_functions.R
# Shared Local Projection toolkit for Lee & Ma (2025) replication
# Owner: Person A. Others: read-only, file issues for changes.
# ==============================================================================
# 依赖: dplyr, sandwich, lmtest
# 用法: source("shared/R/lp_functions.R")
#
# 主函数:
#   run_lp(data, dep, shock, controls, horizons, state_dummies = NULL,
#          hac_lag = NULL, ci_level = 0.68, trend = TRUE)
#
# 返回 data.frame:
#   horizon | spec       | point | se    | ci_low | ci_high
#   0       | linear     | -0.12 | 0.04  | -0.16  | -0.08
#   0       | high       | ...
#   0       | low        | ...
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(sandwich)
  library(lmtest)
})

# ------------------------------------------------------------------------------
# build_lp_data(): 给定原始季度数据 (date, log_gdp, log_cpi, ffr, shock, ...),
# 构造一份带滞后项 + lead 因变量 + trend 的扁平 data.frame.
#
#   data: 必须含 date 列 (Date) 和所有用到的变量
#   max_lag: 最大滞后阶数 (Lee & Ma 用 4)
#   max_h:   最大 horizon (Lee & Ma 用 20)
# ------------------------------------------------------------------------------
build_lp_data <- function(data, max_lag = 4, max_h = 20,
                          vars_to_lag = c("log_gdp", "log_cpi", "ffr"),
                          vars_to_lead = c("log_gdp", "log_cpi", "ffr")) {
  d <- data %>% arrange(date)

  # 加滞后
  for (v in vars_to_lag) {
    for (k in 0:max_lag) {
      if (k == 0) next
      d[[paste0(v, "_l", k)]] <- dplyr::lag(d[[v]], k)
    }
  }
  # 加 lead (因变量)
  for (v in vars_to_lead) {
    for (h in 0:max_h) {
      if (h == 0) next
      d[[paste0(v, "_h", h)]] <- dplyr::lead(d[[v]], h)
    }
  }
  # 加 trend
  d$trend  <- seq_len(nrow(d))
  d$trend2 <- d$trend^2
  d
}

# ------------------------------------------------------------------------------
# nw_se(): Newey-West HAC 标准误.
# Lee & Ma 用 lag = h + 1 (Jorda 2005 标准做法)
# ------------------------------------------------------------------------------
nw_se <- function(model, lag) {
  vc <- sandwich::NeweyWest(model, lag = lag, prewhite = FALSE, adjust = TRUE)
  lmtest::coeftest(model, vcov. = vc)
}

# ------------------------------------------------------------------------------
# run_lp(): 主函数.
#
#   data:     build_lp_data() 的输出
#   dep:      因变量名, e.g. "log_gdp" (函数自己加 _h{h})
#   shock:    冲击变量名, e.g. "ffr" 或 "rrshock"
#   controls: 控制变量名向量 (不含 shock 滞后, 函数会加)
#             e.g. c("log_gdp", "log_cpi") -> 自动加 _l1..l4
#   horizons: integer vector, e.g. 0:20
#   state_dummies: NULL = 线性 LP; list(high="D_H", low="D_L") = 状态 LP
#   hac_lag:  NULL = 用 h+1; 否则用指定值
#   trend:    TRUE 加 trend + trend^2
# ------------------------------------------------------------------------------
run_lp <- function(data, dep, shock, controls, horizons,
                   state_dummies = NULL, hac_lag = NULL,
                   ci_level = 0.68, trend = TRUE,
                   max_lag = 4, policy_var = NULL) {
  # policy_var: 在 state LP 中, 滞后政策利率的变量名 (默认 = shock).
  #   - Fig 2 (FFR shock):       shock="ffr",     policy_var=NULL -> "ffr"
  #   - Fig 5 (Romer-Romer):     shock="rrshock", policy_var="ffr"
  if (is.null(policy_var)) policy_var <- shock

  z <- qnorm(0.5 + ci_level / 2)  # 0.68 -> ~1; 0.95 -> ~1.96
  out <- list()

  for (h in horizons) {
    y <- if (h == 0) dep else paste0(dep, "_h", h)
    lag_for_nw <- if (is.null(hac_lag)) h + 1 else hac_lag

    # =========================================================================
    # 边界情况: h=0 且 dep == shock (e.g. FFR 对 FFR shock 的 impact response)
    #   - y_t = shock_t, 完美共线, R 会 drop shock 这列
    #   - 数学上 impact response = 1 (sign-flipped per paper convention = -1)
    # =========================================================================
    if (h == 0 && dep == shock) {
      if (is.null(state_dummies)) {
        out[[length(out) + 1]] <- data.frame(
          horizon = 0, spec = "linear",
          point = -1, se = 0, ci_low = -1, ci_high = -1
        )
      } else {
        for (sp in c("high", "low")) {
          out[[length(out) + 1]] <- data.frame(
            horizon = 0, spec = sp,
            point = -1, se = 0, ci_low = -1, ci_high = -1
          )
        }
      }
      next
    }

    # ----- 控制变量字符串 -----
    ctrl_terms <- c()
    # shock 同期 + 滞后 1..4
    ctrl_terms <- c(ctrl_terms, shock, paste0(shock, "_l", 1:max_lag))
    # 其他控制变量: 当期 + 滞后 1..4
    # NOTE: 跳过当期 dep (避免 y 出现在 RHS, h=0 时会完美共线)
    for (v in controls) {
      if (v != dep) ctrl_terms <- c(ctrl_terms, v)        # 当期 (跳过 dep)
      ctrl_terms <- c(ctrl_terms, paste0(v, "_l", 1:max_lag))
    }
    if (trend) ctrl_terms <- c(ctrl_terms, "trend", "trend2")

    # =========================================================================
    # CASE 1: 线性 LP
    # =========================================================================
    if (is.null(state_dummies)) {
      f <- as.formula(paste(y, "~", paste(ctrl_terms, collapse = " + ")))
      m <- lm(f, data = data)
      ct <- nw_se(m, lag = lag_for_nw)

      # shock 系数 (paper 取 -1 倍, 因为正向 shock 是紧缩)
      b  <- -1 * ct[shock, "Estimate"]
      se <-      ct[shock, "Std. Error"]

      out[[length(out) + 1]] <- data.frame(
        horizon = h, spec = "linear",
        point = b, se = se,
        ci_low = b - z * se, ci_high = b + z * se
      )

    # =========================================================================
    # CASE 2: 状态 LP (full interaction, EViews adlpm 风格)
    # =========================================================================
    } else {
      DH <- state_dummies$high
      DL <- state_dummies$low

      # 构造交互项: 每个控制变量 × DH, × DL
      build_interactions <- function(d_var) {
        terms <- c(d_var,
                   paste0(d_var, ":", shock))              # 关键: shock×state (当期)
        # 政策利率滞后 (Fig_2.prg 用 ffr 滞后, 不用 shock 滞后)
        for (k in 1:max_lag) {
          terms <- c(terms, paste0(d_var, ":", policy_var, "_l", k))
        }
        for (v in controls) {
          # 当期: 跳过 dep (避免 D_H:dep + D_L:dep = dep, 与 y 完美共线)
          if (v != dep) terms <- c(terms, paste0(d_var, ":", v))
          for (k in 1:max_lag) {
            terms <- c(terms, paste0(d_var, ":", v, "_l", k))
          }
        }
        terms
      }

      rhs <- c(build_interactions(DH), build_interactions(DL))
      if (trend) rhs <- c(rhs, "trend", "trend2")

      f <- as.formula(paste(y, "~ 0 +", paste(rhs, collapse = " + ")))
      m <- lm(f, data = data)
      ct <- nw_se(m, lag = lag_for_nw)

      # 提取系数: DH:shock 和 DL:shock
      # NOTE: R 有时会把交互项变量名顺序反过来 (e.g. "ffr:D_L" 而不是 "D_L:ffr"),
      # 所以两种命名都要试.
      get_interaction <- function(ct, var_state, var_other) {
        n1 <- paste0(var_state, ":", var_other)
        n2 <- paste0(var_other, ":", var_state)
        if (n1 %in% rownames(ct)) return(ct[n1, ])
        if (n2 %in% rownames(ct)) return(ct[n2, ])
        stop(sprintf("h=%d, dep=%s: 找不到系数 '%s' 或 '%s' (可能 R 把变量 drop 了).\n  Coef table: %s",
                     h, dep, n1, n2, paste(rownames(ct), collapse=", ")))
      }

      row_high <- get_interaction(ct, DH, shock)
      row_low  <- get_interaction(ct, DL, shock)

      bH  <- -1 * row_high["Estimate"];  seH <- row_high["Std. Error"]
      bL  <- -1 * row_low["Estimate"];   seL <- row_low["Std. Error"]

      out[[length(out) + 1]] <- data.frame(
        horizon = h, spec = "high",
        point = bH, se = seH,
        ci_low = bH - z * seH, ci_high = bH + z * seH
      )
      out[[length(out) + 1]] <- data.frame(
        horizon = h, spec = "low",
        point = bL, se = seL,
        ci_low = bL - z * seL, ci_high = bL + z * seL
      )
    }
  }

  bind_rows(out)
}

# ------------------------------------------------------------------------------
# Sanity check: 函数签名打印 (source 完会运行一次)
# ------------------------------------------------------------------------------
message("[lp_functions.R] loaded. Exports: build_lp_data(), run_lp(), nw_se()")
