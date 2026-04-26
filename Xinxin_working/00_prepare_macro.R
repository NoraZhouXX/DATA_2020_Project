# ==============================================================================
# 00_prepare_macro.R
# 从 ICPSR 的 All_data.xls (Macro_data sheet) 提取季度宏观变量,
# 输出 data/macro_q.csv 给后续 02_fig2.Rmd 等使用.
# ==============================================================================
# 依赖: readxl, dplyr, lubridate, readr
# 用法: 改下面两个路径, 然后 source() 整个文件即可.
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(readr)
})

# ------------------------------------------------------------------------------
# 路径配置
# ------------------------------------------------------------------------------
ICPSR_XLS  <- "/Users/norazhou/Desktop/DATA2020/Project/ICPSR_216361-V3.1/All_data.xls"
OUT_CSV    <- "/Users/norazhou/Desktop/DATA2020/Project/DATA_2020_Project/Xinxin_working/data/macro_q.csv"

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. 读 Excel 的 Macro_data sheet
# ------------------------------------------------------------------------------
raw <- read_excel(ICPSR_XLS, sheet = "Macro_data")

cat("Raw shape:", dim(raw), "\n")
cat("Columns  :", paste(names(raw), collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------
# 2. 把 "1975Q1" 字符串 -> Date (季度首月的 1 号)
#    Q1 -> Jan, Q2 -> Apr, Q3 -> Jul, Q4 -> Oct
# ------------------------------------------------------------------------------
q_to_date <- function(q) {
  yr  <- as.integer(substr(q, 1, 4))
  qtr <- as.integer(substr(q, 6, 6))
  mo  <- (qtr - 1) * 3 + 1
  as.Date(sprintf("%d-%02d-01", yr, mo))
}

# ------------------------------------------------------------------------------
# 3. 计算 log_gdp, log_cpi (按 paper Fig_2.prg: log(.)*100), 选列, 重命名
# ------------------------------------------------------------------------------
macro_q <- raw %>%
  mutate(
    date         = q_to_date(date),
    log_gdp      = log(RGDP) * 100,        # = log(rgdp) * 100  (paper convention)
    log_cpi      = log(PCE)  * 100,        # PCE deflator
    ffr          = FFR,
    rrshock      = RRSHOCK,                # Romer-Romer (Fig 5-7 用)
    shadow_rate  = SHADOW_RATE             # Wu-Xia (Fig 11 用)
  ) %>%
  select(date, log_gdp, log_cpi, ffr, rrshock, shadow_rate,
         employment   = EMPLOYMENT,        # B 做 Fig 4 用
         real_wage    = REAL_WAGE,
         house_price  = HOUSE_PRICE,
         r_investment = R_INVESTMET,       # paper 拼写就是少 N
         consumption  = CONSUMPTION) %>%
  arrange(date)

# ------------------------------------------------------------------------------
# 4. 保存
# ------------------------------------------------------------------------------
write_csv(macro_q, OUT_CSV)

# ------------------------------------------------------------------------------
# 5. Sanity check
# ------------------------------------------------------------------------------
cat("Saved:", OUT_CSV, "\n")
cat("Shape:", dim(macro_q), "\n")
cat("Date range:", as.character(min(macro_q$date)),
    "to", as.character(max(macro_q$date)), "\n\n")

# Fig 2 sample 检查
fig2_sample <- macro_q %>%
  filter(date >= as.Date("1981-01-01"),
         date <= as.Date("2007-10-01"))
cat("Fig 2 sample (1981Q1-2007Q4):", nrow(fig2_sample), "quarters (expect 108)\n\n")

# 缺失值
cat("Missing values:\n")
print(colSums(is.na(macro_q)))
cat("\n注意: rrshock 在 2008 之后是 NA, 这是 Romer-Romer 原始序列的特性, 正常.\n")
