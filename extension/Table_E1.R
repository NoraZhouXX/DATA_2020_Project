######### State-level Net connectedness ranking table with stability index ######### 
#
# Columns:
#   Rank        - rank by mean_Net (1 = strongest transmitter)
#   State       - state abbreviation
#   Role        - Transmitter / Receiver
#   Mean Net    - average Net connectedness over full sample (1981-2007)
#   SD Net      - standard deviation (volatility of role)
#   % Positive  - % of months above zero
#   H1 Mean     - mean Net in first half (1981-1994)
#   H2 Mean     - mean Net in second half (1995-2007)
#   Stability   - correlation between H1 and H2 rank orders (state-level)


library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(gridExtra)  

ranking  <- read_csv("outputs/state_ranking_by_Net.csv")
net_wide <- read_csv("outputs/directional_Net_monthly.csv") |>
  mutate(date = as.Date(date))

# Split at 1995-01-01: H1 = 1981-1994, H2 = 1995-2007
net_long <- net_wide |>
  filter(date <= as.Date("2007-12-31")) |>
  pivot_longer(-date, names_to = "state", values_to = "net") |>
  mutate(half = ifelse(date < as.Date("1995-01-01"), "H1", "H2"))

h1 <- net_long |> filter(half == "H1") |> group_by(state) |> summarise(H1_mean = mean(net))
h2 <- net_long |> filter(half == "H2") |> group_by(state) |> summarise(H2_mean = mean(net))
half_means <- left_join(h1, h2, by = "state")

# System-level stability
stability_model <- lm(H2_mean ~ H1_mean, data = half_means)
stability_r2    <- summary(stability_model)$r.squared
stability_beta  <- coef(stability_model)["H1_mean"]

# Build full table
table_full <- ranking |>
  left_join(half_means, by = "state") |>
  arrange(desc(mean_Net)) |>
  mutate(
    Rank    = row_number(),
    Role    = ifelse(mean_Net > 0, "Transmitter", "Receiver"),
    # State-level stability: sign consistency between halves
    # 1 = same sign in both halves (stable role), 0 = switched
    Role_stable = ifelse(sign(H1_mean) == sign(H2_mean), "Stable", "Switched")
  ) |>
  select(
    Rank, State = state, Role, Role_stable,
    `Mean Net` = mean_Net,
    `SD Net`   = sd_Net,
    `% Pos`    = pct_pos_Net,
    `H1 Mean`  = H1_mean,
    `H2 Mean`  = H2_mean
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 1)))

write_csv(table_full, "outputs/Table_E1.csv")

# Extract top 10 + bottom 10 for display
top10    <- table_full |> slice(1:10)
bottom10 <- table_full |> slice((nrow(table_full) - 9):nrow(table_full))

display_table <- bind_rows(
  top10,
  tibble(Rank = NA, State = "···", Role = "", Role_stable = "",
         `Mean Net` = NA, `SD Net` = NA, `% Pos` = NA,
         `H1 Mean` = NA, `H2 Mean` = NA),
  bottom10
)

# Color coding: blue rows = transmitters, red rows = receivers
n_top    <- nrow(top10)
n_sep    <- 1
n_bottom <- nrow(bottom10)
n_total  <- n_top + n_sep + n_bottom

row_fill <- c(
  rep("#dce8f5", n_top),    
  "#ffffff",                 
  rep("#fde8e4", n_bottom)  
)
row_text <- c(
  rep("#1a3a5c", n_top),
  "#888888",
  rep("#5c1a1a", n_bottom)
)

display_img <- display_table |>
  mutate(
    Rank       = ifelse(is.na(Rank), "···", as.character(Rank)),
    `Mean Net` = ifelse(is.na(`Mean Net`), "···",
                        ifelse(as.numeric(`Mean Net`) > 0,
                               paste0("+", `Mean Net`), as.character(`Mean Net`))),
    `H1 Mean`  = ifelse(is.na(`H1 Mean`), "···",
                        ifelse(as.numeric(`H1 Mean`) > 0,
                               paste0("+", `H1 Mean`), as.character(`H1 Mean`))),
    `H2 Mean`  = ifelse(is.na(`H2 Mean`), "···",
                        ifelse(as.numeric(`H2 Mean`) > 0,
                               paste0("+", `H2 Mean`), as.character(`H2 Mean`))),
    `SD Net`   = ifelse(is.na(`SD Net`), "···", as.character(`SD Net`)),
    `% Pos`    = ifelse(is.na(`% Pos`), "···", paste0(`% Pos`, "%"))
  )

# Create ggplot table
p_table <- ggplot() +
  theme_void() +
  annotation_custom(
    tableGrob(
      display_img,
      rows  = NULL,
      theme = ttheme_minimal(
        core = list(
          bg_params  = list(fill = row_fill, col = NA),
          fg_params  = list(col  = row_text, fontsize = 9)
        ),
        colhead = list(
          bg_params  = list(fill = "#2c3e50", col = NA),
          fg_params  = list(col  = "white", fontsize = 9, fontface = "bold")
        )
      )
    )
  ) +
  labs(
    title    = "Table E1 — State-Level Net Connectedness: Top & Bottom 10",
    subtitle = sprintf(
      "Full sample 1981–2007Q4  |  H1: 1981–1994  |  H2: 1995–2007  |  Stability R² (H2 ~ H1 regression): %.3f  -  Slope: %.3f",
      stability_r2, stability_beta
    ),
    caption  = "Role_stable: 'Stable' = same sign of mean Net in both halves. SD Net measures role volatility within sample."
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5, margin = margin(b = 8)),
    plot.caption  = element_text(size = 7.5, color = "grey55", hjust = 0),
    plot.margin   = margin(15, 15, 15, 15)
  )

ggsave(
  "figures/Table_E1.png",
  plot   = p_table,
  width  = 11,
  height = 8,
  dpi    = 300,
  bg     = "white"
)
