################ Top transmitters and receivers throughout the time ################ 

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(lubridate)
library(patchwork)  

net_wide <- read_csv("outputs/directional_Net_monthly.csv") |>
  mutate(date = as.Date(date))

transmitters <- c("MI", "MN", "CA", "CO")
receivers    <- c("KS", "KY", "WV")

state_meta <- tibble(
  state  = c("MI",     "MN",      "CA",      "CO",      "KS",      "KY",      "WV"),
  label  = c("MI", "MN", "CA", "CO",
             "KS", "KY", "WV"),
  color  = c("#7b2d8b", "#1d4e89", "#2e86c1", "#85c1e9",
             "#e74c3c", "#922b21", "#641e16")
)

color_vec <- setNames(state_meta$color, state_meta$state)
label_vec <- setNames(state_meta$label, state_meta$state)

# Long format
make_long <- function(states) {
  net_wide |>
    select(date, all_of(states)) |>
    pivot_longer(-date, names_to = "state", values_to = "net") |>
    mutate(state = factor(state, levels = states))
}

tx_long <- make_long(transmitters)
rx_long <- make_long(receivers)

# NBER recessions
recessions <- tibble(
  start = as.Date(c("1981-07-01", "1990-07-01", "2001-03-01", "2007-12-01")),
  end   = as.Date(c("1982-11-30", "1991-03-31", "2001-11-30", "2009-06-30"))
) |>
  mutate(
    start = pmax(start, min(net_wide$date)),
    end   = pmin(end,   max(net_wide$date))
  )

recession_layer <- geom_rect(
  data = recessions,
  aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
  fill = "grey80", alpha = 0.5, inherit.aes = FALSE
)

hline_layer <- geom_hline(
  yintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed"
)

x_scale <- scale_x_date(
  limits      = c(as.Date("1981-01-01"), as.Date("2007-12-31")),
  date_breaks = "4 years",
  date_labels = "%Y",
  expand      = c(0.01, 0)
)

base_theme <- theme_minimal(base_size = 11) 

# Top panel: transmitters
p_tx <- ggplot() +
  recession_layer +
  hline_layer +
  geom_line(
    data = tx_long,
    aes(x = date, y = net, color = state),
    linewidth = 0.75, alpha = 0.9
  ) +
  scale_color_manual(
    values = color_vec[transmitters],
    labels = label_vec[transmitters],
    name   = NULL
  ) +
  x_scale  +
  labs(
    title = "Net Transmitter States",
    y     = "Net Connectedness (To − From)"
  ) +
  base_theme +
  theme(
    plot.title = element_text(face = "bold", size = 11, color = "#1d4e89"),
    plot.background = element_rect(fill = "#f0f4f9", color = NA)
  )

# Bottom panel: receivers
p_rx <- ggplot() +
  recession_layer +
  hline_layer +
  geom_line(
    data = rx_long,
    aes(x = date, y = net, color = state),
    linewidth = 0.75, alpha = 0.9
  ) +
  scale_color_manual(
    values = color_vec[receivers],
    labels = label_vec[receivers],
    name   = NULL
  ) +
  x_scale +
  labs(
    title = "Net Receiver States",
    y     = "Net Connectedness (To − From)"
  ) +
  base_theme +
  theme(
    plot.title = element_text(face = "bold", size = 11, color = "#922b21"),
    plot.background = element_rect(fill = "#fdf3f2", color = NA)
  )

# Combine with patchwork
p_combined <- p_tx / p_rx +
  plot_annotation(
    title    = "Net Directional Connectedness Over Time — Selected States",
    subtitle = "Positive = net transmitter  |  Negative = net receiver  |  Shaded = NBER recessions ",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave(
  "figures/Fig_E2.png",
  plot   = p_combined,
  width  = 13,
  height = 9,
  dpi    = 300,
  bg     = "white"
)