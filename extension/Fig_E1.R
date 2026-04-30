######## U.S. choropleth of mean Net directional connectedness (L2 Extension) ######## 

library(dplyr)
library(ggplot2)
library(readr)
library(usmap)       
library(ggrepel)   
library(sf)

ranking <- read_csv("outputs/state_ranking_by_Net.csv")

ranking <- ranking |>
  rename(values = mean_Net)

# Identify top 5 transmitters and bottom 5 receivers
top5    <- ranking |> slice_max(values, n = 5)   # highest mean_Net
bottom5 <- ranking |> slice_min(values, n = 5)   # lowest mean_Net
labeled <- bind_rows(top5, bottom5)

# Adding DC as a state 
state_centers <- data.frame(
  state = c(state.abb, "DC"),
  lon   = c(state.center$x, -77.0369),
  lat   = c(state.center$y,  38.9072)
)
# Project the coordinates into usmap's coordinate system
transformed <- usmap::usmap_transform(state_centers, input_names = c("lon", "lat"))

# Getting the dataframe with state, lon, lat
coords    <- st_coordinates(transformed)
centroids <- transformed |>
  st_drop_geometry() |>
  mutate(x = coords[, 1], y = coords[, 2]) |>
  select(state, x, y)

# Labeling the top/bottom states according to its coordinates
labeled_coords <- labeled |>
  left_join(centroids, by = "state") |>
  mutate(
    label_color = ifelse(values > 0, "transmitter", "receiver"),
    label = paste0(state, "\n", round(values, 1))
  )

# Symmetric color scale centered at 0 
max_abs <- max(abs(ranking$values))  
limit   <- ceiling(max_abs / 5) * 5  # Round up to the nearest multiple of 5

# Build map
p <- plot_usmap(
  data        = ranking,
  values      = "values",
  color       = "white",    # state border color
  linewidth   = 0.3
) +
  scale_fill_gradientn(
    colours = c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac"),
    limits  = c(-limit, limit),
    breaks  = c(-limit, -limit/2, 0, limit/2, limit),
    name    = "Mean Net\nConnectedness"
  ) +
  geom_label(
    data = labeled_coords |> filter(label_color == "transmitter"),
    aes(x = x, y = y, label = label),
    size      = 2.8,
    fontface  = "bold",
    fill      = alpha("#2166ac", 0.85),
    color     = "white",
    linewidth = 0.2
  ) +
  geom_label(
    data = labeled_coords |> filter(label_color == "receiver"),
    aes(x = x, y = y, label = label),
    size      = 2.8,
    fontface  = "bold",
    fill      = alpha("#b2182b", 0.85),
    color     = "white",
    linewidth = 0.2
  ) +
  labs(
    title    = "State-Level Net Directional Connectedness in U.S. Housing Markets",
    subtitle = "Mean Net connectedness, rolling 120-month windows, 1981Q1–2007Q4\nBlue = net transmitter  |  Red = net receiver  |  Labeled: top 5 and bottom 5 states"
  ) +
  theme(
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 9,  color = "grey40", hjust = 0.5, lineheight = 1.3)
  )

ggsave(
  "figures/Fig_E1.png",
  plot   = p,
  width  = 11,
  height = 7,
  dpi    = 300,
  bg     = "white"
)