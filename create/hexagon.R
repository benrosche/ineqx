library(hexSticker)
library(ggplot2)

# Create a stylized subplot representing variance decomposition
p <- ggplot() +
  # The "Between" gap highlight
  annotate(
    "segment",
    x = -1,
    xend = 1.2,
    y = 0.45,
    yend = 0.45,
    color = "#00F3FF",
    size = 0.6,
    linetype = "dotted"
  ) +
  # Area under the curve
  geom_area(data = df_within1, aes(x, y), fill = "#FF007A", alpha = 0.4) +
  geom_area(data = df_within2, aes(x, y), fill = "#00F3FF", alpha = 0.4) +
  # Distributions
  geom_line(data = df_within1, aes(x, y), color = "#FF007A", size = 0.7) +
  geom_line(data = df_within2, aes(x, y), color = "#00F3FF", size = 0.7) +
  # Dots at the end of distribution
  geom_point(
    data = subset(df_within1, x %in% seq(min(x), max(x), length.out = 5)),
    aes(x, y),
    color = "#FF007A",
    size = 0.7,
    alpha = 0.9
  ) +
  geom_point(
    data = subset(df_within2, x %in% seq(min(x), max(x), length.out = 5)),
    aes(x, y),
    color = "#00F3FF",
    size = 0.7,
    alpha = 0.9
  ) +
  theme_void() +
  theme_transparent()

# Generate the Hex Sticker
sticker(
  p,
  package = "ineqx",
  p_size = 26, # Package name size
  p_y = 1.5, # Position of name
  p_color = "#FFFFFF", # Text color
  s_x = 1,
  s_y = 0.85, # Subplot position
  s_width = 1.3,
  s_height = 0.8, # Subplot dimensions
  h_fill = "#141b2d", # Hexagon background color
  h_color = "#00F3FF", # Hexagon border color
  h_size = 2,
  filename = "ineqx-hexagon.jpg"
)
