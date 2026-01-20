library(hexSticker)
library(ggplot2)
library(tweedie)

## ---- Generate a stylised Tweedie density ----

set.seed(1)

x <- seq(0, 4, length.out = 400)

## Tweedie parameters (chosen for visual shape, not pedagogy)
mu   <- 1.2
phi  <- 1
p    <- 1.05

y <- dtweedie(x, 
              mu = mu, 
              phi = phi, 
              power = p)

df <- data.frame(x = x, y = y)

density_plot <- ggplot(df, aes(x, y)) +
  geom_line(linewidth = 1.2, 
            colour = "#2C3E50") +
  theme_void() +
  coord_cartesian(expand = FALSE)

## ---- Create hex sticker ----

sticker(
  subplot    = density_plot,
  package    = "tweedie",
  p_size     = 18,
  p_color    = "#2C3E50",
  s_x        = 1,
  s_y        = 0.8,
  s_width    = 0.75,
  h_fill     = "white",
  h_color    = "#2C3E50",
  filename   = "tweedie_hex.png"
)
