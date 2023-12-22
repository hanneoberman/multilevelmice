# some functions used in the manuscript

#' Visualize missing data (syst. vs. spor.)
#'
#' @return An object of class `ggplot`
plot_na <- function() {
  dat <- expand.grid(rows = 1:7, cols = 1:6) |>
    cbind(
      text = c("1", "1", "2", "2", "3", "", "N", rep("", 35)),
      miss = c(rep("", 16), "NA", "NA", "", "", "", "NA", "", "", "NA", rep("", 17))
    )
  ggplot(dat, aes(x = cols, y = rows)) +
    geom_tile(fill = "white",
              color = "black",
              linewidth = 0.5) +
    geom_text(aes(label = text), color = "black", size = 3) +
    geom_text(
      aes(label = miss),
      color = mice:::mdc(2),
      family = "mono",
      fontface = "bold"
    ) +
    scale_x_continuous(
      breaks = 1:6,
      labels = c(
        "cluster",
        expression(X[1]),
        expression(X[2]),
        expression(X[3]),
        "...",
        expression(X[p])
      ),
      name = NULL,
      position = "top"
    ) +
    scale_y_continuous(
      breaks = 1:7,
      labels = c(1:5, "...", "n"),
      name = NULL,
      trans = "reverse"
    ) +
    # coord_cartesian(expand = c(0,0)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt"))
}
