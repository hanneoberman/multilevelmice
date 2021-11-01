# missing data pattern
md_pat <- function(dat) {
  # get md pattern and store additional info
  pat <- mice::md.pattern(dat, plot = FALSE)
  vrb <- colnames(pat)[-ncol(pat)]
  vrb_labs <- abbreviate(vrb)
  vrb_caption <- purrr::map2_chr(vrb_labs, vrb, ~{paste(.x, " = ", .y)}) %>% paste(collapse = ", ")
  colnames(pat) <- c(vrb, "NA_per_pat")
  pat_freq <- as.numeric(rownames(pat))[-nrow(pat)]
  NA_per_pat <- as.numeric(pat[, ncol(pat)])[-nrow(pat)]
  NA_per_vrb <- as.numeric(pat[nrow(pat),])[-ncol(pat)]
  NA_total <- pat[nrow(pat), ncol(pat)]
  # make the pattern tidy
  if (NA_total == 0) {
    pat_long <-
      data.frame(
        NA_per_vrb = 0,
        NA_per_pat = 0,
        pat_freq = pat_freq,
        pat_nr = 1,
        obs = 1
      ) %>% cbind(vrb = vrb)
  } 
  if (NA_total != 0) {
    pat_long <- pat[-nrow(pat),] %>%
    cbind(., pat_freq, pat_nr = 1:nrow(.)) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = all_of(vrb),
                        names_to = "vrb",
                        values_to = "obs") %>%
    cbind(NA_per_vrb, .)
  }
  # plot the md pattern
  p <- pat_long %>%
    ggplot2::ggplot() +
    ggplot2::geom_tile(ggplot2::aes(
      x = vrb,
      y = pat_nr,
      fill = as.factor(obs),
      group = NA_per_pat
    ),
    color = "black") +
    # set axes
    ggplot2::scale_x_discrete(
      limits = vrb,
      position = "bottom",
      labels = as.character(NA_per_vrb),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_reverse(
      breaks = 1:max(pat_long$pat_nr),
      labels = as.character(pat_freq),
      expand = c(0, 0),
      sec.axis = ggplot2::dup_axis(labels = as.character(NA_per_pat),
                                   name = "Number of missing entries per pattern")
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    # add labels
    ggplot2::labs(
      x = "Number of missing entries per variable",
      y = "Pattern frequency",
      title = "Missing data pattern",
      subtitle = paste0("Total number of missing entries: ",
                        NA_total,
                        "\n\n"),
      caption = paste("\n Note.", vrb_caption)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = vrb,
        y = -Inf,
        label = vrb_labs
      ),
      data = pat_long[1:length(vrb), ],
      vjust = -0.5
    ) +
    # add styling
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(
        t = 20,
        l = 10,
        b = 10,
        r = 10,
        unit = "pt"
      ),
      axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(l = 10)),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = c("1" = mice:::mdc(1), "0" = mice:::mdc(2)))
  # output
  return(p)
}
