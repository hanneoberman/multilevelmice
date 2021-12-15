# plotting functions for incomplete data
# missingness indicator plot (monotone vs non-monotone)
plot_md_set <- function(...){
  miss_ind <- tibble(
  rownr = 5:1,
  # unit = rep("observed", 5),
  cluster = rep("observed", 5),
  X1 = rep("observed", 5),
  X2 = c(rep("observed", 2), rep("missing", 2), "observed"),
  X3 = c("missing", rep("observed", 2), "missing", "observed")
) %>% pivot_longer(cols = c(cluster, X1, X2, X3))

# add text 
miss_ind$text <- ""
# miss_ind[miss_ind$name == "unit", "text"] <- as.character(1:5)
miss_ind[miss_ind$name == "cluster", "text"] <- as.character(c(1,1,2,2,3))

# plot
miss_ind %>% 
  mutate(name = factor(name, levels = unique(name))) %>% 
  ggplot(aes(
    x = name,
    y = rownr,
    color = value,
    width = 0.8, 
    height = 0.8,
    label = text
  )) +
  geom_tile(fill = "white", size = 2) +
  geom_text(show.legend = FALSE) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = 5:1) +
  scale_color_manual(values = plot_col, name = "") +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    panel.grid.major = element_blank()) +
  labs(x = "", y = "")
}


# missing data pattern
plot_md_pat <- function(dat) {
  # get md pattern and store additional info
  pat <- mice::md.pattern(dat, plot = FALSE)
  vrb <- colnames(pat)[-ncol(pat)]
  vrb_labs <- abbreviate(vrb, minlength = 3)
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
      caption = stringr::str_wrap(paste("\n\nNote.", vrb_caption))
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
      plot.caption.position = "plot",
      axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(l = 10)),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = c("1" = mice:::mdc(1), "0" = mice:::mdc(2)))
  # output
  return(p)
}

# plot conditional distributions
plot_conditional <- function(data, x, z, cluster, ...){
  ggplot(data, aes(x = get(x), color = factor(is.na(get(z)), labels = c("missing", "observed")))) +
    geom_density(data = data %>% filter(!is.na(get(z))), aes(x = get(x), group = get(cluster)), alpha = 0.01, fill = plot_col[2], color = NA) +
    geom_density(data = data %>% filter(is.na(get(z))), aes(x = get(x), group = get(cluster)), alpha = 0.01, fill = plot_col[1], color = NA) +
    geom_density() +
    scale_color_manual(values = plot_col, name = z) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(x = x)
}


#' Title Create a histogram/bar plot of a variable conditional on missingness in another variable
#'
#' @param dat An incomplete dataset of class dataframe
#' @param x A variable to plot
#' @param z A second variable to plot
#'
#' @return A ggplot object with two facets (variable z observed vs. missing)
#' @export
plot_NA_cond <- function(dat, x, z, bins = NULL) {
  # define graphing elements to add to plot
  if (is.numeric(dat[[x]])) {
    # for continuous variables
    geom <- ggplot2::geom_histogram(fill = "white", binwidth = bins)
  } else {
    # for categorical variables
    geom <- ggplot2::geom_bar(fill = "white")
  }
  # create facet labels
  facet_labs <- c(paste(z, "observed"), paste(z, "missing")) %>%
    setNames(c("observed", "missing"))
  # preprocess the data
  d <- dat[!is.na(dat[[x]]),] %>%
    dplyr::mutate(conditional = factor(
      is.na(.data[[z]]),
      levels = c(FALSE, TRUE),
      labels = c("observed", "missing")
    )) 
  
  # plot
  p <- d %>% ggplot2::ggplot(ggplot2::aes(x = .data[[x]],
                                          color = conditional,
                                          fill = "white")) +
    geom +
    # split by conditional variable
    ggplot2::facet_wrap(
      ~ conditional,
      ncol = 1,
      scales = "free_y",
      labeller = ggplot2::labeller(conditional = facet_labs)
    ) +
    # style
    ggplot2::scale_color_manual(
        values = c(
          "observed" = mice:::mdc(1),
          "missing" = mice:::mdc(2),
          "imputed" = mice:::mdc(2)
        )
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "observed" = mice:::mdc(1),
          "missing" = mice:::mdc(2),
          "imputed" = mice:::mdc(2)
        )
      ) +
      ggplot2::scale_size_manual(values = c(
        "observed" = 1,
        "imputed" = 0.5
      )) +
      ggplot2::labs(
        color = NULL,
        fill = NULL,
        size = NULL
      ) +
    theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(size = 0.5)
    ) +
    labs(title = paste("Distribution of", x, "given missingness in", z))
  # output
  return(p)
}
