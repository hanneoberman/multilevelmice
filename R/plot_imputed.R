# plotting functions for imputed data

# basic pre-processing of imputations
plot_imps <- function(imp, x, y = NULL) {
  # parse inputs
  if (is.null(y)) {
    y <- x
  }
  # combine observed and imputed data
  xy_obs <- imp$data %>%
    cbind(datapoint = "observed",
          .imp = 0,
          .id = 1:nrow(.),
          .) %>%
    .[!is.na(imp$data[[x]]) &
        !is.na(imp$data[[y]]),]
  xy_imps <- imp %>%
    mice::complete("long") %>%
    cbind(datapoint = "imputed", .) %>%
    .[.$.id %nin% xy_obs$.id, ]
  xy_dat <- rbind(xy_obs, xy_imps) %>%
    dplyr::mutate(datapoint = factor(datapoint, levels = c("observed", "imputed")))
  # initialize plot
  p <- xy_dat %>%
    ggplot2::ggplot() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::scale_color_manual(
      values = c(
      "observed" = mice:::mdc(1),
      # "missing" = mice:::mdc(2),
      "imputed" = mice:::mdc(2)
      )
  ) +
    ggplot2::scale_fill_manual(
    values = c(
      "observed" = mice:::mdc(1),
      # "missing" = mice:::mdc(2),
      "imputed" = mice:::mdc(2)
    ))
  # output
  return(p)
}

# stripplot
plot_strip <- function(imp, x) {
  # plot individual values (stripplot)
  p <- imp %>% plot_imps(x) +
    ggplot2::geom_jitter(
      ggplot2::aes(
        y = factor(.imp, ordered = TRUE),
        x = .data[[x]],
        color = datapoint
      ),
      height = 0.25,
      width = 0
    ) +
    scale_y_discrete(limits = rev) +
    ggplot2::ylab("Imputation (0 = observed data)")
  # output
  return(p)
}

# boxplot
plot_box <- function(imp, x, cluster = NULL, strip = FALSE) {
  # plot
  if (strip){
    p <- imp %>% plot_imps(x) +
      ggplot2::geom_jitter(
        ggplot2::aes(
          y = factor(.imp, ordered = TRUE),
          x = .data[[x]],
          color = datapoint
        ),
        height = 0.25,
        width = 0,
        alpha = 0.1,
        stroke = 0
      )
  } else {
    p <- imp %>% plot_imps(x)
  }
  # if (is.numeric(imp$data[[x]])) {
    # plot box and whiskers for numeric variables
    p <- p +
      ggplot2::geom_boxplot(ggplot2::aes(
        y = as.factor(.imp),
        x = .data[[x]],
        color = datapoint
      ),
      width = 0.5,
      alpha = 0.5,
      outlier.shape = NA) +
      scale_y_discrete(limits = rev) +
      ggplot2::ylab("Imputation (0 = observed data)")
  # } else {
  #   # plot faceted barplot for categorical variables
  #   p <- imp %>% plot_imps(x) +
  #     ggplot2::geom_bar(ggplot2::aes(y = .data[[x]],
  #                                    color = datapoint),
  #                       width = 0.5,
  #                       fill = "white") +
  #     ggplot2::facet_wrap(
  #       ~ .imp,
  #       scales = "free_x",
  #       ncol = 2,
  #       labeller = ggplot2::labeller(.imp = c(
  #         "Observed data", paste("Imputation", 1:imp$m)
  #       ) %>% setNames(0:imp$m))
  #     )
  # }
  if(!is.null(cluster)){
    p <- p + facet_wrap(get(cluster))
    }
    # output
  return(p)
}

# new plotting function for imputed data
plot_mids <- function(imp, vrb, type = c("bwplot", "stripplot", "densityplot")) {
  
  # process mids object for plotting
  obs <- do.call(rbind, replicate(imp$m, !imp$where, simplify = FALSE)) %>%
    rbind(matrix(FALSE, nrow(imp$data), ncol(imp$data)), .) %>% 
    cbind(.imp = FALSE, .id = FALSE, .) 
  comp <- mice::complete(imp, "long", include = TRUE) %>% 
    mutate(.imp = factor(.imp, ordered = TRUE))
  comp[obs] <- NA
  
  # bwplot or stripplot
  if (type == "stripplot" | type == "bwplot"){
    p <- ggplot(comp, aes(x = .imp, y = get(vrb), color = ifelse(.imp < 1, "Observed", "Imputed"))) +
    labs(x = "Imputation",
         y = vrb,
         color = "",
         fill = "") 
    if(type == "stripplot"){
      p <- p + geom_point(position = position_jitter(width = 0.25, height = 0))
    }
    p <- p + geom_boxplot(size = 1,
                 width = 0.5,
                 alpha = 0.5,
                 outlier.shape = NA)
  }
  
  # densityplot
  if (type == "densityplot"){
  p <- ggplot(comp, aes(x = get(vrb), group = .imp, color = ifelse(.imp < 1, "Observed", "Imputed"))) +
    geom_density() +
    labs(x = vrb,
         color = "",
         fill = "")
  }
  
  # # xyplot
  # ggplot(comp, aes(x = extrav, y = popular, color = ifelse(.imp < 1, "Observed", "Imputed"))) +
  #   geom_point() +
  #   geom_point(aes())
  #   labs(color = "")
  
  # output
  return(p)
}
 
# # test
# plot_mids(imp_ignored, "popular")

