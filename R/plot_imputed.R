# # plotting functions for imputed data
# 
# # basic pre-processing of imputations
# plot_imps <- function(imp, x, y = NULL) {
#   # parse inputs
#   if (is.null(y)) {
#     y <- x
#   }
#   # combine observed and imputed data
#   xy_obs <- imp$data %>%
#     cbind(datapoint = "observed",
#           .imp = 0,
#           .id = 1:nrow(.),
#           .) %>%
#     .[!is.na(imp$data[[x]]) &
#         !is.na(imp$data[[y]]),]
#   xy_imps <- imp %>%
#     mice::complete("long") %>%
#     cbind(datapoint = "imputed", .) %>%
#     .[.$.id %nin% xy_obs$.id, ]
#   xy_dat <- rbind(xy_obs, xy_imps) %>%
#     dplyr::mutate(datapoint = factor(datapoint, levels = c("observed", "imputed")))
#   # initialize plot
#   p <- xy_dat %>%
#     ggplot2::ggplot() +
#     ggplot2::theme_classic() +
#     ggplot2::theme(legend.position = "top") +
#     ggplot2::scale_color_manual(
#       values = c(
#       "observed" = mice:::mdc(1),
#       # "missing" = mice:::mdc(2),
#       "imputed" = mice:::mdc(2)
#       )
#   ) +
#     ggplot2::scale_fill_manual(
#     values = c(
#       "observed" = mice:::mdc(1),
#       # "missing" = mice:::mdc(2),
#       "imputed" = mice:::mdc(2)
#     ))
#   # output
#   return(p)
# }
# 
# # stripplot
# plot_strip <- function(imp, x) {
#   # plot individual values (stripplot)
#   p <- imp %>% plot_imps(x) +
#     ggplot2::geom_jitter(
#       ggplot2::aes(
#         y = factor(.imp, ordered = TRUE),
#         x = .data[[x]],
#         color = datapoint
#       ),
#       height = 0.25,
#       width = 0
#     ) +
#     scale_y_discrete(limits = rev) +
#     ggplot2::ylab("Imputation (0 = observed data)")
#   # output
#   return(p)
# }
# 
# # boxplot
# plot_box <- function(imp, x, cluster = NULL, strip = FALSE) {
#   # plot
#   if (strip){
#     p <- imp %>% plot_imps(x) +
#       ggplot2::geom_jitter(
#         ggplot2::aes(
#           y = factor(.imp, ordered = TRUE),
#           x = .data[[x]],
#           color = datapoint
#         ),
#         height = 0.25,
#         width = 0,
#         alpha = 0.1,
#         stroke = 0
#       )
#   } else {
#     p <- imp %>% plot_imps(x)
#   }
#   # if (is.numeric(imp$data[[x]])) {
#     # plot box and whiskers for numeric variables
#     p <- p +
#       ggplot2::geom_boxplot(ggplot2::aes(
#         y = as.factor(.imp),
#         x = .data[[x]],
#         color = datapoint
#       ),
#       width = 0.5,
#       alpha = 0.5,
#       outlier.shape = NA) +
#       scale_y_discrete(limits = rev) +
#       ggplot2::ylab("Imputation (0 = observed data)")
#   # } else {
#   #   # plot faceted barplot for categorical variables
#   #   p <- imp %>% plot_imps(x) +
#   #     ggplot2::geom_bar(ggplot2::aes(y = .data[[x]],
#   #                                    color = datapoint),
#   #                       width = 0.5,
#   #                       fill = "white") +
#   #     ggplot2::facet_wrap(
#   #       ~ .imp,
#   #       scales = "free_x",
#   #       ncol = 2,
#   #       labeller = ggplot2::labeller(.imp = c(
#   #         "Observed data", paste("Imputation", 1:imp$m)
#   #       ) %>% setNames(0:imp$m))
#   #     )
#   # }
#   if(!is.null(cluster)){
#     p <- p + facet_wrap(get(cluster))
#     }
#     # output
#   return(p)
# }

# option 1: bwplot etc with argument
# option 2: ggmice which overwrites the bw functions ect.

# new plotting function for imputed data
plot_imps <- function(imp, type = c("bwplot", "stripplot", "densityplot"), x, y = NULL) {
  # pre-process mids object for plotting
  if(is.null(y)){y <- x}
  comp <- mice::complete(imp, "long")[is.na(imp$data[, x]) | is.na(imp$data[, y]), ] %>% 
    rbind(cbind(.imp = 0, .id = rownames(imp$data), imp$data), .) %>% 
    mutate(.imp = factor(.imp, ordered = TRUE))
  
  # basic plot object
  p <- ggplot2::ggplot(comp) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(
      values = c(mice:::mdc(2), mice:::mdc(1)),
      labels = c("Imputed", "Observed"))
  
  # bwplot or stripplot
  if (type == "stripplot") {
    p <- p + geom_point(
      aes(x = .imp, y = get(x), group = .imp, color = ifelse(.imp < 1, "Observed", "Imputed")),
      position = position_jitter(width = 0.25, height = 0),
      data = p$data) 
  }  
  if (type == "stripplot" | type == "bwplot"){
    p <- p + geom_boxplot(
      aes(x = .imp, y = get(x), group = .imp, color = ifelse(.imp < 1, "Observed", "Imputed")),
      size = 1,
      width = 0.5,
      alpha = 0.5,
      outlier.shape = NA,
      data = p$data) + 
      labs(x = "Imputation",
           y = x,
           color = "",
           fill = "") 
  }
  
  # densityplot
  if (type == "densityplot"){
  p <-  p + geom_density(
    aes(x = get(x), group = .imp, color = ifelse(.imp < 1, "Observed", "Imputed")),
    data = p$data) +
    labs(x = x,
         color = "",
         fill = "")
  }
  
  # # xyplot
  if (type == "xyplot") {
  p <- p + geom_point(
    aes(x = get(x), y = get(y), color = ifelse(.imp < 1, "Observed", "Imputed")),
    data = p$data) +
    labs(x = x,
         y = y,
         color = "")
  }
  
  # output
  return(p)
}

# plot iterations
plot_chains <- function(imp){
# call <- match.call()
if (!is.mids(imp)) {
  stop("argument 'imp' must be a 'mids' object", call. = FALSE)
}
if (is.null(imp$chainMean)) {
  stop("no convergence diagnostics found", call. = FALSE)
}

# extract chain means and chain variances
mn <- imp$chainMean
sm <- sqrt(imp$chainVar)

# select subset of nonmissing entries
obs <- apply(!(is.nan(mn) | is.na(mn)), 1, all)
varlist <- names(obs)[obs]
p <- length(varlist)
m <- imp$m
it <- imp$iteration
dat <- data.frame(
  .ms = rep(c("mean", "sd"), each = m * it * p),
  vrb = rep(varlist, each = m * it, times = 2),
  val = c(matrix(aperm(mn[varlist, , , drop = FALSE], c(2, 3, 1)), nrow = m * it * p),
          matrix(aperm(sm[varlist, , , drop = FALSE], c(2, 3, 1)), nrow = m * it * p))
) %>% cbind(expand.grid(.it = seq_len(it), .m = seq_len(m)), .)

# ## Dummy to trick R CMD check
# .m <- NULL
# rm(.m)

ggplot(dat, aes(x = .it, y = val, color = as.factor(.m))) +
  geom_line() +
  facet_wrap(vrb~.ms, scales = "free", ncol = 2, strip.position = "left") +
  labs(x = "Iteration",
       y = "",
       color = "Imputation") +
  theme_classic() +
  theme(strip.placement = "outside")
}
