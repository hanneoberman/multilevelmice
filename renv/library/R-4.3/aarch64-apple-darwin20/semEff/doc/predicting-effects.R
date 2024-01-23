## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message = FALSE----------------------------------------------------------
# install.packages(c("semEff", "ggplot2"))
library(semEff)
library(ggplot2)

## -----------------------------------------------------------------------------
knitr::kable(head(shipley))

## ----eval = FALSE-------------------------------------------------------------
#  shipley.sem <- list(
#    DD = lme4::lmer(DD ~ lat + (1 | site) + (1 | tree), data = shipley),
#    Date = lme4::lmer(Date ~ DD + (1 | site) + (1 | tree), data = shipley),
#    Growth = lme4::lmer(Growth ~ Date + (1 | site) + (1 | tree), data = shipley),
#    Live = lme4::glmer(Live ~ Growth + (1 | site) + (1 | tree), data = shipley,
#                       family = binomial)
#  )

## ----eval = FALSE-------------------------------------------------------------
#  shipley.sem.boot <- bootEff(shipley.sem, R = 1000, seed = 13, ran.eff = "site")

## -----------------------------------------------------------------------------
(shipley.sem.eff <- semEff(shipley.sem.boot))

## -----------------------------------------------------------------------------
summary(shipley.sem.eff, "Growth")

## -----------------------------------------------------------------------------
tot <- getTotEff(shipley.sem.eff, "Growth")
tot.b <- getTotEff(shipley.sem.eff, "Growth", type = "boot")

## -----------------------------------------------------------------------------
dat <- na.omit(shipley)
mod <- shipley.sem$Growth
fit <- sapply(c("Date", "DD"), function(i) {
  x <- data.frame(seq(min(dat[i]), max(dat[i]), length = 100)); names(x) <- i
  f <- predEff(mod, newdata = x, effects = tot[i], eff.boot = tot.b)
  c(x, f)
}, simplify = FALSE)

## -----------------------------------------------------------------------------
plotFit <- function(x, y, fit, x.lab = NULL, y.lab = NULL) {
  d <- fit[[1]]
  f <- fit$fit
  ci.l <- fit$ci.lower
  ci.u <- fit$ci.upper
  ggplot () + 
    geom_point(aes(x, y)) +
    geom_ribbon(aes(d, ymin = ci.l, ymax = ci.u, alpha = "0.15"), fill = "blue") +
    geom_line(aes(d, f), color = "blue", size = 1) +
    xlab(x.lab) + 
    ylab(y.lab) +
    theme_bw() + 
    theme(legend.position = "none")
}

## ----warning = FALSE----------------------------------------------------------
plotFit(x = dat$Date, y = dat$Growth, fit = fit$Date, 
        x.lab = "Julian Date of Bud Burst (direct effect)", y.lab = "Stem Growth (mm)")

## ----warning = FALSE----------------------------------------------------------
plotFit(x = dat$DD, y = dat$Growth, fit = fit$DD, 
        x.lab = "Degree Days to Bud Burst (indirect effect)", y.lab = "Stem Growth (mm)")

## -----------------------------------------------------------------------------
c(r2_marg = R2(mod, re.form = NA)[[1]],
  r2_cond = R2(mod)[[1]])


## -----------------------------------------------------------------------------
c(r2_site = R2(mod, re.form = ~ 1 | site)[[1]],
  r2_tree = R2(mod, re.form = ~ 1 | tree)[[1]])

