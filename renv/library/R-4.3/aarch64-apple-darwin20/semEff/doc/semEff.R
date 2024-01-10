## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message = FALSE----------------------------------------------------------
# install.packages(c("semEff", "piecewiseSEM"))
library(semEff)
data(keeley, package = "piecewiseSEM")

## -----------------------------------------------------------------------------
knitr::kable(head(keeley))

## -----------------------------------------------------------------------------
keeley.sem <- list(
  lm(age ~ distance, data = keeley),
  lm(hetero ~ distance, data = keeley),
  lm(abiotic ~ distance, data = keeley),
  lm(firesev ~ age, data = keeley),
  lm(cover ~ firesev, data = keeley),
  lm(rich ~ distance + hetero + abiotic + cover, data = keeley)
)

## ----message = FALSE, warning = FALSE-----------------------------------------
piecewiseSEM:::plot.psem(
  piecewiseSEM::as.psem(keeley.sem),
  node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "grey"),
  layout = "tree"
)

## -----------------------------------------------------------------------------
system.time(
  keeley.sem.boot <- bootEff(keeley.sem, R = 100, seed = 13, parallel = "no")
)

## -----------------------------------------------------------------------------
(keeley.sem.eff <- semEff(keeley.sem.boot))

## -----------------------------------------------------------------------------
summary(keeley.sem.eff, response = "rich")

## -----------------------------------------------------------------------------
summary(
  semEff(keeley.sem.boot, predictor = "distance"),
  response = "rich"
)

## -----------------------------------------------------------------------------
RVIF(keeley.sem[[6]])

