## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##",
  fig.width = 15, fig.align = 'center', out.width = '100%'
)

## ----args-pglmm---------------------------------------------------------------
args(phyr::pglmm_plot_re)

## -----------------------------------------------------------------------------
library(ape)
library(phyr)
suppressPackageStartupMessages(library(dplyr))

set.seed(12345)
nspp <- 7
nsite <- 5
# Simulate a phylogeny that has a lot of phylogenetic signal (power = 1.3)
phy <- compute.brlen(rtree(n = nspp), method = "Grafen", power = 1.3)
# Simulate species means
sd.sp <- 1
mean.sp <- rTraitCont(phy, model = "BM", sigma = sd.sp^2)
Y.sp <- rep(mean.sp, times = nsite)
# Phylogenetically correlated response of species to env
sd.trait <- 1
trait <- rTraitCont(phy, model = "BM", sigma = sd.trait)
trait <- rep(trait, times = nsite)
# Simulate site means
sd.site <- 1
mean.site <- rnorm(nsite, sd = sd.site)
Y.site <- rep(mean.site, each = nspp)
# Site-specific environmental variation
sd.env <- 1
env <- rnorm(nsite, sd = sd.env)
# Generate covariance matrix for phylogenetic attraction
sd.attract <- 1
Vphy <- vcv(phy)
Vphy <- Vphy / (det(Vphy) ^ (1 / nspp))
V.attract <- kronecker(diag(nrow = nsite, ncol = nsite), Vphy)
Y.attract <- array(t(mvtnorm::rmvnorm(n = 1, sigma = sd.attract ^ 2 * V.attract)))
# Residual errors
sd.e <- 1
Y.e <- rnorm(nspp * nsite, sd = sd.e)
# Construct the dataset
d <- data.frame(sp = rep(phy$tip.label, times = nsite), 
                site = rep(1:nsite, each = nspp),
                env = rep(env, each = nspp))
# Simulate abundance data
d$Y <- Y.sp + Y.attract + trait * d$env + Y.e
head(d)

# fit a model
mod_1 = pglmm(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
              data = d, cov_ranef = list(sp = phy))
summary(mod_1)

## ----fig.asp=0.6--------------------------------------------------------------
# plot var-cov matrices of random terms
mod1re = pglmm_plot_re(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
                       data = d, cov_ranef = list(sp = phy), show.image = TRUE, 
                       show.sim.image = FALSE)

## ----fig.asp=0.6--------------------------------------------------------------
# all use color with useAbs = FALSE
pglmm_plot_re(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
              data = d, cov_ranef = list(sp = phy), show.image = TRUE, 
              show.sim.image = FALSE, useAbs = FALSE)

## ----fig.asp=0.6--------------------------------------------------------------
# suppress key with colorkey = FALSE
pglmm_plot_re(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
              data = d, cov_ranef = list(sp = phy), show.image = TRUE, 
              show.sim.image = FALSE, useAbs = FALSE, colorkey = FALSE)

## ----fig.asp=0.6--------------------------------------------------------------
# suppress colorkey, let the function decide whether use color or not
pglmm_plot_re(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
              data = d, cov_ranef = list(sp = phy), show.image = TRUE, 
              show.sim.image = FALSE, colorkey = FALSE)

## ----fig.asp=0.6--------------------------------------------------------------
# all black and white
pglmm_plot_re(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
              data = d, cov_ranef = list(sp = phy), show.image = TRUE, 
              show.sim.image = FALSE, useAbs = TRUE)

## -----------------------------------------------------------------------------
names(mod1re)

## -----------------------------------------------------------------------------
names(mod1re$vcv)

## ----fig.height=6, fig.width=6, out.width='50%'-------------------------------
names(mod1re$plt_re_list)
mod1re$plt_re_list[[6]]

## ----fig.height=6, fig.width=6, out.width='50%'-------------------------------
Matrix::image(mod1re$vcv[[6]], xlab = "", ylab = "", sub = "", main = "1|sp__@site")

## ----fig.asp=0.4--------------------------------------------------------------
gridExtra::grid.arrange(grobs = mod1re$plt_re_list[c(2, 5, 6)], nrow = 1)

## ----fig.asp=0.6--------------------------------------------------------------
# plot simulated matrices of random terms
mod1sim = pglmm_plot_re(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
                        data = d, cov_ranef = list(sp = phy), show.image = FALSE, 
                        show.sim.image = TRUE)

## ----fig.asp=0.6--------------------------------------------------------------
pglmm_plot_re(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
              data = d, cov_ranef = list(sp = phy), show.image = FALSE, 
              show.sim.image = TRUE, add.tree.sp = FALSE)

## ----fig.asp=0.6--------------------------------------------------------------
pglmm_plot_re(Y ~ 1 + env + (1|sp__) + (1|site) + (env|sp__) + (1|sp__@site),
              data = d, cov_ranef = list(sp = phy), show.image = FALSE,
              show.sim.image = TRUE, add.tree.sp = TRUE,
              colorkey = FALSE, useAbs = FALSE)

## -----------------------------------------------------------------------------
names(mod1sim)

## ----fig.width = 12, fig.asp=0.6----------------------------------------------
gridExtra::grid.arrange(grobs = mod1sim$plt_sim_list[c(2, 6)], nrow = 1)
gridExtra::grid.arrange(grobs = lapply(mod1sim$plt_sim_list[c(2, 6)], 
                                       update, 
                                       par.settings = list(layout.heights = 
                                                             list(key.top = 0.3,
                                                                  main = 5))), 
                        nrow = 1)

## ----fig.asp=0.6--------------------------------------------------------------
pglmm_plot_re(x = mod_1, show.image = FALSE, show.sim.image = TRUE, 
              add.tree.sp = TRUE, colorkey = FALSE, useAbs = FALSE)
communityPGLMM.show.re(x = mod_1, show.image = TRUE, show.sim.image = FALSE, useAbs = TRUE)

