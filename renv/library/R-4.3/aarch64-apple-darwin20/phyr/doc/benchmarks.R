## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, eval = FALSE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
#  library(phyr)
#  
#  # simulate data
#  nspp = 500
#  nsite = 100
#  tree_sim = ape::rtree(n = nspp)
#  comm_sim = matrix(rbinom(nspp * nsite, size = 1, prob = 0.6),
#                    nrow = nsite, ncol = nspp)
#  row.names(comm_sim) = paste0("site_", 1:nsite)
#  colnames(comm_sim) = paste0("t", 1:nspp)
#  comm_sim = comm_sim[, tree_sim$tip.label]
#  # about 40 times faster
#  rbenchmark::benchmark(
#    "picante" = {picante::psv(comm_sim, tree_sim)},
#    "phyr R" = {phyr::psv(comm_sim, tree_sim, cpp = FALSE)},
#    "phyr c++" = {phyr::psv(comm_sim, tree_sim, cpp = TRUE)},
#    replications = 10,
#    columns = c("test", "replications", "elapsed",
#                "relative", "user.self", "sys.self"))
#  #>       test replications elapsed relative user.self sys.self
#  #> 3 phyr c++           10   0.339    1.000     0.298    0.030
#  #> 2   phyr R           10   3.287    9.696     2.907    0.303
#  #> 1  picante           10  16.265   47.979    14.824    0.795

## -----------------------------------------------------------------------------
#  comm_sim = matrix(rpois(nspp * nsite, 3), nrow = nsite, ncol = nspp)
#  row.names(comm_sim) = paste0("site_", 1:nsite)
#  colnames(comm_sim) = paste0("t", 1:nspp)
#  comm_sim = comm_sim[, tree_sim$tip.label]
#  # about 2-3 times faster
#  rbenchmark::benchmark(
#    "picante" = {picante::pse(comm_sim, tree_sim)},
#    "phyr R" = {phyr::pse(comm_sim, tree_sim, cpp = FALSE)},
#    "phyr c++" = {phyr::pse(comm_sim, tree_sim, cpp = TRUE)},
#    replications = 20,
#    columns = c("test", "replications", "elapsed",
#                "relative", "user.self", "sys.self"))
#  #>       test replications elapsed relative user.self sys.self
#  #> 3 phyr c++           20   1.456    1.000     1.329    0.105
#  #> 2   phyr R           20   4.233    2.907     3.453    0.555
#  #> 1  picante           20   3.858    2.650     3.319    0.475

## ---- message=FALSE-----------------------------------------------------------
#  # pcd is about 20 times faster
#  rbenchmark::benchmark(
#    "phyr" = {phyr::pcd(comm = comm_a, tree = phylotree, reps = 1000, verbose = FALSE)},
#    "picante" = {picante::pcd(comm = comm_a, tree = phylotree, reps = 1000)},
#    replications = 10,
#    columns = c("test", "replications", "elapsed",
#                "relative", "user.self", "sys.self"))
#  #>      test replications elapsed relative user.self sys.self
#  #> 1    phyr           10   0.214    1.000     0.192    0.012
#  #> 2 picante           10   4.516   21.103     4.043    0.074

## -----------------------------------------------------------------------------
#  library(ape)
#  # Set up parameter values for simulating data
#  n <- 50
#  phy <- rcoal(n, tip.label = 1:n)
#  trt_names <- paste0("par", 1:2)
#  
#  R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2)
#  d <- c(0.3, 0.95)
#  B2 <- 1
#  
#  Se <- c(0.2, 1)
#  M <- matrix(Se, nrow = n, ncol = 2, byrow = TRUE)
#  colnames(M) <- trt_names
#  
#  # Set up needed matrices for the simulations
#  p <- length(d)
#  
#  star <- stree(n)
#  star$edge.length <- array(1, dim = c(n, 1))
#  star$tip.label <- phy$tip.label
#  
#  Vphy <- vcv(phy)
#  Vphy <- Vphy/max(Vphy)
#  Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)
#  
#  tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
#  C <- matrix(0, nrow = p * n, ncol = p * n)
#  for (i in 1:p) for (j in 1:p) {
#    Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
#    C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
#  }
#  MM <- matrix(M^2, ncol = 1)
#  V <- C + diag(as.numeric(MM))
#  
#  iD <- t(chol(V))
#  
#  XX <- iD %*% rnorm(2 * n)
#  X <- matrix(XX, n, p)
#  colnames(X) <- trt_names
#  rownames(X) <- phy$tip.label
#  rownames(M) <- phy$tip.label
#  
#  U <- list(cbind(rnorm(n, mean = 2, sd = 10)))
#  names(U) <- trt_names[2]
#  
#  X[,2] <- X[,2] + B2[1] * U[[1]][,1] - B2[1] * mean(U[[1]][,1])
#  
#  z <- cor_phylo(variates = X,
#                 covariates = U,
#                 meas_errors = M,
#                 phy = phy,
#                 species = phy$tip.label)
#  
#  
#  U2 <- list(NULL, matrix(rnorm(n, mean = 2, sd = 10), nrow = n, ncol = 1))
#  rownames(U2[[2]]) <- phy$tip.label
#  colnames(U2[[2]]) <- "par2"
#  X2 = X
#  X2[,2] <- X2[,2] + B2[1] * U2[[2]][,1] - B2[1] * mean(U2[[2]][,1])
#  
#  z_r <- corphylo(X = X2, SeM = M, U = U2, phy = phy, method = "Nelder-Mead")
#  
#  rbenchmark::benchmark(
#    "cor_phylo" = {cor_phylo(variates = X, covariates = U, meas_errors = M,
#                             phy = phy, species = phy$tip.label)},
#    "corphylo" = {corphylo(X = X2, SeM = M, U = U2, phy = phy, method = "Nelder-Mead")},
#    replications = 5,
#    columns = c("test", "replications", "elapsed",
#                "relative", "user.self", "sys.self")
#  )
#  #>        test replications elapsed relative user.self sys.self
#  #> 1 cor_phylo            5   4.511    1.000     4.329    0.062
#  #> 2  corphylo            5  16.190    3.589    13.863    1.369

