#' spCV: Generalized Spatially Varying Coefficient Modeling  (GSVCM) for the
#' detection of spatial expression patterns for large spatial transcriptomic
#' studies.
#'
#' This function generates datasets for simulation studies.
#'
#' @param gene.size Number of genes to simulate for each simulation setting.
#' @param k Number of cells to be randomly generated from the *coords*.
#' @param coords A matrix of all the cell locations.
#' @param Z A matrix of all the covariates.
#' @return
#' \item{sample}{Randomly simulated gene expression level.}
#' \item{S}{A matrix of the random-sampled locations.}
#' \item{Z}{A matrix of the random-sampled covariates.}
#' @export
data.simulation <- function(gene.size, k, coords, Z) {

  sample.idx <- sample(1:nrow(coords), k)
  Z <- Z[sample.idx, ]
  Z <- cbind(1, Z)
  Z <- as.matrix(Z)
  S <- as.matrix(coords[sample.idx, 1:2])

  # simulation setting 1 ----
  beta0 <- runif(gene.size, -2.5, -2.1)
  mu <- exp(beta0)
  mu.mtx.1 <- matrix(rep(mu, each = k), nrow = k)

  # simulation setting 2 ----
  beta0 <- runif(gene.size, -2.5, -2.1)
  beta.mtx <- matrix(2 * runif(4 * gene.size, 0.5, 0.7), ncol = 4)
  beta.sign <- matrix(sample(c(-1, 1), 4 * gene.size, replace = TRUE), ncol = 4)
  beta.zero <- matrix(rep(0, 4 * gene.size), ncol = 4)
  for(iter in 1:gene.size) {
    beta.zero[iter, sample(x = 1:4, size = 2)] = 1
  }
  beta.mtx <- cbind(beta0, beta.mtx * beta.sign * beta.zero)
  mu.mtx.2 <- exp(tcrossprod(Z, beta.mtx))

  # simulation setting 3 ----
  f0.3 <- sapply(1:gene.size, f0.3, S = S)
  beta0 <- runif(gene.size, -2.5, -2.1)
  beta.mtx <- matrix(2 * runif(4 * gene.size, 0.5, 0.7), ncol = 4)
  beta.sign <- matrix(sample(c(-1, 1), 4 * gene.size, replace = TRUE), ncol = 4)
  beta.zero <- matrix(rep(0, 4 * gene.size), ncol = 4)
  for(iter in 1:gene.size) {
    beta.zero[iter, sample(x = 1:4, size = 2)] = 1
  }
  beta.mtx <- cbind(beta0, beta.mtx * beta.sign * beta.zero)
  mu.mtx.3 <- exp(tcrossprod(Z, beta.mtx) + f0.3)

  # simulation setting 4 ----
  # generate f0
  f0.mtx <- sapply(1:gene.size, f0.4, S = S)
  # generate f1 & f2
  f1.mtx <- sapply(1:gene.size, f1, S = S)
  f2.mtx <- sapply(1:gene.size, f2, S = S)

  beta.mtx <- matrix(2 * runif(4 * gene.size, 0.5, 0.7), ncol = 4)
  beta.sign <- matrix(sample(c(-1, 1), 4 * gene.size, replace = TRUE),
                      ncol = 4)
  beta.zero <- matrix(rep(0, 4 * gene.size), ncol = 4)
  for(iter in 1:gene.size) {
    beta.zero[iter, c(2, 4)] = 1
  }
  beta.mtx <- cbind(beta0, beta.mtx * beta.sign * beta.zero)
  mu.mtx.4 <- exp(tcrossprod(Z, beta.mtx) + f0.mtx + 2 * f1.mtx * Z[, 3] +
                    2 * f2.mtx * Z[, 5])

  # aggregate ----
  mu.mtx.all <- cbind(mu.mtx.1, mu.mtx.2, mu.mtx.3, mu.mtx.4)
  mu.mtx.all <- sweep(mu.mtx.all, 1, rowSums(mu.mtx.all), FUN="/")

  size.factors <- runif(k, 0.6, 1.7)
  sample <- sapply(1:k, function(x) {
    rmultinom(1, size = floor(size.factors[x] * 10000), prob = mu.mtx.all[x, ])
  })
  size.factors <- floor(size.factors * 10000)

  list(sample = sample, S = S, Z = Z)
}
