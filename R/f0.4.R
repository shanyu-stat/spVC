f0.4 <- function(iter, S) {
  set.seed(iter)
  a <- runif(1, -3.5, 4)
  b <- runif(1, -3.5, 4)
  sign <- sample(c(1, -1), 1, replace = TRUE)

  f0.dat <- 2 * exp(- 0.05 * (S[, 1] - a)^2 -
                      0.05 * (S[, 2] - b)^2)
  f0.dat <- f0.dat - mean(f0.dat) # + 0.1 * cos(S[, 1] - 3.5)
  f0.dat * sign
}
