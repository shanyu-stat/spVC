f2 <- function(iter, S){
  set.seed(iter)
  b <- runif(1, 1, 5)
  sign <- sample(c(1, -1), 1, replace = TRUE)
  f2.dat <- cos(S[, 1] + b)
  f2.dat <- f2.dat - mean(f2.dat)
  f2.dat * sign
}
