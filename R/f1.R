f1 <- function(iter, S){
  set.seed(iter)
  a <- runif(1, 1, 5.5)
  sign <- sample(c(1, -1), 1, replace = TRUE)
  f1.dat <- cos(S[, 2] + a)
  f1.dat <- f1.dat - mean(f1.dat)
  f1.dat * sign
}
