BPST.components <- function(V, Tr, d, r, S.est, X.est) {
  basis.cell <- basis(V = V, Tr = Tr, d = d, r = 1, 
                      Z = as.matrix(S.est))
  
  B <- basis.cell$B; dim(B)
  Q2 <- basis.cell$Q2; dim(Q2)
  
  TV.all <- list()
  Mat <- BPST::build(d - 1)
  for (k in 1:nrow(Tr)) {
    TV.all[[k]] <- locTV(V[Tr[k, 1], ], V[Tr[k, 2], ], V[Tr[k, 3], ], Mat, d)
  }
  P <- bdiag(TV.all)
  BQ2 <- as.matrix(B %*% Q2); dim(BQ2)
  BQ2.center <- scale(BQ2, scale = FALSE)
  BQ2.mean <- colMeans(BQ2)
  PQ2 <- as.matrix(crossprod(Q2, P) %*% Q2); dim(PQ2)
  
  BQ2.center.full <- matrix(0, nrow = length(X.est),
                            ncol = ncol(BQ2.center))
  BQ2.center.full[basis.cell$Ind.inside, ] <- BQ2.center
  
  dat.fit <- kr(matrix(X.est, ncol = 1), BQ2.center.full[, -1])
  pen.mtx <- PQ2[2:ncol(Q2), 2:ncol(Q2)]
  BQ2.mtx <- BQ2.center.full[, -1]
  cat(dim(dat.fit), dim(pen.mtx), dim(BQ2.mtx), "\n")
  
  list(dat.fit = dat.fit, pen.mtx = pen.mtx,
       BQ2.mtx = BQ2.mtx, BQ2.mean = BQ2.mean)
}