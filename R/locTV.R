#' Calculate local energy matrix
#'
#' This function computes the local
#' energy matrix with respect to the triangle \code{(V1,V2,V3)}.
#' @param V1,V2,V3
#' three vectors of length two indicating the vertices locations of a triangle.
#' @param Mat inner coefficients matrix. This arguments could be
#' calculated by function \code{build}.
#' @param d degree of polynomials.
#' @return energy matrix with respect to the triangle \code{(V1,V2,V3)}.
#' @export
locTV <- function(V1, V2, V3, Mat, d) {
  Mat <- as.matrix(Mat)
  m <- (d + 1) * (d + 2) / 2
  Id <- diag(m)
  vx <- c(1, 0)
  vy <- c(0, 1)
  lamx <- tcord(V1, V2, V3, vx)
  lamy <- tcord(V1, V2, V3, vy)
  Dx <- dirder(Id, lamx[1], lamx[2], lamx[3])
  Dy <- dirder(Id, lamy[1], lamy[2], lamy[3])

  area <- triarea(V1, V2, V3)
  Mat <- abs(area) * Mat
  K <- (crossprod(Dx, Mat) %*% Dx + crossprod(Dy, Mat) %*% Dy)
  return(K)
}
