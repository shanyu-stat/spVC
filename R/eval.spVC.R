#' spCV: Generalized Spatially Varying Coefficient Modeling  (GSVCM) for the
#' detection of spatial expression patterns for large spatial transcriptomic
#' studies.
#'
#' This function generates estimated gamma functions at given locations
#'
#' @import BPST
#' @param coef.gamma matrix of spline coefficients; each column represents one
#' coefficients for one spline function.
#' @param coef Estimated coefficient matrix of gamma functions from GSVCMs.
#' @param S A matrix \code{n} by two of locations to be evaluated.
#' @param V The \code{nV} by two matrix of verities of a triangulation,
#' where \code{nV} is the number of vertices. Each row is the coordinates
#' for a vertex.
#' @param Tr The triangulation matrix of dimension \code{nTr} by three,
#' where \code{nTr} is the number of triangles in the triangulation.
#' @param center A vector of estimated mean of spline basis functions.
#' @return
#' \item{gamma.value}{estimated gamma functions at points \code{S}}
#' @export
eval.spVC = function(coef.gamma, S, V, Tr, center){

  basis.new <- basis(V = V, Tr = Tr, d = 2, r = 1,
                      Z = as.matrix(S))
  B.new <- basis.new$B
  Q2 <- basis.new$Q2
  BQ2.new <- B.new %*% Q2

  coef.matrix <- rbind(0, coef.gamma)
  gamma.value <- matrix(NA, nrow = nrow(S), ncol = ncol(coef.matrix))
  gamma.value[basis.new$Ind.inside, ] <- as.matrix(sweep(BQ2.new %*% coef.matrix, 2,
                                               center %*% coef.matrix))
  colnames(gamma.value) <- colnames(coef.gamma)

  gamma.value
}
