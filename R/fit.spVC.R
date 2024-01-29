#' spCV: Generalized Spatially Varying Coefficient Modeling  (GSVCM) for the
#' detection of spatial expression patterns for large spatial transcriptomic
#' studies.
#'
#' This function generates the p-value and estimated coefficients of each
#' component in the GSVCM.
#'
#' @import mgcv stats
#' @param formula.spVC A formula of Generalized Partially Spatially Varying
#' Coefficient Model.
#' @param Y.iter A vector of count of RNA-seq reads with length \code{n}.
#' @param dat.fit A list of data containing explanatory variables,
#' spline basis functions, and the kronecker product of explanatory variables and
#' spline basis functions.
#' @param size.factors A vector of size.factors with length \code{n}.
#' @param pen.list A list of roughness penalty in estimating spatial patterns.
#' @return
#' \item{p.value}{p-values of model components}
#' \item{coeff.beta}{estimate beta}
#' \item{coeff.gamma}{estimate basis coefficients of gamma functions}
#' @export

fit.spVC <- function(formula.spVC, Y.iter, dat.fit, size.factors, pen.list){

  # x ="Calb2"
  # Y.iter <- as.vector(Y.est[x, ])
  # dat.fit$Y = NULL
  dat.fit[length(dat.fit) + 1] <- Y.iter
  names(dat.fit)[length(dat.fit)] = "Y"
  mfit.spVC <- gam(formula.spVC, data = dat.fit,
                  family = quasipoisson(), offset = log(size.factors),
                  paraPen = pen.list)
  # R2 <- max(0, 1 - mfit.spVC$deviance/mfit.spVC$null.deviance)
  Deviance <- mfit.spVC$deviance
  spVC.term <- attr(mfit.spVC$terms, "term.labels")
  idx.c <- which(substr(spVC.term, 1, 5) == "beta_")
  idx.v <- which(substr(spVC.term, 1, 6) == "gamma_")

  p.value.c <- anova(mfit.spVC)$pTerms.pv[idx.c]
  coeff.fit <- mfit.spVC$coefficients
  coeff.beta <- mfit.spVC$coefficients[idx.c]
  p.value.v <- c()
  coeff.gamma <- c()

  if(length(idx.v) != 0) {
    coeff.gamma <- matrix(mfit.spVC$coefficients[-idx.c],
                          ncol = length(idx.v))
    colnames(coeff.gamma) <- spVC.term[idx.v]

    p.X <- length(idx.c)
    V.all <- mfit.spVC$Vp
    edf.all <- mfit.spVC$edf
    Xt <- dat.fit$gamma_0
    rdf <- length(mfit.spVC$y) - sum(mfit.spVC$edf)
    dim.BQ2 <- cumsum(rep(ncol(Xt), length(idx.v)))
    if(length(idx.v) == 1) {
      start.all <- p.X + 1
    } else {
      start.all <- c(p.X + 1, dim.BQ2[-length(idx.v)] + p.X + 1)
    }
    stop.all <- dim.BQ2 + p.X

    p.value.v <- varying.test(start.all, stop.all,
                              V.all, coeff.fit, edf.all, Xt, rdf)
    # p.value.v
    names(p.value.v) <- spVC.term[idx.v]
  }

  p.value <- c(p.value.c, p.value.v)

  list(p.value = p.value, coeff.beta = coeff.beta,
       coeff.gamma = coeff.gamma, Deviance = Deviance)
}
