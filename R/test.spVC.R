#' spVC: Generalized Spatially Varying Coefficient Modeling  (GSVCM) for the
#' detection of spatial expression patterns for large spatial transcriptomic
#' studies.
#'
#' This function generates the p-value and estimated coefficients of each
#' component in the GSVCM.
#'
#' @import BPST Matrix MGLM parallel
#' @param Y A \code{m} by \code{n} matrix of count of RNA-seq reads.
#' Each row represents one genes, and each column represents one cell.
#' @param X A matrix \code{n} by \code{p} of explanatory variables.
#' @param S A matrix \code{n} by two of cell locations.
#' @param V The \code{nV} by two matrix of verities of a triangulation,
#' where \code{nV} is the number of vertices. Each row is the coordinates
#' for a vertex.
#' @param Tr The triangulation matrix of dimension \code{nTr} by three,
#' where \code{nTr} is the number of triangles in the triangulation.
#' @param para.cores number of cores in parallel computing.
#' @param scaleX an indicator of whether to standardize the covariates.
#' @param subset a vector of genes to be tested.
#' @param p.adjust.method The p-value adjustment methods.
#' @param p.adjust.thresh significant value for adjusted p-value.
#' @param linear.fit indicator of whether to fit generalized linear model.
#' @param reduced.only indicator of whether to consider the spatially varying
#' coefficients of the covariates.
#' @param twostep indicator of whether to use the two step model estimation
#'  pipeline or directly fit the full model and evaluate the significance of
#'  each individual components in the spVC model.
#' @param filter.min.nonzero filter genes whose number of nonzero counts is
#' larger than \code{filter.min.nonzero}.
#' @param filter.spot.counts filter spots whose total gene expression counts is
#' larger than \code{filter.spot.counts}.
#' @param fix.constant subset of covariates only considering the constant effect
#' in the model.
#' @param fix.varying subset of covariates considering the spatially varying
#' coefficients.
#' @param size.factors a vector of given size.factor for each cell. If
#' \code{size.factor = NULL}, the size.factor is calculated from data.
#' @return
#' \item{results.linear}{A list of spatial pattern testing results for
#' each gene based on the generalized linear model. We conduct the test for
#' genes with more than one hundred nonzero counts. Each element contains
#' p-values of model components and estimate coefficients.}
#' \item{results.constant}{A list of spatial pattern testing results for
#' each gene based on Model 1. We conduct the test for genes with more than one hundred nonzero
#' counts. Each element contains p-values of model components and estimate
#' coefficients.}
#' \item{results.varying}{A list of spatial pattern testing results for
#' each gene based on Model 2. We conduct the test for genes with significant spatial pattern
#' and explanatory variable. Each element contains p-values of model components
#' and estimate coefficients.}
#' \item{BQ2.center.est}{A vector of estimated mean of spline basis functions.}
#' @export

test.spVC <- function(Y, X = NULL, S, V, Tr, para.cores = 1, scaleX = FALSE,
                      subset = 1:nrow(Y), p.adjust.method ="BH",
                      p.adjust.thresh = 0.05, linear.fit = FALSE,
                      reduced.only = FALSE, twostep = TRUE,
                      filter.min.nonzero = 100, filter.spot.counts = 100,
                      fix.constant = NULL, fix.varying = NULL,
                      size.factors = NULL){
  # data prep ----
  # standardize location points and boundary
  min.x <- min(V[, 1]); max.x <- max(V[, 1])
  min.y <- min(V[, 2]); max.y <- max(V[, 2])
  V[, 1] <- (V[, 1] - min.x)/max.x
  V[, 2] <- (V[, 2] - min.y)/max.y
  S[, 1] <- (S[, 1] - min.x)/max.x
  S[, 2] <- (S[, 2] - min.y)/max.y

  # filter spots with low gene expression
  idx.s <- which(colSums(Y) >= filter.spot.counts)

  S <- S[idx.s, ]
  ind <- inVT(V, Tr, S[, 1], S[, 2])$ind.inside
  cat("spVC model will use ", length(ind)/ncol(Y)*100, "% of the original data.\n")
  S.est <- S[ind, ]

  if(is.null(X)) {
    X.est <- matrix(1, nrow = length(ind), ncol = 1)
  } else {
    X <- matrix(X, nrow = ncol(Y))
    X <- X[idx.s, ]
    X <- matrix(X, nrow = length(idx.s))
    if(is.null(colnames(X))) colnames(X) <- paste0("X", 1:ncol(X))
    X.est <- cbind(1, X[ind, ])
    colnames(X.est) <- c("0", colnames(X))
    if(scaleX == TRUE){
      X.est[, -1] = scale(X.est[, -1])
    }
  }
  colnames(X.est)[1] <- "0"

  Y <- Y[, idx.s]
  Y.est <- Y[, ind]

  if(is.null(size.factors)){
    size.factors <- colSums(Y.est)
    size.factors <- size.factors/median(size.factors)
  } else {
    size.factors <- size.factors[ind]
  }

  d = 2
  basis.cell <- basis(V = V, Tr = Tr, d = d, r = 1, Z = as.matrix(S.est))

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
  PQ2 <- as.matrix(crossprod(Q2, P) %*% Q2); dim(PQ2)
  dat.fit <- as.data.frame(X.est)
  pen.list = list()
  p.X = ncol(X.est)
  for(ii in 1:p.X){
    dat.fit[[p.X + ii]] = kr(matrix(X.est[, ii], ncol = 1), BQ2.center[, -1])
    pen.list[[ii]] = list(PQ2[2:ncol(Q2), 2:ncol(Q2)])
  }
  names(dat.fit)[1:p.X] = paste0("beta_", colnames(X.est))
  names(dat.fit)[p.X + 1:p.X] = paste0("gamma_", colnames(X.est))
  names(pen.list) = paste0("gamma_", colnames(X.est))

  if(is.null(rownames(Y.est))) rownames(Y.est) <- paste0("gene", 1:nrow(Y.est))
  idx <- names(which(apply(Y.est[subset, ], 1,
                           function(x) sum(x != 0) > filter.min.nonzero)))
  cat("Conducting tests for", length(idx), " genes.\n")
  print.idx <- 1:length(idx)
  names(print.idx) <- idx
  results.linear <- NULL

  # fit linear models ----
  if(linear.fit == TRUE) {
    if(is.null(X)) warning("The covariate matrix is NULL. The model will fit constants models.")
    formula.linear <- as.formula(
      paste0("Y ~ 0 + ", paste0(names(dat.fit)[c(1:(p.X))], collapse = " + "))
    )

    results.linear <- mclapply(idx, mc.cores = para.cores,
                               FUN = function(x){
                                 if(print.idx[x] %% 500 == 0) {
                                   cat("Fitting Linear Model for Gene", print.idx[x],
                                       "out of", length(idx), "genes.\n")
                                 }
                                 mfit.iter <- fit.spVC(formula.linear, Y.iter = as.vector(Y.est[x, ]),
                                                       dat.fit = dat.fit, size.factors = size.factors, pen.list = pen.list)
                               }
    )
    names(results.linear) <- idx
  }

  # full spVC and evaluate the significance of each individual component ----
  if(twostep == FALSE){
    formula.full <- as.formula(paste0("Y ~ 0 + ",
                                      paste0(names(dat.fit), collapse = " + ")))

    results.full <- mclapply(idx, mc.cores = para.cores,
                                 FUN = function(x){
                                   if(print.idx[x] %% 500 == 0) {
                                     cat("Fitting Full Model for Gene", print.idx[x],
                                         "out of", length(idx), "genes.\n")
                                   }
                                   mfit.iter <- fit.spVC(formula.full, Y.iter = as.vector(Y.est[x, ]),
                                                         dat.fit = dat.fit, size.factors = size.factors,
                                                         pen.list = pen.list)
                                 }
    )

    names(results.full) <- idx
  }

  # fit model with spatial effect ----
  if(twostep == TRUE){
    formula.ggam <- as.formula(
      paste0("Y ~ 0 + ", paste0(names(dat.fit)[c(1:(p.X+1))], collapse = " + "))
    )

    results.constant <- mclapply(idx, mc.cores = para.cores,
                                 FUN = function(x){
                                   if(print.idx[x] %% 500 == 0) {
                                     cat("Fitting Model 1 for Gene", print.idx[x],
                                         "out of", length(idx), "genes.\n")
                                   }
                                   mfit.iter <- fit.spVC(formula.ggam, Y.iter = as.vector(Y.est[x, ]),
                                                         dat.fit = dat.fit, size.factors = size.factors,
                                                         pen.list = pen.list)
                                 }
    )

    names(results.constant) <- idx

    # fit spatially varying coefficient models ----
    if(reduced.only == FALSE) {
      if(!is.null(X) & length(fix.constant) != length(colnames(X))) {
        p.value.all <- lapply(results.constant, "[[", 1)
        p.value.mtx <- do.call("rbind", p.value.all)
        p.adj.name <- colnames(X.est)[-1]
        p.adj <- apply(p.value.mtx[, -1], 2, p.adjust, method = p.adjust.method)
        if(is.null(fix.varying)){
          idx.X <- which(
            apply(
              as.matrix(p.adj[, 1:(p.X - 1)]), 1,
              FUN = function(x){
                varying.set1 <- p.adj.name[x < p.adjust.thresh]
                varying.set2 <- union(setdiff(colnames(X.est)[-1], fix.constant),
                                      fix.varying)
                length(intersect(varying.set1, varying.set2)) > 0
              }
            )
          )
          idx.S <- which(p.adj[, p.X] < p.adjust.thresh)
          idx.test <- intersect(idx.X, idx.S)
        } else {
          idx.test <- 1:length(idx)
        }

        cat("Model 2: Conducting tests for", length(idx.test), "genes.\n")

        print.idx <- 1:length(idx.test)
        names(print.idx) <- idx[idx.test]

        results.varying <-
          mclapply(idx[idx.test], mc.cores = para.cores,
                   FUN = function(x){
                     if(print.idx[x] %% 500 == 0) {
                       cat("Fitting Model 2 for Gene", print.idx[x],
                           "out of", length(idx), "genes.\n")
                     }
                     varying.set1 <- p.adj.name[p.adj[x, -p.X] < p.adjust.thresh]
                     varying.set2 <- union(setdiff(colnames(X.est)[-1], fix.constant),
                                           fix.varying)
                     v.set <- paste0("gamma_",
                                     c("0", intersect(varying.set1, varying.set2)))
                     c.set <- paste0("beta_",  c("0", p.adj.name))
                     formula.iter <- as.formula(
                       paste0("Y ~ 0 + ", paste0(c.set, collapse = " + "),
                              " + ", paste0(v.set, collapse = " + ")))
                     mfit.iter <- fit.spVC(formula.iter, Y.iter = Y.est[x, ],
                                           dat.fit = dat.fit,
                                           size.factors = size.factors,
                                           pen.list = pen.list)
                   }
          )
        names(results.varying) <- idx[idx.test]
      }
    }
  }

  # prepare output ----
  if (is.null(X) | length(fix.constant) == length(colnames(X)) | reduced.only == TRUE) {
    results <- list(results.constant = results.constant,
                    BQ2.center.est = colMeans(BQ2))
  }

  if (twostep == FALSE) {
    results <- list(results.full = results.full,
                    BQ2.center.est = colMeans(BQ2))
  }
  if (twostep == TRUE) {
    if(linear.fit == TRUE) {
      results <- list(results.linear = results.linear,
                      results.constant = results.constant,
                      results.varying = results.varying,
                      BQ2.center.est = colMeans(BQ2))
    }

    if(linear.fit == FALSE) {
      results <- list(results.constant = results.constant,
                      results.varying = results.varying,
                      BQ2.center.est = colMeans(BQ2))
    }
  }

  return(results)
}
