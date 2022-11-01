#' spCV: Generalized Spatially Varying Coefficient Modeling  (GSVCM) for the 
#' detection of spatial expression patterns for large spatial transcriptomic
#' studies.
#'  
#' This function generates the p-value and estimated coefficients of each 
#' component in the GSVCM.
#'
#' @import BPST MGLM parallel
#' @param Y A \code{m} by \code{n} matrix of count of RNA-seq reads. 
#' Each row represents one genes, and each column represents one cell.
#' @param X A matrix \code{n} by \code{p} of explanatory variables.
#' @param S A matrix \code{n} by two of cell locations.
#' @param V The \code{nV} by two matrix of verities of a triangulation,
#' where \code{nV} is the number of vertices. Each row is the coordinates 
#' for a vertex.
#' @param Tr The triangulation matrix of dimension \code{nTr} by three,
#' where \code{nTr} is the number of triangles in the triangulation.
#' @param para.cores Number of cores in parallel computing.
#' @param p.adjust.method The p-value adjustment methods.
#' @param p.adjust.thresh Significant value for adjusted p-value.
#' @param linear.fit indicator of whether to fit generalized linear model.
#' @param filter.min.nonzero filter genes whose number of nonzero counts is larger than \code{filter.min.nonzero}.
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
library(BPST)
library(MGLM)
library(parallel)
source("fit.spVC.R")

test.spVC = function(Y, X, S, V, Tr, para.cores, scaleX = FALSE, 
                     subset = 1:nrow(Y), p.adjust.method ="BH",
                     p.adjust.thresh = 0.05, linear.fit = FALSE,
                     filter.min.nonzero = 100, fix.contant = NULL){
  
  # standardize location points and boundary
  min.x <- min(V[, 1]); max.x <- max(V[, 1])
  min.y <- min(V[, 2]); max.y <- max(V[, 2])
  V[, 1] <- (V[, 1] - min.x)/max.x
  V[, 2] <- (V[, 2] - min.y)/max.y
  S[, 1] <- (S[, 1] - min.x)/max.x
  S[, 2] <- (S[, 2] - min.y)/max.y
  
  ind <- inVT(V, Tr, S[, 1], S[, 2])$ind.inside
  cat("spCV model will use ", length(ind)/nrow(S)*100,
      "% of the original data.\n")
  d = 2
  
  S.est <- S[ind, ]
  if(scaleX == TRUE){
    X.est = scale(X.est)
  }
  X.est <- cbind(1, X[ind, ])
  colnames(X.est)[1] <- "0"
  Y.est <- Y[, ind]
  size.factors <- colSums(Y.est)
  size.factors <- size.factors/median(size.factors)
  
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
  
  formula.ggam <- as.formula(
    paste0("Y ~ 0 + ", paste0(names(dat.fit)[c(1:(p.X+1))], collapse = " + "))
    )

  idx <- names(which(apply(Y[subset, ], 1, 
                           function(x) sum(x != 0) > filter.min.nonzero)))
  cat("Linear Model & Model 1: Conducting tests for", length(idx), " genes.\n")
  print.idx <- 1:length(idx)
  names(print.idx) <- idx
  results.linear <- NULL
  
  if(linear.fit == TRUE) {
    formula.linear <- as.formula(
      paste0("Y ~ 0 + ", paste0(names(dat.fit)[c(1:(p.X))], collapse = " + "))
    )
    
    results.linear <- mclapply(idx, mc.cores = para.cores,
      FUN = function(x){
        cat("Fitting Linear Model for Gene", print.idx[x], "out of",
            length(idx), "genes.\n")
        mfit.iter <- fit.spVC(formula.linear, Y.iter = as.vector(Y.est[x, ]),
        dat.fit = dat.fit, size.factors = size.factors, pen.list = pen.list)
      }
    )
    names(results.linear) <- idx
  }
  
  
  results.constant <- mclapply(idx, mc.cores = para.cores,
    FUN = function(x){
      cat("Fitting Model 1 for Gene", print.idx[x], "out of", length(idx), "genes.\n")
      mfit.iter <- fit.spVC(formula.ggam, Y.iter = as.vector(Y.est[x, ]),
        dat.fit = dat.fit, size.factors = size.factors, 
        pen.list = pen.list)
      }
    )
                               
  names(results.constant) <- idx
  
  p.value.all <- lapply(results.constant, "[[", 1)
  p.value.mtx <- do.call("rbind", p.value.all)
    
  p.adj <- apply(p.value.mtx[, -1], 2, p.adjust, method = p.adjust.method)
  idx.X <- which(apply(p.adj[, 1:(p.X - 1)], 1, 
                       FUN = function(x) any(x < p.adjust.thresh)))
  idx.S <- which(p.adj[, p.X] < p.adjust.thresh)
  idx.test <- intersect(idx.X, idx.S)
  p.adj.name <- colnames(X.est[, -1])
  cat("Model 2: Conducting tests for", length(idx.test), "genes.\n")
  
  print.idx <- 1:length(idx.test)
  names(print.idx) <- idx[idx.test]
  
  results.varying <- mclapply(idx[idx.test], mc.cores = para.cores,
    FUN = function(x){
      
          # idx.x <- which(idx[idx.test] == x)
          cat("Fitting Model 2 for Gene", print.idx[x], "out of", 
              length(idx.test), "genes.\n")
          varying.set1 <- p.adj.name[p.adj[x, -p.X] < p.adjust.thresh]
          varying.set2 <- setdiff(colnames(X.est)[-1], fix.contant)
          
          v.set <- paste0("gamma_", c("0", intersect(varying.set1,
                                                     varying.set2)))
          c.set <- paste0("beta_",  c("0", p.adj.name))
          formula.iter <- as.formula(
            paste0("Y ~ 0 + ", paste0(c.set, collapse = " + "),
            " + ", paste0(v.set, collapse = " + "))
            )
          mfit.iter <- fit.spVC(formula.iter, Y.iter = Y.est[x, ],
            dat.fit = dat.fit, size.factors = size.factors, 
            pen.list = pen.list)
    }
  )
  
  names(results.varying) <- idx[idx.test]
  
  if(linear.fit == FALSE) {
    results <- list(results.linear = results.linear, results.constant = results.constant,
                results.varying = results.varying,
                BQ2.center.est = colMeans(BQ2))
  }
  
  if(linear.fit == FALSE) {
    results <- list(results.constant = results.constant,
       results.varying = results.varying,
       BQ2.center.est = colMeans(BQ2))
  }
  
  return(results)
}