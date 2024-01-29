rm(list = ls())

# required packages & source files ----
# for model estimation
# library(spVC)
source("R/eval.spVC.R")
source("R/fit.spVC.R")
source("R/locTV.R")
source("R/test.spVC.R")
source("R/testStat.R")
source("R/varying.test.R")
library(parallel)
library(sp)
# install_github("FIRST-Data-Lab/BPST")
library(BPST)
# install_github("FIRST-Data-Lab/Triangulation")
library(Triangulation)
library(MGLM)
# for plots
library(ggplot2)
library(colorRamps)

# input ----
dir = "/Users/sy5jx/Dropbox/Research/Papers/Collaboration/SingleCellData/Code/numeric_results/application/Cortex-10x/"


Y = readRDS(paste0(dir, "data/Y.rds"))
S = readRDS(paste0(dir, "data/S.rds"))
X = readRDS(paste0(dir, "data/X.rds"))
sub.ind <- sample(1:ncol(Y), 500)
Y <- Y[, sub.ind]
S <- S[sub.ind, ]
X <- X[sub.ind, ]

boundary = read.csv(paste0(dir, "data/boundary-Cortex-10x.csv"))

Tr.cell <- TriMesh(boundary, n = 1.8)
points(S[, 1], S[, 2], pch = ".")
V <- as.matrix(Tr.cell$V)
Tr <- as.matrix(Tr.cell$Tr)

# model fitting ----
pvalues <- c()
for (iter in 1:500) {
  # shuffle the locations
  idx.random <- sample(sample(1:nrow(S)))

  res.model1 <- test.spVC(Y = Y,
                          X = X[idx.random, ],
                          S = S[idx.random, ], V, Tr,
                          para.cores = 1,
                          subset = rownames(Y)[1:29],
                          filter.min.nonzero = 100,
                          fix.constant = colnames(X),
                          linear.fit = TRUE)

  pvalues.iter <- res.model1$results.constant[[2]]$p.value[8]
  max.gamma <- max(abs(res.model1$results.constant[[2]]$coeff.gamma))
  quantile.gamma <- quantile(abs(res.model1$results.constant[[2]]$coeff.gamma), 0.75)
  cat(pvalues.iter, max.gamma, quantile.gamma, "\n")
  pvalues <- c( pvalues, pvalues.iter)
}

pvalues.orginal <- c()
for(iter in 1:length(res.model1$results.constant)) {
  pvalues.orginal <- c(pvalues.orginal,
                       res.model1$results.constant[[iter]]$p.value[8])
}
hist(pvalues[pvalues > 0.001])
hist(pvalues)
res.model2 <- test.spVC(Y = Y,
                        X = X, S[idx.random, ], V, Tr,
                        para.cores = 2,
                        filter.min.nonzero = 100,
                        fix.constant = colnames(X))



