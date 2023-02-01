rm(list = ls())

# required packages & source files ----
# for model estimation
# install_github("FIRST-Data-Lab/BPST")
library(BPST)
library(sp)
# install_github("FIRST-Data-Lab/Triangulation")
library(Triangulation)
library(MGLM)
library(mgcv)
# for plots
library(ggplot2)
library(colorRamps)

source("locTV.R")
source("varying.test.R")
source("test.spVC.R")
source("eval.spVC.R")

# input ----
k <- 5000
# file.name <- paste0("../../Code/numeric_results/simulation/data-simu/11-2022/simulation-data-v4-",
# k, ".RData")
file.name <- "~/Dropbox/Research/Papers/Collaboration/SingleCellData/Code/simulated.RData"
load(file.name)
# load("/Users/sy5jx/Dropbox/Research/Papers/Collaboration/SingleCellData/Code/numeric_results/simulation/data-simu/03-2022/simulation-data-v3-8000.RData")
# rm(case2_data, case3_data, case1_data)
X <- Z

# boundary and triangulation for bivariate spline over triangulation
boundary <- read.csv("cell_boundary_v1.csv")/1000
boundary <- boundary[-c(1, 6, 8, 11, 12, 14, 17), ]
Tr.cell <- TriMesh(boundary, n = 2)
points(S[, 1], S[, 2], pch = ".")
V = as.matrix(Tr.cell$V)
Tr = as.matrix(Tr.cell$Tr)

# model fitting ----
res.model1 <- test.spVC(Y = sample[1:5000, ], X = NULL, S, V, Tr, 
                     para.cores = 1, 
                     filter.min.nonzero = -1)
res.model2 <- test.spVC(Y = sample[5001:5100, ], X = NULL, S, V, Tr, 
                        para.cores = 1, 
                        filter.min.nonzero = -1)
res.model3 <- test.spVC(Y = sample[10001:10100, ], X = NULL, S, V, Tr, 
                        para.cores = 1, 
                        filter.min.nonzero = -1)
res.model4 <- test.spVC(Y = sample[15001:15100, ], X = NULL, S, V, Tr, 
                        para.cores = 1, 
                        filter.min.nonzero = -1)

save(res.model1, res.model2, res.model3, res.model4, file = paste0("test.RData"))



