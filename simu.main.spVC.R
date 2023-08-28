rm(list = ls())

# required packages & source files ----
# for model estimation
library(spVC)
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
k <- 8000
gene.size <- 5000

# generate random dataset ----
data(coords)
coords$x <- coords$x/1000
coords$y <- coords$y/1000
data(weight.cell)
df <- data.frame(x = coords$x, y = coords$y, as.matrix(weight.cell))
Z <- df[, c("Bergmann", "Granule", "MLI2", "Oligodendrocytes")]
Z <- as.matrix(Z)
data.sample <- data.simulation(gene.size, k, coords, Z)

X <- data.sample$Z[, -1]
S <- data.sample$S
sample <- data.sample$sample

# boundary and triangulation for bivariate spline over triangulation
boundary <- data(cell.boundary)
boundary <- cell.boundary/1000
boundary <- boundary[-c(1, 6, 8, 11, 12, 14, 17), ]
Tr.cell <- TriMesh(boundary, n = 2)
points(S[, 1], S[, 2], pch = ".")
V <- as.matrix(Tr.cell$V)
Tr <- as.matrix(Tr.cell$Tr)

# model fitting ----
res.model1 <- test.spVC(Y = sample, X = X, S, V, Tr,
                        para.cores = 1, filter.min.nonzero = -1)




