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
# rows are genes, columns are cells
Y = readRDS("counts.rds")
S = readRDS("coords.rds")
weights = readRDS("weight.cell.rds")

type = c("Granule", "Bergmann", "Oligodendrocytes", "Purkinje", "MLI2",
         "Astrocytes")
X = as.matrix(weights)[, type]
X = X/rowSums(as.matrix(weights))

# boundary and triangulation for bivariate spline over triangulation
boundary <- read.csv("cell_boundary_v2.csv")
Tr.cell <- TriMesh(boundary, n = 2.1)
points(S[, 1], S[, 2], pch = ".")
V = as.matrix(Tr.cell$V)
Tr = as.matrix(Tr.cell$Tr)

# spatial pattern testing ----
results <- test.spVC(Y = Y, X, S, V, Tr, para.cores = 2, 
                     subset = 1:2000)
save(results, file = "Cerebellum-BPST-results.RData")

# plot the estimated gamma functions ----
# extract coefficients for one specific gene
coef <- results[[2]]$`6430548M08Rik`$coeff.fit
# extract p-values
# p values for constant effect starts with "c_"
# p values for varying effect starts with "v_"
results[[2]]$`6430548M08Rik`$p.value

# scatter plots based on the original data ----
gamma.points <- eval.spVC(coef, S, V, Tr, results$BQ2.center.est)
df <- data.frame(x = S[, 1], y = S[, 2], value = gamma.points[, 1])

ggplot(data = df, aes_string (x = "x", y = "y", col = "value")) +
  geom_point(size = 0.3) + scico::scale_colour_scico(palette = "lajolla")

# raster plots based on the grid points ----
x.grid <- seq(min(boundary[, 1]), max(boundary[, 1]), length.out = 50)
y.grid <- seq(min(boundary[, 2]), max(boundary[, 2]), length.out = 50)
coord.grid <- expand.grid(x.grid, y.grid)
gamma.grid <- eval.spVC(coef, coord.grid, V, Tr, results$BQ2.center.est)

data <- data.frame(coord.grid, gamma.grid[, 1])
names(data) <- c("x", "y", "value")
names(boundary) = c("V1", "V2")

value.range <- range(data$value, na.rm = TRUE)
coord.ratio <- 1
title <- "gamma_0"
# generate plots
p <- ggplot(data = data, aes(x = x, y = y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = matlab.like(104),
                       limits = value.range, na.value = 'transparent') +
  coord_fixed(ratio = coord.ratio) + labs(title = title) +
  theme_bw()+theme(axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(), legend.title=element_blank(),
                   text=element_text(size = 10),legend.key.size = unit(1.5, "cm"),
                   legend.key.width = unit(0.6,"cm"),
                   legend.box.spacing = unit(0.4, "cm"),
                   legend.spacing.x = unit(0.2, 'cm')) +
  geom_polygon(data = boundary, aes(V1, V2), fill=NA, colour="black", size = 1)

p

