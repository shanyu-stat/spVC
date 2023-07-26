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
setwd("~/Dropbox/Research/Papers/Collaboration/SingleCellData/Code/numeric_results/application/Cortex-10x/permutation/")
Y = readRDS("counts.rds")
S = readRDS("S-permuted.rds")
weights = readRDS("weight.cell.rds")

# boundary and triangulation for bivariate spline over triangulation
boundary <- read.csv(paste0(path.data, "boundary-Cortex-10x.csv"))
Tr.cell <- TriMesh(boundary, n = 1.8)
points(S[, 1], S[, 2], pch = ".")
points(boundary)
V = as.matrix(Tr.cell$V)
Tr = as.matrix(Tr.cell$Tr)

# spatial pattern testing ----
results <- test.spVC(Y = Y, X, S, V, Tr, para.cores = 1)

p.value.all <- sapply(1:length(results$results.constant), 
                      function(iter) results$results.constant[[iter]]$p.value[8])
hist(p.value.all)
p.adjust.all <- p.adjust(p.value.all, method = "BH")
mean(p.adjust.all < 0.05)
save(results, file = "Cerebellum-BPST-results.RData")

# plot the estimated gamma functions ----
# extract coefficients for one specific gene
results[[2]]$`Apoe`$coeff.beta
# extract p-values
# p values for constant effect starts with "beta_"
# p values for varying effect starts with "gamma_"
results[[2]]$`Apoe`$p.value

# scatter plots based on the original data ----
coef.matrix <- results[[2]]$`Apoe`$coeff.gamma
gamma.points <- eval.spVC(coef.matrix, S, V, Tr, results$BQ2.center.est)
dim(gamma.points)
colnames(gamma.points)

for(iter in 1:ncol(gamma.points)){
  df <- data.frame(x = S[, 1], y = S[, 2], value = gamma.points[, iter])
  title <- colnames(gamma.points)[iter]
  p <- ggplot(data = df, aes_string (x = "x", y = "y", col = "value")) +
    geom_point(size = 0.3) + 
    # scico::scale_colour_scico(palette = "lajolla") +
    scale_color_gradient2(low = "blue", high = "red", mid = "white",
                          midpoint = 0, na.value = 'transparent') +
    theme_classic() +
    labs(title = title)
  print(p)
}


# raster plots based on the grid points ----
x.grid <- seq(min(boundary[, 1]), max(boundary[, 1]), length.out = 50)
y.grid <- seq(min(boundary[, 2]), max(boundary[, 2]), length.out = 50)
coord.grid <- expand.grid(x.grid, y.grid)
gamma.grid <- eval.spVC(coef.matrix, coord.grid, V, Tr, results$BQ2.center.est)
names(boundary) = c("V1", "V2")

for(iter in 1:ncol(gamma.points)){
  data <- data.frame(coord.grid, gamma.grid[, iter])
  names(data) <- c("x", "y", "value")
  value.range <- range(data$value, na.rm = TRUE)
  title <- colnames(gamma.grid)[iter]
  
  # generate plots
  p.iter <- ggplot(data = data, aes(x = x, y = y)) +
    geom_raster(aes(fill = value)) +
    # scale_fill_gradientn(colours = matlab.like(104),
    #                      limits = value.range, na.value = 'transparent') +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, na.value = 'transparent') +
    coord_fixed(ratio = 1) + labs(title = title) +
    theme_bw()+theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), legend.title=element_blank(),
                     text=element_text(size = 10),legend.key.size = unit(1.5, "cm"),
                     legend.key.width = unit(0.6,"cm"),
                     legend.box.spacing = unit(0.4, "cm"),
                     legend.spacing.x = unit(0.2, 'cm')) +
    geom_polygon(data = boundary, aes(V1, V2), fill=NA, 
                 colour="black", size = 1)
  
  print(p.iter)
}



