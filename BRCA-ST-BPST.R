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
path.data <- "/Users/sy5jx/Dropbox/Research/Papers/Collaboration/SingleCellData/Code/numeric_results/application/BRCA-ST/data/"
Y = readRDS(paste0(path.data, "Y.rds"))
S = readRDS(paste0(path.data, "S.rds"))
X = readRDS(paste0(path.data, "X.rds"))

# plot(S)
# # identify(S, order = TRUE)
# idx <- c(9,11,23,38,84,177,331,416,431,549,566,572,574)
# idx.order <- c(6,5,4,7,8,9,10,11,3,12,2,13,1)
# idx[idx.order] <- idx
# points(S[idx, ], type = "l")
# 
# # boundary and triangulation for bivariate spline over triangulation
# boundary <- S[idx, ]
# boundary[7, ] <- c(31, 12)
# boundary[2, ] <- c(3, 34)
# boundary <- boundary[-c(1, 8, 10), ]
# boundary[boundary[, 1] == 3, 1] <- 2.5
# boundary[boundary[, 1] == 31, 1] <- 31.5
# boundary[boundary[, 2] == 10, 2] <- 9.5
# boundary[boundary[, 2] == 34, 2] <- 34.5
# boundary[7, ] <- c(31.5, 20) 
# boundary[8, ] <- c(29.5,27.5) 
# write.csv(boundary, file = paste0(path.data, "boundary-BRCA-ST.csv"), row.names = FALSE)
boundary <- read.csv(paste0(path.data, "boundary-BRCA-ST.csv"))
boundary <- as.matrix(boundary)
Tr.cell <- TriMesh(boundary, n = 2.5)
points(S[, 1], S[, 2], pch = ".")
points(boundary)
V = as.matrix(Tr.cell$V)
Tr = as.matrix(Tr.cell$Tr)

# spatial pattern testing ----
results <- test.spVC(Y = Y, X, S, V, Tr, para.cores = 1,
                     subset = 1:1000)
# save(results, file = "BRCA-ST-results.RData")

# plot the estimated gamma functions ----
# extract coefficients for one specific gene
results[[2]]$`S100A9`$coeff.beta
# extract p-values
# p values for constant effect starts with "beta_"
# p values for varying effect starts with "gamma_"
results[[2]]$`S100A9`$p.value

# scatter plots based on the original data ----
coef.matrix <- results[[2]]$`S100A9`$coeff.gamma
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
colnames(boundary) <- c("V1", "V2")
boundary <- data.frame(boundary)

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

# summary ----
length(results[[1]])
length(results[[2]])
length(results[[3]])
