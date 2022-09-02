library(ggplot2)
library(reshape2)

load("Cerebellum-BPST-results.RData")
terms.name <- c("Int", "Granule", "Bergmann", "Oligodendrocytes", 
                "Purkinje", "MLI2", "Astrocytes")

result.matrix <- lapply(results$results.interact, function(x){
  p.value.all = matrix(NA, ncol = 14, nrow = 1)
  colnames(p.value.all) = c(paste0("c_", terms.name), paste0("v_", terms.name))
  p.value.all[1, names(x$p.value)] = x$p.value
  p.value.all
})
result.matrix <- do.call("rbind", result.matrix)
rownames(result.matrix) <- names(results$results.interact)

Y = readRDS("counts.rds")[names(results$results.interact), ]
count.all = rowSums(Y)
gene.order = order(count.all, decreasing = TRUE)

result.df <- melt(result.matrix[gene.order[1:50], 8:14])
colnames(result.df) <- c("gene", "cell_type", "Pvalue")
result.df <- result.df[complete.cases(result.df), ]

ggplot(result.df, aes(x = cell_type, y = gene)) + 
  geom_raster(aes(fill = log(Pvalue + 10^(-6)))) + 
  scale_fill_gradient(low="red", high="grey90") +
  labs(y="gene", x="cell type") +
  theme_bw() + 
  theme(axis.text.x=element_text(size=8, angle=20, vjust=0.3),
        axis.text.y=element_text(size=8),
        plot.title=element_text(size=11))

