varying.test <- function(start.all, stop.all, V.all, coeff.all,
                         edf.all, Xt, rdf){
  p.value.all <- c()
  for(iter in 1:length(start.all)) {
    start <- start.all[iter]
    stop <- stop.all[iter]
    V.i <- V.all[start:stop, start:stop, drop = FALSE]
    p.i <- coeff.all[start:stop]
    edf1i <- edfi <- sum(edf.all[start:stop])
    p.value <- testStat(p = p.i, X = Xt, V = V.i, rank = min(ncol(Xt), edf1i),
                               type = 0, res.df = rdf)$pval
    p.value <- max(2e-17, p.value)
    p.value.all <- c(p.value.all, p.value)
  }
  p.value.all
}
