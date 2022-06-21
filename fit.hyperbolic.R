###This is a modification from slope.bipartite taken from Bastazini et al. 2019. Environmental Conservation 46.1 (2019): 52-58.

fit.hyperbolica=function (object, plot.it = TRUE, ...) 
{
  if (class(object) != "bipartite") 
    stop("This function cannot be meaningfully applied to objects of this class.")
  N <- colSums(object)
  if (all(object[-nrow(object), 2] == 1)) 
    y <- -object[, 3]
  else y <- -object[, 2]
  y <- (sum(y) - cumsum(y))/sum(y)
  x <- (object[, "no"]/max(object[, "no"]))
  fit <- try(nls(y ~ 1 - x^a, start = list(a = 1)))
  if (class(fit) == "try-error") 
    fit <- nls((y + rnorm(length(y), s = 0.01)) ~ 1 - x^a, 
               start = list(a = 1))
  if (plot.it) {
    par(mar = c(5, 5, 1, 1))
    plot(x, y, xlab = "Fração de manchas eliminadas", 
         ylab = "Fração de espécies de aves sobreviventes", 
         axes = TRUE, type = "n", cex.lab = 1)
    
    
    points(x, y, pch=16, ...)
    lines(seq(0, 1, 0.1), predict(fit, newdata = data.frame(x = seq(0, 
                                                                    1, 0.1))), col = "red", lwd = 2)
  }
  return(c(exponent = as.numeric(coef(fit)[1])))
}
