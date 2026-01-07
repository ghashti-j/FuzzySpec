######## the following script allows us to the run the LDAs algorithm
######## Xianchao Zhang, Jingwei Li, Hong Yu,

# Xianchao Zhang, Jingwei Li, Hong Yu,
# Local density adaptive similarity measurement for spectral clustering,
# Pattern Recognition Letters,
# Volume 32, Issue 2, 2011,
# Pages 352-358,
# https://doi.org/10.1016/j.patrec.2010.09.014.

ldasm <- function(X, sigsq = 0.15) {
  X  <- as.matrix(X)
  n  <- nrow(X)
  dv <- dist(X)                    
  DM <- as.matrix(dv)              
  sigsq <- 0.15*(max(dv)-min(dv))
  epsilon <- {
    min.d  <- min(dv)
    max.d  <- max(dv)
    mean.d <- mean(dv)
    nn     <- apply(DM, 1, function(r) sort(r)[2])  
    max.n  <- max(nn)
    mean.n <- mean(nn)
    20*mean.d+54*min.d+13-max.n-6*max.d-65*mean.n
  }
  within <- (DM <= epsilon)*1L
  C <- tcrossprod(within)+1           
  W <- exp(-DM/(2*sigsq^2*C))
  diag(W) <- 1
  return(W)
}