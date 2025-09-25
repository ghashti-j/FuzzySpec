find.radius <- function(D) {
  n <- nrow(D)
  r <- 1
  nn <- ceiling(sqrt(n))
  while(TRUE) {
    if(r > nn) nn <- r + 10
    dist.obj <- rNN.dist(D, nn)
    NN.Num <- tabulate(dist.obj[, 1:r], nbins = n)
    NN.Num.0 <- length(NN.Num[NN.Num == 0])

    if(r == 1) NN.Num.0.Upd <- NN.Num.0
    if(r > 1 & NN.Num.0.Upd == NN.Num.0) break

    NN.Num.0.Upd <- length(NN.Num[NN.Num == 0])
    r <- r + 1
  }
  return(r)
}
