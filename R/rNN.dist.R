rNN.dist <- function(D, r) {
  n <- nrow(D)
  id <- t(apply(D, 1, function(d) order(d)[-1][1:r]))
  return(id)
}
