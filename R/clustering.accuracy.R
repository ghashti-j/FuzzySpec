clustering.accuracy <- function(A,B){
  if(length(A) != length(B)) stop("Both vectors must be of same length")
  if(nrow(table(A,B)) != ncol(table(A,B))) stop("Number of unique values in each vector must match")
  acc <- sum(diag(Thresher::matchLabels(table(A,B))))/length(A)
  return(acc)
}

