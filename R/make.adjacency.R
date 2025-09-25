make.adjacency <- function(data,
                            method = "vw",
                            isLocWeighted = FALSE,
                            isModWeighted = FALSE,
                            isSparse = FALSE,
                            ModMethod = NULL,
                            scale = FALSE,
                            sig = 1,
                            radius = NULL,
                            cv.method = "cv.ls"){

  if(!(is.matrix(data) || (is.data.frame(data) && all(sapply(data, is.numeric)))))
    stop("Argument 'data' must be a numeric dataframe or matrix")
  if(!(method %in% c("vw","eu")))
    stop("Argument 'method' must be either 'vw' or 'eu'")
  for(nm in c("isLocWeighted","isModWeighted","isSparse","scale")){
    v <- get(nm, inherits = FALSE)
    if(!(is.logical(v) && length(v) == 1L)) stop(sprintf("Argument '%s' must be a single logical", nm))
  }
  if(isTRUE(isModWeighted)){
    if(is.null(ModMethod) || !(ModMethod %in% c("snn","sim","both")))
      stop("Argument 'ModMethod' must be one of 'snn', 'sim', or 'both' when 'isModWeighted' is TRUE")
  }
  if(!(is.numeric(sig) && length(sig) == 1L && is.finite(sig) && sig > 0))
    stop("Argument 'sig' must be a single positive numeric")
  if(!(cv.method %in% c("cv.ml","cv.ls")))
    stop("Argument 'cv.method' must be either 'cv.ml' or 'cv.ls'")
  if(any(!is.finite(data))) stop("'data' has non-finite values")
  if(scale) data <- scale(data)
  data <- data.matrix(data)
  if(any(!is.finite(data))) stop("scaled 'data' has non-finite values")
  n <- nrow(data)

  if(method == "vw") {
    h <- np::npudensbw(data, bwmethod = cv.method, nmulti = 3)$bw
    w <- 1/h^2
    X.sc <- sweep(data, 2, sqrt(w), "*")
    D <- as.matrix(stats::dist(X.sc))
  } else {
    D <- as.matrix(stats::dist(data))
  }
  if(any(!is.finite(D))) stop("Distance matrix contains non-finite values")

  if(is.null(radius)) r <- find.radius(D) else r <- as.integer(radius)
  if(!is.finite(r)) stop("'radius' must be finite")
  if(r < 1L) stop("'radius' must be at least 1")
  if(r >= n) stop("'radius' must be max(nrow(data)) - 1")

  if(isLocWeighted) {
    cs <- compute.sigma(D, r)
    sigma <- cs$sigma
    if(length(sigma) != n || any(!is.finite(sigma)) || any(sigma <= 0))
      stop("'sigma' from compute.sigma must be length n, positive, and finite")
    S <- exp(-(D^2)/outer(sigma, sigma))
    diag(S) <- 1
    if(isModWeighted){
      if(ModMethod == "snn") {
        SNN <- compute.SNN(S, r)
        M <- if(!isSparse) 0.5*(1+SNN) else SNN
      } else if(ModMethod == "sim") {
        rho <- numeric(n)
        for(i in 1:n) {
          ss <- sort(S[i,], decreasing = TRUE)
          rho[i] <- ss[r + 1L]
        }
        SIM <- sqrt(outer(rho, rho))
        diag(SIM) <- 1
        M <- if(!isSparse) 0.5*(1+SIM) else SIM
      } else if(ModMethod == "both"){
        rho <- numeric(n)
        for(i in 1:n) {
          ss <- sort(S[i,], decreasing = TRUE)
          rho[i] <- ss[r + 1L]
        }
        SIM <- sqrt(outer(rho, rho))
        SNN <- compute.SNN(S, r)
        diag(SIM) <- 1
        M <- if(!isSparse) 0.25*(1+SIM)*(1+SNN) else SIM*SNN
      }
      return(S*M)
    } else {
      return(S)
    }
  } else {
    S <- exp(-(D^2)/(sig^2))
    diag(S) <- 1
    if(isModWeighted){
      if(ModMethod == "snn") {
        SNN <- compute.SNN(S, r)
        M <- if(!isSparse) 0.5*(1+SNN) else SNN
      } else if(ModMethod == "sim") {
        rho <- numeric(n)
        for(i in 1:n) {
          ss <- sort(S[i,], decreasing = TRUE)
          rho[i] <- ss[r + 1L]
        }
        SIM <- sqrt(outer(rho, rho))
        diag(SIM) <- 1
        M <- if(!isSparse) 0.5*(1+SIM) else SIM
      } else if(ModMethod == "both"){
        rho <- numeric(n)
        for(i in 1:n) {
          ss <- sort(S[i,], decreasing = TRUE)
          rho[i] <- ss[r + 1L]
        }
        SIM <- sqrt(outer(rho, rho))
        SNN <- compute.SNN(S, r)
        diag(SIM) <- 1
        M <- if(!isSparse) 0.25*(1+SIM)*(1+SNN) else SIM*SNN
      }
      return(S*M)
    } else {
      return(S)
    }
  }
}
