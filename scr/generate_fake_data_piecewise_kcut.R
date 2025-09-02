
# Individual: Segmented index sampling event duration (without PH)
rpexp_1 <- function(lambda, cuts){
  c0 <- c(0, head(cuts, -1)); c1 <- cuts
  t <- 0
  for (j in seq_along(cuts)) {
    rate <- lambda[j]
    L    <- c1[j] - c0[j]
    if (!is.finite(rate) || rate <= 0) { t <- t + L; next }
    w <- rexp(1, rate)
    if (w < L) return(t + w)   
    t <- t + L                 
  }
  t
}

#Batch generation: Observations within window A (administrative deletion)
simulate_pwe <- function(n, lambda, cuts, A, seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  T <- replicate(n, rpexp_1(lambda, cuts))
  s <- runif(n, 0, A)                 
  C <- A - s                         
  event <- as.integer(T <= C)
  time  <- ifelse(event==1, T, C)
  data.frame(time=time, event=event)
}
