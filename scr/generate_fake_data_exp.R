generate_data_fixed_lambda <- function(n = 1000,
                                       lambda = 0.05,  
                                       a      = -100,  # starting time T lower bound
                                       seed   = 42) {
  
  set.seed(seed)
  
  # latent duration
  y  <- rexp(n, rate = lambda)
  
  # entry time relative to "now = 0"
  T  <- runif(n, min = a, max = 0)
  
  # leaving time
  t  <- T + y
  
  # censoring rule
  event <- as.integer(t < 0)        #T<0 generates logical vectors: TRUE (event observation)/False (censoring).
  #As. intege() converts True/False to 1/0,，
  time  <- ifelse(event == 1, y, -T)   # observed time
  
  data.frame(time = time, event = event)
}