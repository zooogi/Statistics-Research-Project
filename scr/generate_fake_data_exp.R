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
  event <- as.integer(t < 0)          # t < 0 产生逻辑向量：TRUE(事件观测) / FALSE(删失)。
  #as.integer() 把 TRUE / FALSE 转成 1 / 0，
  time  <- ifelse(event == 1, y, -T)   # observed time
  
  data.frame(time = time, event = event)
}