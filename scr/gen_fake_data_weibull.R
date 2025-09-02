# Fake data generation in Weibull (rate λ, shape k) version
# The rweibull of R is expressed as (shape=k, scale=θ), 
# and the conversion to rate is θ=λ ^ (-1/k)
generate_data_fixed_weibull <- function(n = 1000,
                                        lambda = 0.05,   # rate
                                        k      = 1.4,    # shape
                                        a      = -100,  
                                        seed   = 42) {
  set.seed(seed)
  
  # 1)  T ~ Weibull(k, scale = λ^{-1/k})
  theta <- lambda^(-1/k)
  y <- rweibull(n, shape = k, scale = theta)
  
  # 2) Sample injection time (relative to "now=0"): uniform entry
  T_entry <- runif(n, min = a, max = 0) 
  
  # 3) Departure time=sample injection+potential duration
  t_leave <- T_entry + y
  
  # 4) Administrative deletion: up to 0
  event <- as.integer(t_leave < 0)     
  time  <- ifelse(event == 1, y, -T_entry) 
  
  data.frame(time = time, event = event)
}
