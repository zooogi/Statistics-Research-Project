# Weibull(rate λ, shape k) 版的假数据生成
# R 的 rweibull 用 (shape=k, scale=θ)，与 rate 的换算是 θ = λ^(-1/k)
generate_data_fixed_weibull <- function(n = 1000,
                                        lambda = 0.05,   # rate
                                        k      = 1.4,    # shape
                                        a      = -100,   # 进入时间下界，窗口是 [a, 0]
                                        seed   = 42) {
  set.seed(seed)
  
  # 1) 潜在事件时长 T ~ Weibull(k, scale = λ^{-1/k})
  theta <- lambda^(-1/k)
  y <- rweibull(n, shape = k, scale = theta)
  
  # 2) 入样时间（相对“现在=0”）：均匀进入
  T_entry <- runif(n, min = a, max = 0)  # 负数
  
  # 3) 离开时间 = 入样 + 潜在时长
  t_leave <- T_entry + y
  
  # 4) 行政删失：到 0 为止
  event <- as.integer(t_leave < 0)       # 1=事件被观测；0=删失
  time  <- ifelse(event == 1, y, -T_entry)  # 观测到的时间
  
  data.frame(time = time, event = event)
}
