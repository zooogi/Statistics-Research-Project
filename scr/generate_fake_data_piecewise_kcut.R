
# 单个个体：分段指数采样事件时长（无PH）
rpexp_1 <- function(lambda, cuts){
  c0 <- c(0, head(cuts, -1)); c1 <- cuts
  t <- 0
  for (j in seq_along(cuts)) {
    rate <- lambda[j]
    L    <- c1[j] - c0[j]
    if (!is.finite(rate) || rate <= 0) { t <- t + L; next }
    w <- rexp(1, rate)
    if (w < L) return(t + w)   # 在该段内发生
    t <- t + L                 # 跨段
  }
  t
}

# 批量生成：窗口 A 内的观测 (行政删失)
simulate_pwe <- function(n, lambda, cuts, A, seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  T <- replicate(n, rpexp_1(lambda, cuts))
  s <- runif(n, 0, A)                 # 入样时间
  C <- A - s                          # 行政删失阈值
  event <- as.integer(T <= C)
  time  <- ifelse(event==1, T, C)
  data.frame(time=time, event=event)
}
