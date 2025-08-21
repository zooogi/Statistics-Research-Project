## ----- 0) 准备数据：你已有 df -----
df <- read.csv("data/raw_data_turnover.csv")
## df$stag = time (months), df$event = 1(event) / 0(censored)
y     <- df$stag
delta <- df$event
n     <- length(y)
ymax  <- max(y)



##########################似然函数！######################


## ----- 1) 对数似然函数 logLik(lambda, A) -----
loglik_exp_uniform <- function(lambda, A, y, delta) {
  y      <- as.numeric(y)
  delta  <- as.integer(delta)
  ymax   <- max(y)
  
  # 基本约束：参数域 + 指示函数支持集
  if (lambda <= 0) return(-Inf)
  if (A < ymax)    return(-Inf)               # 保证 0 <= y_i <= A
  #if (any(delta == 1 & y >= A)) return(-Inf)  # 严格小于
  
  # 对数似然
  # 事件项: log[ lambda * exp(-lambda y) * (1 - y/A) ]
  # 删失项: log[ exp(-lambda y) / A ]
  log_event <- delta       * (log(lambda) - lambda * y + log1p(- y / A))  # log(1 - y/A)
  log_cens  <- (1 - delta) * (- lambda * y - log(A))
  sum(log_event + log_cens)
}





#############################先验分布




#设置先验分布 lambda是gamma，A是uniform


## 弱信息先验
## lambda ~ Gamma(shape, rate)
logprior_lambda_scalar <- function(l, shape = 0.001, rate = 0.001) {
  if (l <= 0) return(-Inf)
  dgamma(l, shape = shape, rate = rate, log = TRUE)  # 示例：弱信息 Gamma
}
## A ~ Uniform(a_min, a_max)，下界必须 ≥ ymax
logprior_A_scalar <- function(A, a_min = ymax + 1e-6, a_max = ymax + 500) {
  if (a_min >= a_max) stop("a_min must be < a_max")
  if (A < a_min || A > a_max) return(-Inf)
  -log(a_max - a_min)  # Uniform(a_min,a_max) 的对数密度
}

##用 Vectorize 包装成“向量化版”
logprior_lambda <- Vectorize(logprior_lambda_scalar)
logprior_A <- Vectorize(function(A) logprior_A_scalar(A))


##############################################################