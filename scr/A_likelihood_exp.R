## ----- 0) 准备数据：你已有 df -----
df <- read.csv("data/raw_data_turnover.csv")
## df$stag = time (months), df$event = 1(event) / 0(censored)
y     <- df$stag
delta <- df$event
n     <- length(y)
ymax  <- max(y)

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


#--------给等高线图准备数据---------



#逐点计算 l(lambda_i,A_j)，构成一个矩阵 Z。 z[i, j] 就是 l(lambda_i, A_j)
vec_loglik <- function(lam_seq, A_seq, y, delta) {
  #先建一个空矩阵 z 来装结果
  z <- matrix(NA_real_, nrow = length(lam_seq), ncol = length(A_seq))
  #双层循环，逐格点计算
  for (i in seq_along(lam_seq)) {
    for (j in seq_along(A_seq)) {
      z[i, j] <- loglik_exp_uniform(lam_seq[i], A_seq[j], y, delta)
    }
  }
  z
}



#-------------设置图的x，y轴范围----------------

## 找峰值 (λ̂, Â)，目的是让画图的范围可以自动缩小一点，好控制Aseq和lam seq的范围选择

###负对数似然
negloglik_uv <- function(uv){
  lambda <- exp(uv[1])
  A      <- ymax + exp(uv[2])
  -loglik_exp_uniform(lambda, A, y, delta)
}

lam0 <- sum(delta)/sum(y)
###得到极大似然点
fit0 <- optim(c(log(lam0), log(30)), negloglik_uv,
              method="BFGS", control=list(reltol=1e-10))
lam_hat <- exp(fit0$par[1])
A_hat   <- ymax + exp(fit0$par[2])
###设置lambda的画图范围
lam_seq <- exp(seq(log(lam_hat*0.5), log(lam_hat*1.6), length.out=120))
eps     <- 1e-3
A_low   <- ymax + eps
A_high  <- min(A_hat*1.5, ymax + 80)
A_high  <- max(A_high, A_low + 10)
###A的画图范围
A_seq   <- seq(A_low, A_high, length.out=150)

###----------run 计算log似然-----------
Z <- vec_loglik(lam_seq, A_seq, y, delta)
Z_rel <- Z - max(Z, na.rm=TRUE)

##-------------画图plot——————————————————
levels <- -c(0.5, 2, 4.5, 8, 12, 20)
contour(lam_seq, A_seq, Z_rel, levels = levels,
        xlab = expression(lambda), ylab="A (months)",
        main="Likelihood contours for (lambda, A)")
abline(h=ymax, lty=3, col="grey50")
points(lam_hat, A_hat, pch=19, col="red")
legend("topright",
       legend=sprintf("peak: lambda=%.4f, A=%.1f", lam_hat, A_hat),
       bty="n", text.col="red", cex=0.9)




#----_-_-------可选（ggplot）--------------------
grid_df <- expand.grid(lambda = lam_seq, A = A_seq)
grid_df$Z <- as.vector(Z_rel)


levels <- -c(0.5, 2, 4.5, 8, 12, 20)
ggplot(grid_df, aes(x=lambda, y=A, z=Z)) +
  geom_contour(breaks = levels, colour = "black") +
  geom_point(aes(x=lam_hat, y=A_hat), colour="red", size=2) +
  geom_hline(yintercept = ymax, linetype=3, colour="grey50") +
  annotate("text", x = lam_hat, y = A_hat + 4,
           label = sprintf("peak: λ=%.4f, A=%.1f", lam_hat, A_hat),
           colour="red", hjust=0.5, size=3.5) +
  labs(title="Likelihood contours for (lambda, A)",
       x=expression(lambda), y="A (months)") +
  theme_bw()
#-----------------------------------------------








##--------check数值--------
ind <- which(Z == max(Z, na.rm=TRUE), arr.ind=TRUE)
cat("grid peak at  lambda=", lam_seq[ind[1]], 
    ",  A=", A_seq[ind[2]], "\n")
summary(as.vector(Z))







## ----- 5) 画 profile likelihood（剖面似然），辅助解读可识别性 -----
## (a) A 的 profile：对每个 A，最大化 lambda
Sdelta <- sum(delta);  Sy <- sum(y)
lam_hat_closed <- Sdelta / Sy  # 闭式 MLE

profile_A <- sapply(A_seq, function(A) {
  loglik_exp_uniform(lam_hat_closed, A, y, delta)
})
profile_A <- profile_A - max(profile_A)

## (b) lambda 的 profile：对每个 lambda，最大化 A
eps <- 1e-8
profile_lam <- sapply(lam_seq, function(lam) {
  obj <- function(logdA) {
    A <- ymax + exp(logdA) + eps   # A>ymax
    -loglik_exp_uniform(lam, A, y, delta)
  }
  hi <- log((max(A_seq) - ymax) + 1)  # 稍宽的上界
  opt <- optimize(obj, interval = c(log(1e-8), hi))
  -opt$objective
})
profile_lam <- profile_lam - max(profile_lam)


oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(A_seq, profile_A, type = "l",
     xlab = "A (months)", ylab = "profile log-likelihood (rel.)",
     main = "Profile for A")

plot(lam_seq, profile_lam, type = "l",
     xlab = expression(lambda), ylab = "profile log-likelihood (rel.)",
     main = "Profile for lambda")
par(oldpar)
