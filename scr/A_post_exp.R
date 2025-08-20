library(ggplot2)
library(tidyr)
library(ggnewscale)
library(dplyr)
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



#-----------------------选择网格法的网格范围！------------

## ==== lambda的范围，使用精确解的分位数对应的数值====
alpha <- 0.001; beta <- 0.001
#公式20（paper）
shape_post <- sum(delta) + alpha
rate_post  <- sum(y) + beta
lam_lo <- qgamma(0.001, shape_post, rate_post)   # ~ 0.1% 分位
lam_hi <- qgamma(0.999, shape_post, rate_post)   # ~ 99.9% 分位

## A 的范围
A_lo <- ymax + 1e-6
A_hi <- ymax + 500

## 网格尺寸
nL <- 220
nA <- 220
lam_grid <- seq(lam_lo, lam_hi, length.out = nL)
A_grid   <- seq(A_lo,   A_hi,   length.out = nA)





## ==== ####################计算对数后验 ====


## 注意：把先验和似然分开写，便于检查
lp_mat <- matrix(NA_real_, nrow = nL, ncol = nA)
for (i in seq_len(nL)) {
  l <- lam_grid[i]
  ## 先把 λ 的先验和与 λ 有关的似然部分求出来以复用
  for (j in seq_len(nA)) {
    A <- A_grid[j]
    lp <- loglik_exp_uniform(l, A, y, delta) +
      logprior_lambda_scalar(l, shape=alpha, rate=beta) +
      logprior_A_scalar(A, a_min=A_lo, a_max=A_hi)
    lp_mat[i, j] <- lp
  }
}





## ==== 归一化离散后验密度（数值稳定的 log-sum-exp）====
dlam <- diff(lam_grid)[1]
dA   <- diff(A_grid)[1]
cell_area <- dlam * dA   # 等距网格：单元面积常数

lp_vec <- as.vector(lp_mat)

m      <- max(lp_vec)
post_unnorm <- exp(lp_vec - m)          # 避免下溢
Z <- sum(post_unnorm) * cell_area       # 近似积分的归一化常数
post_vec <- post_unnorm / Z             # 离散化后的“密度”（带面积后已归一）






#-----------为画图做一个数据框--------------------


## 恢复到矩阵 & 数据框
post_mat <- matrix(post_vec, nrow = nL, ncol = nA)
grid_df <- expand.grid(lambda = lam_grid, A = A_grid) %>%
  mutate(p = as.vector(post_mat))



## ==== 求 2D HPD 阈值：39.3% & 86.5% ====
hpd_thresholds <- function(pmat, probs = c(0.393, 0.865)){
  v <- sort(as.vector(pmat), decreasing = TRUE)
  cs <- cumsum(v) / sum(v)
  sapply(probs, function(q){
    idx <- which(cs >= q)[1]
    v[idx]
  })
}
levs <- hpd_thresholds(post_mat, probs = c(0.393, 0.865))
names(levs) <- c("HPD39.3", "HPD86.5")





## ====MAP 点 ====
which_max <- which.max(post_vec)
i_map <- ((which_max - 1) %% nL) + 1
j_map <- ((which_max - 1) %/% nL) + 1
map_point <- data.frame(lambda = lam_grid[i_map], A = A_grid[j_map])






## ====plot（ HPD 等高线 + MAP）==================

#下面的这个是为A的画图范围好看点控制的！
# A 的边际（对 lambda 求和）
pA <- colSums(post_mat) * dlam        

# 找 0.1% 到 99.9% 分位的 A 区间（自动“放大”到有质量的区域）
cdfA <- cumsum(pA) * dA
A_lo_zoom <- approx(cdfA, A_grid, xout = 0.001)$y
A_hi_zoom <- approx(cdfA, A_grid, xout = 0.999)$y
A_lo_zoom; A_hi_zoom
#---------------------------------------

############# plot ##################

post_hpd <- ggplot(grid_df, aes(x = lambda, y = A)) +
  new_scale_fill() +
  # 填充 HPD 区域
  geom_contour_filled(aes(z = p), breaks = c(sort(levs), Inf), alpha = 0.4) +
  scale_fill_brewer(palette = "Greys", name = "HPD bands", direction = -1) +
  ## 只画 contour 线
  geom_contour(aes(z = p), breaks = sort(levs), colour = "black", linewidth = 0.6) +
  #geom_hline(yintercept = ymax, linetype = "dashed", colour = "grey40") +
  # MAP 点
  geom_point(data = map_point, aes(x = lambda, y = A),
             inherit.aes = FALSE, colour = "white", fill = "black",
             shape = 21, size = 2.2, stroke = 0.4) +
  annotate("text", x = map_point$lambda, y = map_point$A,
           label = "MAP", vjust = 1.8, size = 3.2) +
  scale_fill_manual(
    name   = "HPD bands",
    values = c("grey75", "grey30"),          # 外圈浅、内圈深
    labels = c("86.5% HPD", "39.3% HPD")
  ) +
  coord_cartesian(ylim = c(A_lo_zoom, A_hi_zoom)) +
  labs(x = expression(lambda), y = "A")+
  theme_bw()


ggsave("images/post_contour.pdf",
       plot   = post_hpd,
       width  = 6,   
       height = 4,   
       device = "pdf") 

#——-------------——————————————————----——————————






#-------------------lambda边际后验的对比-------------------

#没有A的模型的lambda后验解析解
d <- sum(df$event == 1)
alpha<-0.001
sum_y  <- sum(df$stag)
beta<-0.001
analytic_dens <- dgamma(lam_grid, shape =  d + alpha, rate = sum_y+beta)



## ==== 数值边际：沿 A 方向求和 ====
p_lambda <- rowSums(post_mat) * dA                   # Riemann 求和
# 规范化检查（应接近 1）
sum_lambda <- sum(p_lambda) * dlam
cat("Integral of p_lambda ≈", sum_lambda, "\n")


#变成数据框用于画图
lambda_df <- tibble(
  lambda        = lam_grid,
  `Grid posterior (with A)` = p_lambda,
  `Analytic posterior (no A)` = analytic_dens
) |>
  pivot_longer(cols = c(`Grid posterior (with A)`, `Analytic posterior (no A)`),
               names_to = "source", values_to = "density")


#plot
p_lambda_cmp <- ggplot(lambda_df, aes(x = lambda, y = density,
                                      colour = source, linetype = source)) +
  geom_line(linewidth = 0.9) +
  labs(x = expression(lambda), y = "posterior density") +
  geom_line(aes(linewidth = source)) +
  scale_linewidth_manual(values = c("Grid posterior (with A)" = 0.5,
                                    "Analytic posterior (no A)" = 1.5))+
  scale_color_manual(values = c("Grid posterior (with A)" = "#377EB8",
                                "Analytic posterior (no A)" = "#E41A1C")) +
  scale_linetype_manual(values = c("Grid posterior (with A)" = "solid",
                                   "Analytic posterior (no A)" = "dashed")) +
  theme_bw(base_size = 17) + theme(legend.title = element_blank())+
  theme(
    legend.position = "top",       # 或者 "bottom"
    legend.justification = "center"
  )
p_lambda_cmp

ggsave("images/lambda_marginal_compare.pdf",
       plot = p_lambda_cmp, width = 6, height = 4, device = "pdf")


#-----------------------------------------------





###计算95%CrI-------------------

wquantile <- function(grid, pmass, probs){
  o <- order(grid) #将x轴范围升序排序
  g <- grid[o];  #升序后的网格点
  w <- pmass[o] / sum(pmass) #对应的归一化权重
  cdf <- cumsum(w) #累积分布函数（CDF）
  approx(x = cdf, y = g, xout = probs, ties = "ordered", rule = 2)$y
}


## 4) 计算等尾 95% CrI 与中位数
lambda_ci   <- wquantile(lam_seq, post_lambda, c(0.025, 0.975))
lambda_medi <- wquantile(lam_seq, post_lambda, 0.5)
A_ci        <- wquantile(A_seq,   post_A,      c(0.025, 0.975))
A_medi      <- wquantile(A_seq,   post_A,      0.5)

lambda_ci
lambda_medi


