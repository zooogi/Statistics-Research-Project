source("scr/A_post_contour&marginal.R")


#### A prior ###############

## ---- Truncated Normal prior:  A ~ N(mu, sigma^2) truncated at [a_min, +inf) ----
## 默认把截断点设为 ymax，均值锁定在 ymax + 50，标准差 20
logprior_A_truncnorm_scalar <- function(A,
                                        mu    = ymax + 50,
                                        sigma = 20,
                                        a_min = ymax) {
  if (sigma <= 0) stop("sigma must be > 0")
  # 截断支持：A < a_min 时密度为 0（log 密度为 -Inf）
  if (A < a_min) return(-Inf)
  
  # 分子：未截断正态密度 log φ((A-mu)/sigma) - log sigma
  log_num <- dnorm(A, mean = mu, sd = sigma, log = TRUE)
  
  # 分母：P(Z >= a_min) = 1 - Φ((a_min-mu)/sigma)，用 log.p 避免精度问题
  # pnorm(lower.tail = FALSE) 直接给出 log[1 - Φ(.)]
  log_den <- pnorm(a_min, mean = mu, sd = sigma,
                   lower.tail = FALSE, log.p = TRUE)
  
  # 截断密度：log_num - log_den
  return(log_num - log_den)
}

## 向量化版本，方便在网格上调用
logprior_A_truncnorm <- Vectorize(function(A, mu = ymax + 50, sigma = 20, a_min = ymax)
  logprior_A_truncnorm_scalar(A, mu = mu, sigma = sigma, a_min = a_min)
)



lp_mat_tn <- matrix(NA_real_, nrow = nL, ncol = nA)
for (i in seq_len(nL)) {
  l <- lam_grid[i]
  for (j in seq_len(nA)) {
    A <- A_grid[j]
    lp <- loglik_exp_uniform(l, A, y, delta) +
      logprior_lambda_scalar(l, shape = alpha, rate = beta) +
      logprior_A_truncnorm_scalar(A, mu = ymax + 50, sigma = 20, a_min = ymax)
    lp_mat_tn[i, j] <- lp
  }
}

## 归一化
lp_vec_tn     <- as.vector(lp_mat_tn)
m_tn          <- max(lp_vec_tn)
post_unn_tn   <- exp(lp_vec_tn - m_tn)
Z_tn          <- sum(post_unn_tn) * cell_area
post_vec_tn   <- post_unn_tn / Z_tn
post_mat_tn   <- matrix(post_vec_tn, nrow = nL, ncol = nA)




## =========================
## 2) 计算 lambda 的边际
## =========================
p_lambda_tn  <- rowSums(post_mat_tn) * dA   # Trunc-Normal(A) 的结果

## 轻微数值漂移校正（可选）：保证积分≈1
p_lambda_tn  <- p_lambda_tn  / (sum(p_lambda_tn)  * dlam)
cat(" Trunc-Normal ≈", sum(p_lambda_tn) * dlam, "\n")


## =========================
## 3) 组装到同一张图并绘制
## =========================

lambda_df2 <- tibble(
  lambda = lam_grid,
  `Grid posterior (with A, Uniform)` = p_lambda_uni,
  `Grid posterior (with A, Trunc-Normal)` = p_lambda_tn
) |>
  pivot_longer(
    cols = c(`Grid posterior (with A, Uniform)`,
             `Grid posterior (with A, Trunc-Normal)`),
    names_to = "source", values_to = "density"
  ) |>
  mutate(
    source = factor(
      source,
      levels = c("Grid posterior (with A, Uniform)",
                 "Grid posterior (with A, Trunc-Normal)")
    )
  )


#plot
p_lambda_cmp2 <- ggplot(lambda_df2, aes(x = lambda, y = density,
                                        colour = source, linetype = source)) +
  geom_line(aes(linewidth = source)) +
  scale_linewidth_manual(values = c(
    "Grid posterior (with A, Uniform)"      = 1.6,
    "Grid posterior (with A, Trunc-Normal)" = 0.7
  )) +
  scale_color_manual(values = c(
    "Grid posterior (with A, Uniform)"      = "#377EB8",  # 蓝（均匀先验）
    "Grid posterior (with A, Trunc-Normal)" = "#4DAF4A"   # 绿（截断正态）
  )) +
  scale_linetype_manual(values = c(
    "Grid posterior (with A, Uniform)"      = "dashed",
    "Grid posterior (with A, Trunc-Normal)" = "solid"
  )) +
  labs(x = expression(lambda), y = "posterior density") +
  theme_bw(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.justification = "center")

p_lambda_cmp2

ggsave("images/diff_A_prior_marginal_compare.pdf",
       plot = p_lambda_cmp2, width = 7, height = 4, device = "pdf")
