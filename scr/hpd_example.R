# ==== Packages ====
library(ggplot2)
library(mvtnorm)   # dmvnorm, dmvt
library(sn)        # dmvSN (skew-normal)
library(patchwork) # 拼图

# ==== 辅助函数：给定覆盖率 alpha，计算 HPD 阈值（密度水平）====
hpd_threshold <- function(z, dx, dy, alpha) {
  ord <- order(z, decreasing = TRUE)
  cum <- cumsum(z[ord]) * dx * dy
  z[ord][which(cum >= alpha)[1]]
}

# 生成网格 + 密度 + HPD 等高线图
make_hpd_plot <- function(dens_fun, pars, xlim = c(-4, 4), ylim = c(-4, 4),
                          n = 300, title = "", subtitle = "") {
  x <- seq(xlim[1], xlim[2], length.out = n)
  y <- seq(ylim[1], ylim[2], length.out = n)
  dx <- diff(x)[1]; dy <- diff(y)[1]
  
  # 计算密度（dens_fun 需返回向量化的密度）
  XY <- cbind(rep(x, each = n), rep(y, times = n))
  z  <- dens_fun(XY, pars)
  df <- data.frame(x = XY[,1], y = XY[,2], z = z)
  
  # 计算 HPD 阈值：39.3% 与 86.5%
  level_0393 <- hpd_threshold(z, dx, dy, 0.393)
  level_0865 <- hpd_threshold(z, dx, dy, 0.865)

  # 画图：用等高线标出 HPD 边界
  p <- ggplot(df, aes(x, y, z = z)) +
    ## 先画两个 HPD 边界：86.5%（外圈，虚线），39.3%（内圈，实线）
    geom_contour(breaks = c(level_0865), linewidth = 0.6, linetype = "dashed") +
    geom_contour(breaks = c(level_0393), linewidth = 0.9) + 
    # 中心区域微填充
    stat_contour_filled(breaks = c(level_0865, level_0393, max(z) + 1e99),
                        show.legend =  TRUE) +
    ## 两个 band 的颜色：外圈浅灰，内圈深灰；
    scale_fill_manual(
      name   = "HPD bands",
      values = c("grey80", "grey40"),       # 浅灰(外) → 深灰(内)
      labels = c("86.5% HPD", "39.3% HPD")
    ) +
    labs(x = expression(x[1]), y = expression(x[2])) +
    theme_minimal(base_size = 16) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.position = "right" )
  return(p)
}

# ==== 密度函数：二维正态 ====
dens_bvnormal <- function(XY, pars) {
  mvtnorm::dmvnorm(XY, mean = pars$mu, sigma = pars$Sigma)
}

# ==== 密度函数：二维 skew-normal（偏态）====
dens_bvskewnorm <- function(XY, pars) {
  # sn::dmsn 需要 xi (位置), Omega (协方差), alpha (形状)
  sn::dmsn(XY, xi = pars$xi, Omega = pars$Omega, alpha = pars$alpha)
}

# ==== 密度函数：二维 t（重尾）====
dens_bvt <- function(XY, pars) {
  nu <- pars$df
  # 目标协方差：pars$Sigma
  # 传给 dmvt 的 sigma 需要缩放到：Sigma * ((nu - 2)/nu)
  sigma_mvtnorm <- pars$Sigma * ((nu - 2) / nu)
  mvtnorm::dmvt(XY, delta = pars$mu, sigma = sigma_mvtnorm, df = nu, log = FALSE)
}


# -------------------------------
# 左图：二维正态（近似正态后验的“基准”）
# -------------------------------
pars_norm <- list(
  mu    = c(0, 0),
  Sigma = matrix(c(1, 0,
                   0, 1.0), 2, 2, byrow = TRUE)
)

p_left <- make_hpd_plot(
  dens_fun = dens_bvnormal,
  pars     = pars_norm,
  xlim     = c(-4, 4),
  ylim     = c(-4, 4),
  n        = 300,
  title    = "Normal posterior (benchmark)",
  subtitle = "HPD contours at 39.3% (solid) and 86.5% (dashed)"
)

# ---------------------------------
# 右图：偏态/重尾后验（任选其一）
# 这里默认用“偏态” skew-normal；若想改用重尾 t，请见下方替代
# ---------------------------------
use_skew <- FALSE  # TRUE 用 skew-normal；FALSE 用 t 重尾

if (use_skew) {
  pars_skew <- list(
    xi    = c(0, 0),
    Omega = matrix(c(1, 0,
                     0, 1), 2, 2, byrow = TRUE),
    alpha = c(6, 0)  # 形状参数，产生明显偏态
  )
  p_right <- make_hpd_plot(
    dens_fun = dens_bvskewnorm,
    pars     = pars_skew,
    xlim     = c(-4, 4),
    ylim     = c(-4, 4),
    n        = 300,
    title    = "Skewed posterior (illustrative)",
    subtitle = "HPD contours at 39.3% (solid) and 86.5% (dashed)"
  )
} else {
  pars_t <- list(
    mu    = c(0, 0),
    # 直接用与正态相同的协方差作为“目标协方差”
    Sigma = pars_norm$Sigma,
    df    = 2.5
  )
  p_right <- make_hpd_plot(
    dens_fun = dens_bvt,
    pars     = pars_t,
    xlim     = c(-4, 4),
    ylim     = c(-4, 4),
    n        = 300,
    title    = "Heavy-tailed posterior (t, df=2.5)",
    subtitle = "HPD contours at 39.3% (solid) and 86.5% (dashed)"
  )
}

p <- p_left + p_right +
  plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "right")
print(p)


# 保存
ggsave("images/hpd_normal_vs_student-t.pdf",
       p, width = 8.5, height = 4.2, device = "pdf")
