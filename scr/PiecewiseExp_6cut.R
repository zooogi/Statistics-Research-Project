library(rstan)
library(dplyr)
library(ggplot2)
library(patchwork)  
library(bayesplot)
source("scr/generate_fake_data_piecewise_kcut.R")

# ---- data: time=stage, event indicator=event (1=event, 0=deletion) ----
df <- read.csv("data/raw_data_turnover.csv")




stan_pwe_nocov <- "
data {
  int<lower=1> N;              // 个体数
  int<lower=1> J;              // 区间数
  matrix[N, J] tau;            // 暴露时间: overlap([0,y_i], (c_{j-1},c_j])
  int<lower=0,upper=1> d[N,J]; // 事件指示: d[i,j]=1 若 i 在区间 j 发生事件
}
parameters {
  vector[J] log_lambda;        // 每段基线 hazard 的 log
}
transformed parameters {
  vector[J] lambda = exp(log_lambda);
}
model {
  // 弱信息先验
  log_lambda ~ normal(0, 2);

  // 对数似然：sum_{i,j} [ d_ij * log λ_j  -  τ_ij * λ_j ]
  for (i in 1:N){
    for (j in 1:J){
      target += d[i,j] * log_lambda[j]
                - tau[i,j] * lambda[j];
    }
  }
}
generated quantities {
  // 个体 log-lik（供 LOO/WAIC）
  vector[N] log_lik;
  for (i in 1:N){
    real s = 0;
    real e = 0;
    for (j in 1:J){
      s += - tau[i,j] * exp(log_lambda[j]);
      e +=   d[i,j] * log_lambda[j];
    }
    log_lik[i] = e + s;
  }
}
"

# 分段切点（最后一个建议 ≥ max(stag) 或 ≥ A）
cuts <- c(40, 75, 105, 120, 160, 200)
J <- length(cuts)
c0 <- c(0, head(cuts, -1))
c1 <- cuts


y   <- df$stag
del <- df$event
N   <- length(y)

tau <- matrix(0, N, J)
d   <- matrix(0L, N, J)

for (j in 1:J){
  a <- c0[j]; b <- c1[j]
  tau[, j] <- pmax(0, pmin(y, b) - a)                     # 暴露时间
  d[, j]   <- as.integer(del == 1 & y > a & y <= b)       # 事件落在哪段
}

stan_data <- list(N=N, J=J, tau=tau, d=d)


sm  <- stan_model(model_code = stan_pwe_nocov)
fit <- sampling(sm, data=stan_data, chains=4, iter=2000, warmup=1000, seed=1)

# 简要summary
print(fit, pars="lambda", probs=c(.025,.5,.975))

# 详细summary矩阵
sum_all <- summary(fit)$summary
lambda_sum <- sum_all[grep("^lambda\\[", rownames(sum_all)),
                      c("mean","sd","2.5%","50%","97.5%","Rhat","n_eff")]
lambda_sum

# traceplot / 密度
posterior <- as.array(fit)
mcmc_trace(posterior, pars = grep("^lambda\\[", dimnames(posterior)$parameters, value = TRUE))
mcmc_dens_overlay(posterior, pars = grep("^lambda\\[", dimnames(posterior)$parameters, value = TRUE))





#------------ppc model checking------------------------------------

# 取一组后验样本
post_lambda <- rstan::extract(fit, pars="lambda")$lambda   # draws x J
A_fixed <- 179
id <- sample(seq_len(nrow(post_lambda)), 1)
lam_draw <- as.numeric(post_lambda[id, ])
stopifnot(all(is.finite(lam_draw)), all(lam_draw>0))

# 生成假数据
fake <- simulate_pwe(N, lam_draw, cuts, A_fixed, seed=253)

# ECDF 对比
df_plot <- bind_rows(
  data.frame(time=df$stag, event=df$event, source="real"),
  data.frame(time=fake$time, event=fake$event, source="fake")
)


p_event <- ggplot(df_plot %>% filter(event==1),
                  aes(x=time, color=source)) +
  stat_ecdf(size=1) +
  scale_color_manual(values=c("real"="blue", "fake"="red")) +
  labs( y="ECDF", x="Months (Event)") +
  theme_minimal(base_size = 16) +
  theme(legend.position="top")

p_cens <- ggplot(df_plot %>% filter(event==0),
                 aes(x=time, color=source)) +
  stat_ecdf(size=1) +
  scale_color_manual(values=c("real"="blue", "fake"="red")) +
  labs( y="ECDF", x="Months (Censored)") +
  theme_minimal(base_size = 16) +
  theme(legend.position="top")

p_event; p_cens
comb_piece6 <- (p_event + p_cens ) +
  plot_layout(guides = "collect", widths = c(1,1)) &
  theme(legend.position = "top")
comb_piece6
ggsave("images/piece6_ecdf_comparison.pdf",
       plot   = comb_piece6,
       width  = 7,   
       height = 4,   
       device = "pdf") 

