df <- read.csv("data/raw_data_turnover.csv")
table(df$event)
summary(df$stag)
#是否有缺失值
sum(is.na(df))  
library(bayesplot)
library(ggplot2)
ggplot(df, aes(x = stag, fill = factor(event))) +
  geom_histogram(position = "identity", alpha = 0.6, binwidth = 3) +
  scale_fill_manual(values = c("orange", "dodgerblue"),
                    labels = c("Censored", "Event Occurred")) +
  labs(title = "Duration Histogram: Censored vs. Event",
       x = "Duration",
       y = "Count",
       fill = "Event Status") +
  theme_minimal()
#save
ggsave("images/separate_hist.png", width = 6, height = 4)



#————————————fit model————————————————
stan_code <- "
data {
  int<lower=0> N;
  vector<lower=0>[N] y;
  int<lower=0,upper=1> status[N];
}
parameters {
  real<lower=0> lambda;
}
model {
// Weakly informative Gamma prior
  lambda ~ gamma(0.001, 0.001);
  // Likelihood
  for (i in 1:N) {
    if (status[i] == 1) {
      target += exponential_lpdf(y[i] | lambda);
    } else {
      target += exponential_lccdf(y[i] | lambda);
    }
  }
}
"
stan_turnover <- list(
  N     = nrow(df),
  y     = df$stag,       
  status = df$event           # 0/1 向量
)

#fit
fit_stan_turnover <- stan(
  model_code = stan_code,
  data       = stan_turnover,
  chains     = 4,
  iter       = 1000,
  warmup     = 50,
  refresh    = 0
)

#————————————————————————————————————————
##看看fit的结果
summary(fit_stan_turnover)$summary  # 查看后验均值、标准差、Rhat 等
# —— traceplot —— 
bayesplot::mcmc_trace(
  as.array(fit_stan_turnover, pars = "lambda"),
  pars = "lambda"
)
ggsave("images/turnover_exp_model_traceplot.png", width = 6, height = 4)
#取后验采样的lambda
post_lam_turnover_stan <- extract(fit_stan_turnover,"lambda")$lambda
#画出后验采样的lambda直方图
hist(post_lam_turnover_stan,
     breaks = 30,           
     prob   = TRUE,       
     col    = "#1f77b4",    
     border = "white",
     xlab   = expression(lambda),
     main   = "posterior of lambda",
     ylab   = "Density")
ggsave("images/turnover post lambda.png", width = 6, height = 4)

#计算后验预测
d <- sum(df$event == 1)  # 真实的 event 个数
# 对每一个 lambda 后验样本，生成一个预测样本
yrep <- replicate(n = length(post_lam_turnover_stan), rexp(n = d, rate = post_lam_turnover_stan))
summary(as.vector(yrep))
summary(df$stag)

#画出对比图（事件数据的）
bayesplot::ppc_dens_overlay(y = df$stag[df$event == 1], yrep = t(yrep[ , 1:100])) + 
  ggplot2::xlim(0, 500)+  
  ggplot2::scale_color_manual(
  name = "Data Type", 
  values = c("y" = "#000000", "yrep" = "#3182bd"),
  labels = c("y" = "Observed event", "yrep" = "Posterior Predictive")
)+ 
  ggplot2::theme(
    plot.margin = margin(10, 0, 10, 10),  # 左右边距，右边多留一点空
    legend.position.inside = c(3, 0.75),
    legend.text = element_text(size = 9), # 图例字体变小
    legend.title = element_text(size = 8)
  ) 
ggsave("images/turnover_post_predic_event.png", width = 6, height = 4)
