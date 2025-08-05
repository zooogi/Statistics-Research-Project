df <- read.csv("data/raw_data_turnover.csv")
table(df$event)
summary(df$stag)
length(df$event)
#是否有缺失值
sum(is.na(df))  
library(bayesplot)
library(rstan)
library(ggplot2)

ggplot(df, aes(x = stag)) +
  geom_histogram(binwidth = 3, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ event, labeller = labeller(event = c("0" = "Censored", "1" = "Event Occurred"))) +
  labs(title = "Duration Histogram by Event Status",
       x = "Duration",
       y = "Count") +
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
  seed = 2,
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

d <- sum(df$event == 1)
alpha<-0.001
sum_y  <- sum(df$stag)
beta<-0.001
# 用一个网格来评估解析密度
lambda_grid  <- seq(min(post_lam_turnover_stan), max(post_lam_turnover_stan), length = 1000)
analytic_dens <- dgamma(lambda_grid, shape =  d + alpha, rate = sum_y+beta)

#画出后验采样的lambda直方图
hist(post_lam_turnover_stan,
     breaks = 30,           
     prob   = TRUE,       
     col    = "#1f77b4",    
     border = "white",
     xlab   = expression(lambda),
     main   = "posterior of lambda",
     ylab   = "Density",
     ylim   = c(0, max(density(post_lam_turnover_stan)$y, analytic_dens)))
lines(lambda_grid, analytic_dens,
      lty = 2, lwd = 2, col = "red")
legend("topright",
       legend = c("MCMC posterior", "Analytic"),
       lty    = c(NA,2),
       pch    = c(15, NA),
       pt.cex = c(2.5, NA),
       col    = c("#1f77b4","red"),
       bty    = "n",
       inset  = c(-0.07, 0)
)
ggsave("images/turnover post lambda.png", width = 6, height = 4)




######下面可以删了，因为下面是事件数据的后验预测cdf

#计算后验预测

#------------ecdf vs cdf 事件数据event------------------
set.seed(123)
# 使用Monte Carlo近似
y_grid <- seq(0, 400, length.out = 500)
cdf_vals <- sapply(y_grid, function(t) {
  mean(pexp(t, rate = post_lam_turnover_stan))  # pexp 是 CDF
})

# 画图：event 数据的 ECDF + 后验预测的 CDF
ggplot() +
  stat_ecdf(data = df[df$event == 1, ],
            aes(x = stag), 
            color = "black", size = 1, geom = "step") +
  geom_line(data = data.frame(t = y_grid, cdf = cdf_vals),
            aes(x = t, y = cdf),
            color = "blue", size = 1.2) +
  labs(title = "Posterior Predictive CDF vs Observed ECDF",
       x = "Time", y = "Cumulative Probability") +
  theme_minimal()
ggsave("images/turnover_post_predic_event_ecdf.png", width = 6, height = 4)


