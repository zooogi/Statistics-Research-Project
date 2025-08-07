library(survival) 
library(rstan)
library(ggplot2)
# 就已经把 veteran 自动加载进来了
#already automatically loaded in data:veteran
veteran<-veteran
head(veteran)


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

# 准备 stan_data，一定要包含 status
stan_data <- list(
  N      = nrow(veteran),
  y      = veteran$time,
  status = veteran$status
)

# 运行 MCMC
fit_stan <- stan(
  model_code = stan_code,
  data       = stan_data,
  chains     = 4,
  iter       = 1000,
  warmup     = 50,
  refresh    = 0
)

# 提取后验
post_lambda_stan <- extract(fit_stan,"lambda")$lambda
y      <- veteran$time
d <- sum(veteran$status)
alpha<-0.001
sum_y  <- sum(y)
beta<-0.001
# 用一个网格来评估解析密度
lambda_grid_stan      <- seq(min(post_lambda_stan), max(post_lambda_stan), length = 1000)
analytic_density <- dgamma(lambda_grid_stan, shape =  d + alpha, rate = sum_y+beta)


# —— traceplot —— 
bayesplot::mcmc_trace(
  as.array(fit_stan, pars = "lambda"),
  pars = "lambda"
)
ggsave("images/veteran_traceplot.png", width = 6, height = 4)

hist(post_lambda_stan,
     breaks = 30,          
     prob   = TRUE,         
     col    = "#A0CBE8",   
     border = "white",
     xlab   = expression(lambda),
     ylab   = "Density",
     ylim   = c(0, max(density(post_lambda_stan)$y, analytic_density))
)
lines(lambda_grid_stan, analytic_density,
      lty = 2, lwd = 2, col = "red"
)
legend("topright",
       legend = c("MCMC histogram", "Analytic post."),
       lty    = c(NA,2),
       pch    = c(15, NA),
       pt.cex = c(2.5, NA),
       col    = c("#A0CBE8","red"),
       bty    = "n",
       inset  = c(-0.1, 0)
)
ggsave("images/veteran_post_lam.png", width = 6, height = 4)
