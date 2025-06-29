library(survival) 
# 就已经把 veteran 自动加载进来了
#already automatically loaded in data:veteran
veteran<-veteran
head(veteran)

library(brms)
get_prior(time | cens(status) ~ 1, data = veteran, family = exponential())
fit <- brm(time | cens(status) ~ 1, data = veteran, family = exponential())
summary(fit)
plot(fit)
posterior_samples(fit)
#查看 Stan 代码（理解黑箱）
stancode(fit)

##加入所有的协变量尝试了一下
get_prior(time | cens(status) ~ ., data = veteran, family = exponential())
fit_all <- brm(time | cens(status) ~ ., data = veteran, family = exponential())
summary(fit_all)
plot(fit_all)
posterior_samples(fit_all)
#查看 Stan 代码（理解黑箱）
stancode(fit_all)
