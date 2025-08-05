#-----------run function----------
source("scr/generate_fake_data_exp.R") 
source("scr/exponential model.R")
#-------------------package---------
library(dplyr)
library(ggplot2)
library(patchwork)  
#------get post mean-------------------
sum_post_lam_turnover <- summary(fit_stan_turnover, pars = "lambda")$summary
post_mean_lam_turnover <- sum_post_lam_turnover["lambda", "mean"]
#---------------- generate_data (via post mean)-------
fake_postpred_modelcheck_exp<- generate_data_fixed_lambda(n = 1000, lambda = post_mean_lam_turnover,
                                                          a = -160, seed = 25)
#----------------summary------------------
summary(fake_postpred_modelcheck_exp)

#----------------变成数据框合并
real <- df %>%
  transmute(source = "real",
            Y = stag,            # 观测时长 Y = min(T, C)
            delta = event)       # 事件指示

fake <- fake_postpred_modelcheck_exp %>%
  transmute(source = "fake",
            Y = time,  # 统一成正的观测时长
            delta = event)
# 合并
dat_all <- bind_rows(real, fake)

#----------plot----------------
# ECDF（事件子样本）
ggplot(filter(dat_all, delta == 1),
       aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x = "Y | delta = 1 (event times)",
       y = "ECDF",
       title = "ECDF of event times (conditional on being uncensored)") +
  theme_bw()
#save
ggsave("images/ppc_event_ecdf_exp.png", width = 6, height = 4)

#ECDF（删失子样本）
ggplot(filter(dat_all, delta == 0),
        aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x="Y | delta = 0 (censored durations)",
       y = "ECDF",
       title = "ECDF of censored durations (conditional on censoring)") +
  theme_bw()
#save
ggsave("images/ppc_censored_ecdf_exp.png", width = 6, height = 4)


