#-----------run function----------
source("scr/generate_fake_data_exp.R") 
source("scr/exponential model.R")
#-------------------package---------
library(dplyr)
#----------------
# 真实数据：df 有 event (1=事件/未删失, 0=删失), stag (观测时长, 已为正)
real <- df %>%
  transmute(source = "real",
            Y = stag,            # 观测时长 Y = min(T, C)
            delta = event)       # 事件指示

# 模拟数据：postpred_modelcheck_exp_fake 有 event, time
# 你算法里：event=1 时 time = y_i（正）；event=0 时 time = -T_i（负号存“已工作时长”）
fake <- fake_postpred_modelcheck_exp %>%
  transmute(source = "fake",
            Y = time,  # 统一成正的观测时长
            delta = event)
# 合并
dat_all <- bind_rows(real, fake)
#----------
ggplot(dat_all, aes(x = Y, colour = source)) +
  stat_ecdf(size = 0.9) +
  labs(x = "observed time Y = min(T, C)", y = "ECDF",
       title = "Observed-level comparison (ECDF of all Y)") +
  theme_bw()

dat_all %>%
  group_by(source) %>%
  summarise(censor_rate = mean(delta == 0),
            n = n())

# 直方图（删失子样本）
ggplot(filter(dat_all, delta == 0),
       aes(x = Y, fill = source)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 100) +
  labs(x = "Y given censored (delta=0)", y = "count",
       title = "Right-censored durations (histogram)") +
  theme_bw()

# ECDF（删失子样本）
ggplot(filter(dat_all, delta == 0),
       aes(x = Y, colour = source)) +
  stat_ecdf(size = 0.9) +
  labs(x = "Y given event (delta=1)", y = "ECDF",
       title = "event durations (ECDF)") +
  theme_bw()



#------get post mean-------------------
sum_post_lam_turnover <- summary(fit_stan_turnover, pars = "lambda")$summary
post_mean_lam_turnover <- sum_post_lam_turnover["lambda", "mean"]
#---------------- generate_data (via post mean)-------
fake_postpred_modelcheck_exp<- generate_data_fixed_lambda(n = 1000, lambda = post_mean_lam_turnover,
                                                          a = -160, seed = 25)
#----------------summary------------------
summary(fake_postpred_modelcheck_exp)

#---------------------plot---------------------------
# 1. 变成数据框 --------------------------------------------
fake_event  <- fake_postpred_modelcheck_exp$time[fake_postpred_modelcheck_exp$event == 1]
fake_censor <- fake_postpred_modelcheck_exp$time[fake_postpred_modelcheck_exp$event == 0]
real_event  <- df$stag[df$event == 1]
real_censor <- df$stag[df$event == 0]
# 2. 画 2×2 面板 ---------------------------------------------------
par(mfrow = c(2, 2),
    mar   = c(4, 4, 2, 1))   # 调窄边距

## (1) RealEvent --------------------------------------------------
hist(real_event, breaks = 100,
     col = "steelblue", border = NA,
     main = "Real Event", xlab = "time")

## (2) FakeEvent --------------------------------------------------
dummy <- hist(real_event, breaks = 100, plot = FALSE)  # 固定坐标
hist(fake_event, breaks = dummy$breaks,
     col = rgb(0,0,1,0.4), border = NA,
     main = "Fake Event", xlab = "time", ylab = "count",
     ylim = c(0, max(dummy$counts, hist(fake_event, breaks = dummy$breaks, plot = FALSE)$counts)*1.2))

## (3) RealCensored ----------------------------------------------
hist(real_censor, breaks = 100,
     col = "steelblue", border = NA,
     main = "Real Censored", xlab = "time")

## (4) FakeCensored ----------------------------------------------
dummy2 <- hist(real_censor, breaks = 100, plot = FALSE)
hist(fake_censor, breaks = dummy2$breaks,
     col = rgb(0,0,1,0.4), border = NA,
     main = "Fake Censored", xlab = "time", ylab = "count",
     ylim = c(0, max(dummy2$counts, hist(fake_censor, breaks = dummy2$breaks, plot = FALSE)$counts)*1.2))
par(mfrow = c(1,1))      # 还原布局


