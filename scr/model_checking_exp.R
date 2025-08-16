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
N<- length(df$event)
fake_postpred_modelcheck_exp<- generate_data_fixed_lambda(n = N, lambda = post_mean_lam_turnover,
                                                          a = -200, seed = 25)
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


#ECDF（删失子样本）
ggplot(filter(dat_all, delta == 0),
        aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x="Y | delta = 0 (censored durations)",
       y = "ECDF",
       title = "ECDF of censored durations (conditional on censoring)") +
  theme_bw()
#save


ggplot(fake_postpred_modelcheck_exp, aes(x = time)) +
  geom_histogram(binwidth = 3, fill = "#D55E00", alpha = 0.7) +
  facet_wrap(~ event, labeller = labeller(event = c("0" = "Censored", "1" = "Event Occurred"))) +
  labs(
    x = "Months",
    y = "Count") +
  theme_minimal(base_size = 20)


######################################
###### try different a----------------
#---------------a=-30-------------------
fake_ppc_a30_exp<- generate_data_fixed_lambda(n = N, lambda = post_mean_lam_turnover,
                                               a = -30, seed = 295)
summary(fake_ppc_a30_exp)
table(fake_ppc_a30_exp$event)

fake_a30 <- fake_ppc_a30_exp %>%
  transmute(source = "fake",
            Y = time,  # 统一成正的观测时长
            delta = event)
# 合并
dat_a30 <- bind_rows(real, fake_a30)

#########-----------a=-1000----------------------
fake_ppc_a1000_exp<- generate_data_fixed_lambda(n = N, lambda = post_mean_lam_turnover,
                                              a = -1000, seed = 295)
summary(fake_ppc_a1000_exp)
table(fake_ppc_a1000_exp$event)

fake_a1000 <- fake_ppc_a1000_exp %>%
  transmute(source = "fake",
            Y = time,  # 统一成正的观测时长
            delta = event)
# 合并
dat_a1000 <- bind_rows(real, fake_a1000)




#----------plot----------------

#----------plot----------------
#a=-30-------------------
# ECDF（事件子样本）
base <- 15
p31<-ggplot(filter(dat_a30, delta == 1),
       aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x = "Y | delta = 1 (event times)",
       y = "ECDF",
       colour = "source"
  ) +
  theme_bw(base_size = base)
#save
ggsave("images/ppc_event_ecdf_A30.png", width = 6, height = 4, dpi = 300)

#ECDF（删失子样本）
p32<-ggplot(filter(dat_a30, delta == 0),
       aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x="Y | delta = 0 (censored durations)",
       y = "ECDF"
       ) +
  theme_bw(base_size = base)
#save
ggsave("images/ppc_censored_ecdf_A30.png", width = 6, height = 4, dpi = 300)


p33<-ggplot(fake_ppc_a30_exp, aes(x = time)) +
  geom_histogram(binwidth = 3, fill = "#D55E00", alpha = 0.7) +
  facet_wrap(~ event, labeller = labeller(event = c("0" = "Censored", "1" = "Event Occurred"))) +
  labs(
       x = "Months",
       y = "Count") +
  theme_minimal(base_size = 21)
p33
ggsave("images/fake_duration_hist_a30.pdf",plot = p33,device = "pdf", 
       width = 6, height = 4)

combo <- (p31 + p32 ) +
  plot_layout(guides = "collect", widths = c(1,1)) &
  theme(legend.position = "right")
combo
ggsave("images/ppc_two_a30.pdf",
       plot   = combo,
       width  = 7,   # 目标在论文中的实际宽度（英寸）
       height = 2.8,   # 合理高度，别太矮
       device = "pdf") 




#------------------------a=-1000-----------------
# ECDF（事件子样本）
p1k1<-ggplot(filter(dat_a1000, delta == 1),
       aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x = "Y | delta = 1 (event times)",
       y = "ECDF") +
  theme_bw(base_size = base)
#save
ggsave("images/ppc_event_ecdf_A1000.png", width = 6, height = 4, dpi = 300)


#ECDF（删失子样本）
p1k2<-ggplot(filter(dat_a1000, delta == 0),
       aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x="Y | delta = 0 (censored durations)",
       y = "ECDF") +
  theme_bw(base_size = base)
#save
ggsave("images/ppc_censored_ecdf_A1000.png", width = 6, height = 4, dpi = 300)

p1k3<-ggplot(fake_ppc_a1000_exp, aes(x = time)) +
  geom_histogram(binwidth = 3, fill = "#D55E00", alpha = 0.7) +
  facet_wrap(~ event, labeller = labeller(event = c("0" = "Censored", "1" = "Event Occurred"))) +
  labs(x = "Months",
       y = "Count") +
  theme_minimal(base_size = 21)
p1k3
ggsave("images/fake_duration_hist_a1000.pdf",plot = p1k3,device = "pdf", 
       width = 6, height = 4)


comb <- (p1k1 + p1k2 ) +
  plot_layout(guides = "collect", widths = c(1,1)) &
  theme(legend.position = "right")
comb
ggsave("images/ppc_two_a1000.pdf",
       plot   = comb,
       width  = 7,   # 目标在论文中的实际宽度（英寸）
       height = 2.8,   # 合理高度，别太矮
       device = "pdf") 
