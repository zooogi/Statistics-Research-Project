source("scr/A_post_exp.R")
source("scr/generate_fake_data_exp.R") 
N<- length(df$event)
fake_A_post<- generate_data_fixed_lambda(n = N, lambda = lam_mean,
                                                     a = -A_mean, seed = 25)
#----------------summary------------------
summary(fake_A_post)

#----------------变成数据框合并
real_A <- df %>%
  transmute(source = "real",
            Y = stag,            # 观测时长 Y = min(T, C)
            delta = event)       # 事件指示

fake_A <- fake_A_post %>%
  transmute(source = "fake",
            Y = time,  # 统一成正的观测时长
            delta = event)
# 合并
dat_all_A <- bind_rows(real_A, fake_A)

#----------plot----------------
# ECDF（事件子样本）
pA1<- ggplot(filter(dat_all_A, delta == 1),
              aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x = "Y | delta = 1 (event times)",
       y = "ECDF")+
  theme_bw(14)
#save


#ECDF（删失子样本）
pA2<- ggplot(filter(dat_all_A, delta == 0),
              aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x="Y | delta = 0 (censored durations)",
       y = "ECDF") +
  theme_bw(14)
#save


pA3<- ggplot(fake_A_post, aes(x = time)) +
  geom_histogram(binwidth = 3, fill = "#FF0000"  , alpha = 0.7) +
  facet_wrap(~ event, labeller = labeller(event = c("0" = "Censored", "1" = "Event Occurred"))) +
  labs(
    x = "Months",
    y = "Count") +
  theme_minimal(base_size = 20)

ggsave("images/fake_duration_hist_Apost.pdf",plot =pA3,device = "pdf", 
       width = 6, height = 4)
ggsave("images/ecdf_event_Apost.pdf",plot =pA1,device = "pdf", 
       width = 6, height = 4)
ggsave("images/ecdf_censored_Apost.pdf",plot =pA2,device = "pdf", 
       width = 6, height = 4)


#---------------------------------------------
combo2b <- (p2b1 + p2b2 ) +
  plot_layout(guides = "collect", widths = c(1,1)) &
  theme(legend.position = "right")
combo2b
ggsave("images/ppc_two_a200.pdf",
       plot   = combo2b,
       width  = 7,   # 目标在论文中的实际宽度（英寸）
       height = 2.8,   # 合理高度，别太矮
       device = "pdf") 
