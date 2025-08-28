#-----------run function----------
source("scr/gen_fake_data_weibull.R") 
#------------data----------------
df <- read.csv("data/raw_data_turnover.csv")

#-------------------package---------
library(dplyr)
library(ggplot2)


# ===== Fixed A consistent with the extended index model ----
A_fixed <- 179;
##### Extract posterior mean from the fit of rstan:: sampling
lambda_hat <- 0.01;
k_hat <-  1.06;
## Same sample size as real data
n <- nrow(df)
#generate fake data
df_fake_weibull <-generate_data_fixed_weibull(n = n ,
                           lambda =lambda_hat,   # rate
                           k      = k_hat,    # shape
                           a      = -A_fixed,   
                                seed   = 452)


#--------------Create a data frame drawing----------------
df_plot_real <- df %>%
  transmute(time = stag,
            status = ifelse(event == 1, "event times", "censored durations"),
            source = "real")

df_plot_fake <- df_fake_weibull %>%
  transmute(time = time,
            status = ifelse(event == 1, "event times", "censored durations"),
            source = "fake")

df_plot <- bind_rows(df_plot_real, df_plot_fake)

p <- ggplot(df_plot, aes(x = time, color = source, group = source)) +
  stat_ecdf(linewidth = 0.8) +
  facet_wrap(~ status, scales = "free_x") +
  labs(x = "Months", y = "ECDF") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right")
print(p)
ggsave("images/weibull_ecdf.png", plot = p, width = 7, height = 5, dpi = 300)
