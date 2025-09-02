source("scr/A_post_contour&marginal.R")
source("scr/generate_fake_data_exp.R") 
library(ggplot2)
library(patchwork)
N<- length(df$event)
lambda_map <- map_point[[1]] 
A_map <- map_point[[2]]
fake_A_post<- generate_data_fixed_lambda(n = N, lambda = lambda_map,
                                                     a = -A_map, seed = 25)
#----------------summary------------------
summary(fake_A_post)
table(fake_A_post$event)

#----------------Merge into data frame
real_A <- df %>%
  transmute(source = "real",
            Y = stag,            # Observation duration Y = min(T, C)
            delta = event)    

fake_A <- fake_A_post %>%
  transmute(source = "fake",
            Y = time,  
            delta = event)
# merge
dat_all_A <- bind_rows(real_A, fake_A)

#----------plot----------------
# ECDF (Event subsample)
pA1<- ggplot(filter(dat_all_A, delta == 1),
              aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x = "Y | delta = 1 (event times)",
       y = "ECDF")+
  theme_bw(14)

pA1

#ECDF (censored subsample)
pA2<- ggplot(filter(dat_all_A, delta == 0),
              aes(x = Y, colour = source)) +
  stat_ecdf(size = 1) +
  labs(x="Y | delta = 0 (censored durations)",
       y = "ECDF") +
  theme_bw(14)
pA2


pA3<- ggplot(fake_A_post, aes(x = time)) +
  geom_histogram(binwidth = 3, fill = "#FF0000"  , alpha = 0.7) +
  facet_wrap(~ event, labeller = labeller(event = c("0" = "Censored", "1" = "Event Occurred"))) +
  labs(
    x = "Months",
    y = "Count") +
  theme_minimal(base_size = 20)
pA3
ggsave("images/fake_duration_hist_Apost.pdf",plot =pA3,device = "pdf", 
       width = 6, height = 4)
ggsave("images/ecdf_event_Apost.pdf",plot =pA1,device = "pdf", 
       width = 6, height = 4)
ggsave("images/ecdf_censored_Apost.pdf",plot =pA2,device = "pdf", 
       width = 6, height = 4)


#---------------------------------------------
comb_map <- (pA1 + pA2 ) +
  plot_layout(guides = "collect", widths = c(1,1)) &
  theme(legend.position = "right")
comb_map
ggsave("images/ppc_two_map.pdf",
       plot   = comb_map,
       width  = 7,   
       height = 2.8,   
       device = "pdf") 
