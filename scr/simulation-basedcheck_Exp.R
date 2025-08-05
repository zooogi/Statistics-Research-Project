#------————————-run function------- 
source("src/generate_fake_data_exp.R")  

# -------------------------------------------------------------------------
fake_df <- generate_data_fixed_lambda(n = 1000, lambda = 0.05, a = -150, seed = 204)
head(fake_df)
table(fake_df$event)
write.csv(fake_df, "data/fake_data_sanity.csv", row.names = FALSE)

ggplot(fake_df, aes(x = time)) +
  geom_histogram(binwidth = 3, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ event, labeller = labeller(event = c("0" = "Censored", "1" = "Event Occurred"))) +
  labs(title = "Duration Histogram by Event Status",
       x = "Duration",
       y = "Count") +
  theme_minimal()


#---------------------------fit model--------------
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
stan_sanity_data <- list(
  N     = nrow(fake_df),
  y     = fake_df$time,       
  status = fake_df$event           # 0/1 向量
)

#fit
fit_exp_sanity <- stan(
  model_code = stan_code,
  data       = stan_sanity_data,
  chains     = 4,
  iter       = 1000,
  warmup     = 50,
  seed       = 204,
  refresh    = 0
)


#----------result------------------------
summary(fit_exp_sanity)$summary 

# —— traceplot —— 
bayesplot::mcmc_trace(
  as.array(fit_exp_sanity, pars = "lambda"),
  pars = "lambda"
)
ggsave("images/sanitycheck_exp_model_traceplot.png", width = 6, height = 4)

#取后验采样的lambda
post_sanity_exp <- extract(fit_exp_sanity,"lambda")$lambda


#画出后验采样的lambda直方图
hist(post_sanity_exp,
     breaks = 100,           
     prob   = TRUE,       
     col    = "#1f77b4",    
     border = "white",
     xlab   = expression(lambda),
     main   = "posterior of lambda",
     ylab   = "Density",
     ylim   = c(0, max(density(post_sanity_exp)$y)))
ggsave("images/sanitycheck_exp_post_lambda.png", width = 6, height = 4)
