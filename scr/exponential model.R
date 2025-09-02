df <- read.csv("data/raw_data_turnover.csv")
table(df$event)
summary(df$stag)
length(df$event)
#Are there any missing values
sum(is.na(df))  
library(bayesplot)
library(rstan)
library(ggplot2)

ggplot(df, aes(x = stag)) +
  geom_histogram(binwidth = 3, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ event, labeller = labeller(event = c("0" = "Censored", "1" = "Event Occurred"))) +
  labs(x = "Months",
       y = "Count") +
  theme_minimal(base_size = 19)

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
  status = df$event          
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
##Take a look at the fit results
summary(fit_stan_turnover)$summary  
# —— traceplot —— 
bayesplot::mcmc_trace(
  as.array(fit_stan_turnover, pars = "lambda"),
  pars = "lambda"
)
ggsave("images/turnover_exp_model_traceplot.png", width = 6, height = 4)

#Take the lambda of posterior sampling
post_lam_turnover_stan <- rstan::extract(fit_stan_turnover,"lambda")$lambda

d <- sum(df$event == 1)
alpha<-0.001
sum_y  <- sum(df$stag)
beta<-0.001
# Use a grid to evaluate analytical density
lambda_grid  <- seq(min(post_lam_turnover_stan), max(post_lam_turnover_stan), length = 1000)
analytic_dens <- dgamma(lambda_grid, shape =  d + alpha, rate = sum_y+beta)

#Draw a lambda histogram of posterior sampling
png("images/turnover post lambda.png", width = 6, height = 4, units = "in", res = 300)
hist(post_lam_turnover_stan,
     breaks = 30,           
     prob   = TRUE,       
     col    = "#1f77b4",    
     border = "white",
     xlab   = expression(lambda),
     main   = "",
     cex.lab = 1.5,  
     cex.axis = 1.2,
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
       inset  = c(0.05, 0)
)
dev.off()



