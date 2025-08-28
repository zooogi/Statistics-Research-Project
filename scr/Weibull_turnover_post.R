library(rstan)
library(bayesplot)
#data
df <- read.csv("data/raw_data_turnover.csv")



# =======  Stan model code (rate-parameterized Weibull with right censoring)  =======
stan_weibull <- "
data {
  int<lower=0> Nevent;
  vector<lower=0>[Nevent] y_event;
  int<lower=0> Ncens;
  vector<lower=0>[Ncens] y_cens;
}
transformed data {
  int N = Nevent + Ncens;   // 这里定义的 N 属于 transformed data
}
parameters {
  real log_lambda;    // log(rate)
  real log_k;         // log(shape)
}
transformed parameters {
  real<lower=0> lambda = exp(log_lambda);
  real<lower=0> k      = exp(log_k);
  // Stan uses (shape=k, scale=theta); convert rate->scale by theta = lambda^{-1/k}
  real<lower=0> theta  = exp(-log_lambda / k);
}
model {
  // weakly-informative priors 
  log_lambda ~ normal(0, 2);
  log_k      ~ normal(0, 1.5);

  // likelihood: events use lpdf, censorings use survival (lccdf)
  if (Nevent > 0) target += weibull_lpdf(y_event | k, theta);
  if (Ncens  > 0) target += weibull_lccdf(y_cens  | k, theta);
}
generated quantities {
  // for diagnostics / model comparison (LOO/WAIC)
  vector[N] log_lik;
  {
    int idx = 1;
    for (i in 1:Nevent) {
      log_lik[idx] = weibull_lpdf(y_event[i] | k, theta);
      idx += 1;
    }
    for (i in 1:Ncens) {
      log_lik[idx] = weibull_lccdf(y_cens[i] | k, theta);
      idx += 1;
    }
  }
}
"

sm_w <- stan_model(model_code = stan_weibull)


####
y_event <- df$stag[df$event == 1]
y_cens  <- df$stag[df$event == 0]

stan_weibull <- list(
  Nevent = length(y_event),
  y_event = if (length(y_event)>0) as.vector(y_event) else numeric(0),
  Ncens  = length(y_cens),
  y_cens = if (length(y_cens)>0)  as.vector(y_cens)  else numeric(0)
)

# ---- sample ----
fit_wei <- sampling(
  sm_w , data = stan_weibull,
  chains = 4, iter = 2000, warmup = 1000, seed = 27,refresh = 0,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

#summary
print(fit_wei, pars = c("log_lambda","log_k","lambda","k"), probs = c(0.025,0.5,0.975))
quantile(post$k, c(.025,.5,.975))
mean(post$k > 1) 
##
pdf("images/weibull_traceplot.pdf", width = 7, height = 5)
traceplot(fit_wei, pars = c("lambda", "k"))
dev.off()

# Compare with different chains, posterior density plot
posterior <- as.array(fit_wei)
p <- mcmc_dens_overlay(posterior, pars = c("lambda", "k"))
ggsave("images/weibull_posterior_density.png", plot = p, width = 7, height = 5, dpi = 300)

#Posterior histogram
post <- as.data.frame(rstan::extract(fit_wei, pars = c("lambda", "k")))
# Posterior histogram: lambda
png("images/weibull_hist_lambda.png", width = 800, height = 600)
hist(post$lambda, breaks = 50, col = "skyblue",
     main = "Posterior of lambda (Weibull)",
     xlab = "lambda")
dev.off()

# Posterior histogram: k
png("images/weibull_hist_k.png", width = 800, height = 600)
hist(post$k, breaks = 50, col = "salmon",
     main = "Posterior of k (Weibull)",
     xlab = "k")
dev.off()

