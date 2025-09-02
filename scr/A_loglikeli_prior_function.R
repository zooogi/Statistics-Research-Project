## ----- Prepare data: You already have df -----
df <- read.csv("data/raw_data_turnover.csv")

y     <- df$stag
delta <- df$event
n     <- length(y)
ymax  <- max(y)



##########################似然函数！likelihood######################


## ----- Log likelihood logLik(lambda, A) -----
loglik_exp_uniform <- function(lambda, A, y, delta) {
  y      <- as.numeric(y)
  delta  <- as.integer(delta)
  ymax   <- max(y)
  
  # Basic constraint: parameter domain+indicator function support set
  if (lambda <= 0) return(-Inf)
  if (A < ymax)    return(-Inf)               #  0 <= y_i <= A
  
  # event: log[ lambda * exp(-lambda y) * (1 - y/A) ]
  # censored: log[ exp(-lambda y) / A ]
  log_event <- delta       * (log(lambda) - lambda * y + log1p(- y / A))  # log(1 - y/A)
  log_cens  <- (1 - delta) * (- lambda * y - log(A))
  sum(log_event + log_cens)
}





#############################Prior 先验分布




#Set prior distribution: lambda是gamma，A是uniform


## Weak information prior
## lambda ~ Gamma(shape, rate)
logprior_lambda_scalar <- function(l, shape = 0.001, rate = 0.001) {
  if (l <= 0) return(-Inf)
  dgamma(l, shape = shape, rate = rate, log = TRUE)  
}
## A ~ Uniform(a_min, a_max)，The lower bound must be ≥ ymax
logprior_A_scalar <- function(A, a_min = ymax + 1e-6, a_max = ymax + 500) {
  if (a_min >= a_max) stop("a_min must be < a_max")
  if (A < a_min || A > a_max) return(-Inf)
  -log(a_max - a_min)  # log density Uniform(a_min,a_max) 
}

##Package as a "vectorized version" using Vectorize
logprior_lambda <- Vectorize(logprior_lambda_scalar)
logprior_A <- Vectorize(function(A) logprior_A_scalar(A))


##############################################################