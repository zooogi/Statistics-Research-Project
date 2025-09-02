library(ggplot2)
library(tidyr)
library(ggnewscale)
library(dplyr)
source("scr/A_loglikeli_prior_function.R")

#-----------------Choose the grid range for the grid method!------------

#The range of lambda, using the numerical value corresponding to
#the quantile of the exact solution
alpha <- 0.001; beta <- 0.001
#formula 20（paper）
shape_post <- sum(delta) + alpha
rate_post  <- sum(y) + beta
lam_lo <- qgamma(0.001, shape_post, rate_post)  
lam_hi <- qgamma(0.999, shape_post, rate_post)  

## The range of A
A_lo <- ymax + 1e-6
A_hi <- ymax + 500

## grid size
nL <- 220
nA <- 220
lam_grid <- seq(lam_lo, lam_hi, length.out = nL)
A_grid   <- seq(A_lo,   A_hi,   length.out = nA)





## ==== ######Calculate log posterior ====


lp_mat <- matrix(NA_real_, nrow = nL, ncol = nA)
for (i in seq_len(nL)) {
  l <- lam_grid[i]
  for (j in seq_len(nA)) {
    A <- A_grid[j]
    lp <- loglik_exp_uniform(l, A, y, delta) +
      logprior_lambda_scalar(l, shape=alpha, rate=beta) +
      logprior_A_scalar(A, a_min=A_lo, a_max=A_hi)
    lp_mat[i, j] <- lp
  }
}





## ==== Normalized discrete posterior density
#-------（Numerical stability, log-sum-exp）====
dlam <- diff(lam_grid)[1]
dA   <- diff(A_grid)[1]
cell_area <- dlam * dA   # 等距网格：单元面积常数

lp_vec <- as.vector(lp_mat)

m      <- max(lp_vec)
post_unnorm <- exp(lp_vec - m)          # Avoid underflow
Z <- sum(post_unnorm) * cell_area    # Normalization constant for approximate integration
post_vec <- post_unnorm / Z             # Discretized 'density' (normalized with area)






#-----------Create a data frame for drawing--------------------


## Restore to Matrix&Data 
post_mat <- matrix(post_vec, nrow = nL, ncol = nA)
grid_df <- expand.grid(lambda = lam_grid, A = A_grid) %>%
  mutate(p = as.vector(post_mat))



## ==== Calculate 2D HPD threshold: 39.3% & 86.5% ====
hpd_thresholds <- function(pmat, probs = c(0.393, 0.865)){
  v <- sort(as.vector(pmat), decreasing = TRUE)
  cs <- cumsum(v) / sum(v)
  sapply(probs, function(q){
    idx <- which(cs >= q)[1]
    v[idx]
  })
}
levs <- hpd_thresholds(post_mat, probs = c(0.393, 0.865))
names(levs) <- c("HPD39.3", "HPD86.5")





## ====MAP  ====
which_max <- which.max(post_vec)
i_map <- ((which_max - 1) %% nL) + 1
j_map <- ((which_max - 1) %/% nL) + 1
map_point <- data.frame(lambda = lam_grid[i_map], A = A_grid[j_map])






## ====plot（ HPD  + MAP）==================

############# plot ##################

post_hpd <- ggplot(grid_df, aes(x = lambda, y = A)) +
  new_scale_fill() +
  # Fill HPD area
  geom_contour_filled(aes(z = p), breaks = c(sort(levs), Inf), alpha = 0.4) +
  ## Only draw contour lines
  geom_contour(aes(z = p), breaks = sort(levs), colour = "black", linewidth = 0.6) +
  # MAP 
  geom_point(data = map_point, aes(x = lambda, y = A),
             inherit.aes = FALSE, colour = "white", fill = "black",
             shape = 21, size = 2.2, stroke = 0.4) +
  annotate("text", x = map_point$lambda, y = map_point$A,
           label = "MAP", vjust = 1.8, size = 4.2) +
  geom_hline(yintercept = ymax, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(
    name   = "HPD bands",
    values = c("grey75", "grey30"),         
    labels = c("86.5% HPD", "39.3% HPD")
  ) +
  coord_cartesian(ylim = c(179, 182)) +
  labs(x = expression(lambda), y = "A")+
  theme_bw(base_size = 18)+
  theme(
    legend.position = "top",      
    legend.justification = "center"
  )

post_hpd
ggsave("images/post_contour.pdf",
       plot   = post_hpd,
       width  = 6,   
       height = 6,   
       device = "pdf") 

#——-------------——————————————————----——————————





#####################################################################


#----------------Comparison of lambda marginal a posteriori---------------

#The lambda posterior analytical solution of the model without A
d <- sum(df$event == 1)
alpha<-0.001
sum_y  <- sum(df$stag)
beta<-0.001
analytic_dens <- dgamma(lam_grid, shape =  d + alpha, rate = sum_y+beta)



## ==== Numerical margin: Sum along direction A ====
p_lambda_uni <- rowSums(post_mat) * dA                #Riemann summation
# Standardized inspection (should be close to 1)
sum_lambda <- sum(p_lambda_uni) * dlam
cat("Integral of p_lambda_uni ≈", sum_lambda, "\n")


#Transform into a data box for drawing
lambda_df <- tibble(
  lambda        = lam_grid,
  `Grid posterior (with A)` = p_lambda_uni,
  `Analytic posterior (no A)` = analytic_dens
) |>
  pivot_longer(cols = c(`Grid posterior (with A)`, `Analytic posterior (no A)`),
               names_to = "source", values_to = "density")


#plot
p_lambda_cmp <- ggplot(lambda_df, aes(x = lambda, y = density,
                                      colour = source, linetype = source)) +
  geom_line(linewidth = 0.9) +
  labs(x = expression(lambda), y = "posterior density") +
  geom_line(aes(linewidth = source)) +
  scale_linewidth_manual(values = c("Grid posterior (with A)" = 0.5,
                                    "Analytic posterior (no A)" = 1.5))+
  scale_color_manual(values = c("Grid posterior (with A)" = "#377EB8",
                                "Analytic posterior (no A)" = "#E41A1C")) +
  scale_linetype_manual(values = c("Grid posterior (with A)" = "solid",
                                   "Analytic posterior (no A)" = "dashed")) +
  theme_bw(base_size = 17) + theme(legend.title = element_blank())+
  theme(
    legend.position = "top",      
    legend.justification = "center"
  )
p_lambda_cmp

ggsave("images/lambda_marginal_compare.pdf",
       plot = p_lambda_cmp, width = 6, height = 4, device = "pdf")


#-----------------------------------------------



