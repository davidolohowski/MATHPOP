library(tidyverse)
library(sf)
library(sp)
library(raster)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(posterior)
library(coda)
library(loo)
library(xtable)

sourceCpp('source/cpp_help_func.cpp')
source('source/fit_mod_MCMC.R')
Y_obs <- read_csv('data/v6acs_pGC.csv')

# set up neceesary parameters and data for data analysis

X <- c(0, 76, 76, 0)
Y <- c(0, 0, 76, 76)
S <- Polygon(cbind(X,Y))
S <- SpatialPolygons(list(Polygons(list(S),'region')))
S <- SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))

grid <- sp::makegrid(S, n = 100000)
grid <- sp::SpatialPoints(grid, proj4string = CRS(proj4string(S)))
grid <- raster::crop(grid, S)
gridded(grid) <- T

dx1 <- unname(grid@grid@cellsize[1])
dx2 <- unname(grid@grid@cellsize[2])

grid <- as.data.frame(grid)
names(grid) <- c('x', 'y')
grid <- as.matrix(grid)

grid <- list(grid, c(dx1, dx2))

Y_obs <- as.data.frame(Y_obs)
Y_obs[,c('x','y')] <- 76*Y_obs[,c('x','y')]/4300

Lim <- 25.75

c <- 76*rbind(c(909, 2904), c(1999, 2947), c(3385, 1726))/4300
e <- c(0.5, 0.65, 1)
theta <- c(0, 0, 0)

# read in MCMC chains
res_np <- readRDS('data/res_np_v6acs.RDS')
res_cp <- readRDS('data/res_cp_v6acs.RDS')

# posterior diagnostics

res_np_df1 <- as_draws_df(res_np[[1]])
res_np_df2 <- as_draws_df(res_np[[2]])
res_np_df <- bind_draws(res_np_df1, res_np_df2, along = 'chain')

res_cp_df1 <- as_draws_df(res_cp[[1]])
res_cp_df2 <- as_draws_df(res_cp[[2]])
res_cp_df <- bind_draws(res_cp_df1, res_cp_df2, along = 'chain')

summarise_draws(res_np_df)
summarise_draws(res_cp_df)

# sensitivity check with power-scaling
b_ratio <- function(b0, logmean = log(0.05), alpha){
  return(dlnorm(b0, logmean, 0.4*alpha^(-0.5), log = T) - dlnorm(b0, log(0.05), 0.4, log = T))
}

n_ratio <- function(n, logmean = log(1), alpha){
  return(dlnorm(n, logmean, 0.75*alpha^(-0.5), log = T) - dlnorm(n, log(1), 0.75, log = T))
}

#no-pooling model
log_br_1_np <- psis(b_ratio(c(res_np[[1]][,1], res_np[[2]][,1]), alpha = 1.01), r_eff = 1080/170001)
log_br_2_np <- psis(b_ratio(c(res_np[[1]][,1], res_np[[2]][,1]), alpha = 1/1.01), r_eff = 1080/170001)

log_n1r_1_np <- psis(n_ratio(c(res_np[[1]][,8], res_np[[2]][,8]), alpha = 1.01), r_eff = 884/170001)
log_n1r_2_np <- psis(n_ratio(c(res_np[[1]][,8], res_np[[2]][,8]), alpha = 1/1.01), r_eff = 884/170001)

log_n2r_1_np <- psis(n_ratio(c(res_np[[1]][,9], res_np[[2]][,9]), alpha = 1.01), r_eff = 947/170001)
log_n2r_2_np <- psis(n_ratio(c(res_np[[1]][,9], res_np[[2]][,9]), alpha = 1/1.01), r_eff = 947/170001)

log_n3r_1_np <- psis(n_ratio(c(res_np[[1]][,10], res_np[[2]][,10]), alpha = 1.01), r_eff = 962/170001)
log_n3r_2_np <- psis(n_ratio(c(res_np[[1]][,10], res_np[[2]][,10]), alpha = 1/1.01), r_eff = 962/170001)

log_r_1_np <- log_br_1_np$log_weights + log_n1r_1_np$log_weights + log_n2r_1_np$log_weights + log_n3r_1_np$log_weights
log_r_2_np <- log_br_2_np$log_weights + log_n1r_2_np$log_weights + log_n2r_2_np$log_weights + log_n3r_1_np$log_weights

dat_np <- rbind(res_np[[1]], res_np[[2]])
  
cjs_1_np <- apply(dat_np, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 340002), y_weights = exp(log_r_1_np)))
cjs_2_np <- apply(dat_np, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 340002), y_weights = exp(log_r_2_np)))

#complete-pooling model
log_br_1_cp <- psis(b_ratio(c(res_cp[[1]][,1], res_cp[[2]][,1]), alpha = 1.01), r_eff = 1341/170001)
log_br_2_cp <- psis(b_ratio(c(res_cp[[1]][,1], res_cp[[2]][,1]), alpha = 1/1.01), r_eff = 1341/170001)

log_n1r_1_cp <- psis(n_ratio(c(res_cp[[1]][,8], res_cp[[2]][,8]), alpha = 1.01), r_eff = 1341/170001)
log_n1r_2_cp <- psis(n_ratio(c(res_cp[[1]][,8], res_cp[[2]][,8]), alpha = 1/1.01), r_eff = 1341/170001)

log_n2r_1_cp <- psis(n_ratio(c(res_cp[[1]][,9], res_cp[[2]][,9]), alpha = 1.01), r_eff = 1422/170001)
log_n2r_2_cp <- psis(n_ratio(c(res_cp[[1]][,9], res_cp[[2]][,9]), alpha = 1/1.01), r_eff = 1422/170001)

log_n3r_1_cp <- psis(n_ratio(c(res_cp[[1]][,10], res_cp[[2]][,10]), alpha = 1.01), r_eff = 1188/170001)
log_n3r_2_cp <- psis(n_ratio(c(res_cp[[1]][,10], res_cp[[2]][,10]), alpha = 1/1.01), r_eff = 1188/170001)

log_r_1_cp <- log_br_1_cp$log_weights + log_n1r_1_cp$log_weights + log_n2r_1_cp$log_weights + log_n3r_1_cp$log_weights
log_r_2_cp <- log_br_2_cp$log_weights + log_n1r_2_cp$log_weights + log_n2r_2_cp$log_weights + log_n3r_2_cp$log_weights

dat_cp <- rbind(res_cp[[1]], res_cp[[2]])

cjs_1_cp <- apply(dat_cp, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 340002), y_weights = exp(log_r_1_cp)))
cjs_2_cp <- apply(dat_cp, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 340002), y_weights = exp(log_r_2_cp)))

# prior sensitivity
ps_np <- (cjs_1_np + cjs_2_np)/2/log(1.01, base = 2)
ps_cp <- (cjs_1_cp + cjs_2_cp)/2/log(1.01, base = 2)

log_lik_mod_np_PS <- function(sample, data, delta){
  Y_obs <- data$Y_obs
  c <- data$c
  S <- data$S
  grid <- data$grid
  e <- data$e
  theta <- data$theta
  Lim <- data$Lim
  p <- data$p
  
  K <- nrow(c)
  b0 <- sample[1]
  N <- sample[2:(K+1)]
  R_eff <- sample[(K+2):(2*K+1)]
  n <- sample[(2*K+2):(3*K+1)]
  mu <- sample[(3*K+2):(4*K+2)]
  sigma <- sample[(4*K+3):(5*K+3)]
  
  Y <- Y_obs[rbernoulli(nrow(Y_obs), p = Y_obs$GC_prob),]
  return(delta*log_lik_mod_np(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y, p, Lim))
}

# log-posterior for complete-pooling model
log_lik_mod_cp_PS <- function(sample, data, delta){
  Y_obs <- data$Y_obs
  c <- data$c
  S <- data$S
  grid <- data$grid
  e <- data$e
  theta <- data$theta
  Lim <- data$Lim
  p <- data$p
  
  K <- nrow(c)
  b0 <- sample[1]
  N <- sample[2:(K+1)]
  R_eff <- sample[(K+2):(2*K+1)]
  n <- sample[(2*K+2):(3*K+1)]
  mu <- sample[(3*K+2)]
  sigma <- sample[(3*K+3)]
  
  Y <- Y_obs[rbernoulli(nrow(Y_obs), p = Y_obs$GC_prob),]
  return(delta*log_lik_mod_cp(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y, p, Lim))
}

# no-pooling likelihood sensitivity
data_np_PS <- list(Y_obs = Y_obs, c = c, S = S, grid = grid, e = e, theta = theta, p = 1, Lim = 25.75)
lik_ratio_np <- apply(dat_np, 1, function(x){log_lik_mod_np_PS(x, data_np_PS, 0.01)})
psis_lik_r_1_np <- psis(lik_ratio_np, r_eff = 562/170001)
psis_lik_r_2_np <- psis(lik_ratio_np*100*(1/1.01 - 1), r_eff = 562/170001)

cjs_lik_1_np <- apply(dat_np, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 340002), y_weights = exp(psis_lik_r_1_np$log_weights)))
cjs_lik_2_np <- apply(dat_np, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 340002), y_weights = exp(psis_lik_r_2_np$log_weights)))

lik_sens_np <- (cjs_lik_1_np + cjs_lik_2_np)/2/log(1.01, base = 2)

sense_table_np <- data.frame(Parameters = c('$\\beta_0$', '$N_1$', '$N_2$', '$N_3$', '$R_\\mathrm{eff}^1$', '$R_\\mathrm{eff}^2$', '$R_\\mathrm{eff}^3$',
                                            '$n_1$', '$n_2$', '$n_3$', '$\\mu_0$', '$\\mu_1$', '$\\mu_2$', '$\\mu_3$', '$\\sigma_0$', '$\\sigma_1$', '$\\sigma_2$', '$\\sigma_3$'),
                             `Prior Sensitivity` = ps_np, `Likelihood Sensitivity` = lik_sens_np)

print(xtable(sense_table_np, digits=rep(3, 4), caption = "Prior and Likelihood Sensitivity under No-Pooling Model for V6-ACS Data", align = rep("c",4)), 
      caption.placement = "top", include.rownames = FALSE, type = "latex", 
      sanitize.text.function = function(x) {x})

# complete-pooling likelihood sensitivity
data_cp_PS <- list(Y_obs = Y_obs, c = c, S = S, grid = grid, e = e, theta = theta, p = 1, Lim = 25.75)
lik_ratio_cp <- apply(dat_cp, 1, function(x){log_lik_mod_cp_PS(x, data_cp_PS, 0.01)})
psis_lik_r_1_cp <- psis(lik_ratio_cp, r_eff = 907/170001)
psis_lik_r_2_cp <- psis(lik_ratio_cp*100*(1/1.01 - 1), r_eff = 907/170001)

cjs_lik_1_cp <- apply(dat_cp, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 340002), y_weights = exp(psis_lik_r_1_cp$log_weights)))
cjs_lik_2_cp <- apply(dat_cp, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 340002), y_weights = exp(psis_lik_r_2_cp$log_weights)))

median(dat_cp[,2])
x <- sample(dat_cp[,2], size = 340002, prob = exp(psis_lik_r_1_cp$log_weights), replace = T)
median(x)

lik_sens_cp <- (cjs_lik_1_cp + cjs_lik_2_cp)/2/log(1.01, base = 2)

sense_table_cp <- data.frame(Parameters = c('$\\beta_0$', '$N_1$', '$N_2$', '$N_3$', '$R_\\mathrm{eff}^1$', '$R_\\mathrm{eff}^2$', '$R_\\mathrm{eff}^3$',
                                            '$n_1$', '$n_2$', '$n_3$', '$\\mu_0$', '$\\sigma_0$'),
                             `Prior Sensitivity` = ps_cp, `Likelihood Sensitivity` = lik_sens_cp)

print(xtable(sense_table_cp, digits=rep(3, 4), caption = "Prior and Likelihood Sensitivity under Complete-Pooling Model for V6-ACS Data", align = rep("c",4)), 
      caption.placement = "top", include.rownames = FALSE, type = "latex", 
      sanitize.text.function = function(x) {x})




