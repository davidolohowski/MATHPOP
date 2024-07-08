library(tidyverse)
library(sf)
library(sp)
library(raster)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(posterior)
library(bridgesampling)
library(coda)

sourceCpp('source/cpp_help_func.cpp')
source('source/fit_mod_MCMC.R')
Y_obs <- read_csv('data/v9acs_pGC.csv')

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

c <- matrix(76*c(716, 3014)/4300, nrow=1)
e <- c(0.886)
theta <- c(0)

prior_np <- list(b0 = c(log(0.03), 0.4), 
                 N = data.frame(a = rep(0, 1), b = rep(150, 1)), 
                 R_eff = data.frame(a = rep(0.001, 1), b = c(8)),
                 n = data.frame(a = log(rep(1, 1)), b = rep(0.75, 1)),
                 mu = data.frame(a = rep(23, 2), b = rep(27, 2)),
                 sigma = data.frame(a = rep(0.5, 2), b = rep(1.9, 2)))

Theta_np <- list(b0 = 0.03, 
                 N = matrix(rep(logit(10/150), 1),nrow = 1), 
                 R_eff = matrix(rep(1.5, 1), nrow = 1), 
                 n = matrix(log(rep(1, 1)), nrow = 1),
                 mu = matrix(rep(26.2, 2), nrow = 1), 
                 sigma = matrix(rep(1, 2), nrow = 1))

prior_cp <- list(b0 = c(log(0.035), 0.4), 
                 N = data.frame(a = rep(0, 1), b = rep(150, 1)), 
                 R_eff = data.frame(a = rep(0.001, 1), b = c(8)),
                 n = data.frame(a = log(rep(1, 1)), b = rep(0.75, 1)),
                 mu = data.frame(a = 23, b = 27),
                 sigma = data.frame(a = 0.5, b = 1.9))

Theta_cp <- list(b0 = 0.03, 
                 N = matrix(rep(logit(10/150), 1),nrow = 1), 
                 R_eff = matrix(rep(1.5, 1), nrow = 1), 
                 n = matrix(log(rep(1, 1)), nrow = 1),
                 mu = 26, 
                 sigma = 1)

# read in MCMC chains
res_np <- readRDS('data/res_np_v9acs.RDS')
res_cp <- readRDS('data/res_cp_v9acs.RDS')

# log-posterior for no-pooling model
log_post_mod_np_BF <- function(sample, data){
  Y_obs <- data$Y_obs
  c <- data$c
  S <- data$S
  grid <- data$grid
  e <- data$e
  theta <- data$theta
  prior <- data$prior
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
  return(log_post_mod_np(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y, p, Lim, prior))
}

# log-posterior for complete-pooling model
log_post_mod_cp_BF <- function(sample, data){
  Y_obs <- data$Y_obs
  c <- data$c
  S <- data$S
  grid <- data$grid
  e <- data$e
  theta <- data$theta
  prior <- data$prior
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
  return(log_post_mod_cp(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y, p, Lim, prior))
}

data_np <- list(Y_obs = Y_obs, c = c, S = S, grid = grid, e = e, theta = theta, p = 1, Lim = 25.75, prior = prior_np)
samp_np <- res_np
colnames(samp_np[[1]]) <- colnames(samp_np[[2]]) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8')

lb_np <- c(0, 0, 0.001, 0, 23, 23, 0.5, 0.5)
ub_np <- c(Inf, 150, 8, Inf, 27, 27, 1.9, 1.9)
names(lb_np) <- names(ub_np) <- colnames(samp_np[[1]])

data_cp <- list(Y_obs = Y_obs, c = c, S = S, grid = grid, e = e, theta = theta, p = 1, Lim = 25.75, prior = prior_cp)
samp_cp <- res_cp
colnames(samp_cp[[1]]) <- colnames(samp_cp[[2]]) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6')

lb_cp <- c(0, 0, 0.001, 0, 23, 0.5)
ub_cp <- c(Inf, 150, 8, Inf, 27, 1.9)
names(lb_cp) <- names(ub_cp) <- colnames(samp_cp[[1]])

set.seed(378091)
# Marginal likelihood for no-pooling model
M_np <- bridge_sampler(samples = samp_np, 
                       log_posterior = log_post_mod_np_BF,
                       data = data_np,
                       lb = lb_np, ub = ub_np, verbose = T, 
                       cores = 4, repetitions = 30, method = 'warp3')
set.seed(378091)
# Marginal likelihood for complete-pooling model
M_cp <- bridge_sampler(samples = samp_cp, 
                       log_posterior = log_post_mod_cp_BF,
                       data = data_cp,
                       lb = lb_cp, ub = ub_cp, verbose = T, 
                       cores = 4, repetitions = 30, method = 'warp3')

saveRDS(list(M_np, M_cp), 'data/v9acs_BF.RDS')
# posterior diagnostics

res_np_df1 <- as_draws_df(res_np[[1]])
res_np_df2 <- as_draws_df(res_np[[2]])
res_np_df <- bind_draws(res_np_df1, res_np_df2, along = 'chain')

res_cp_df1 <- as_draws_df(res_cp[[1]])
res_cp_df2 <- as_draws_df(res_cp[[2]])
res_cp_df <- bind_draws(res_cp_df1, res_cp_df2, along = 'chain')

summarise_draws(res_np_df)
summarise_draws(res_cp_df)

