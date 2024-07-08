library(tidyverse)
library(sf)
library(sp)
library(raster)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(spatstat)
library(posterior)
library(bridgesampling)

sourceCpp('source/cpp_help_func.cpp')
source('source/fit_mod_MCMC.R')
Y_obs <- read_csv('data/v10acs_pGC.csv')

ggplot(Y_obs, aes(x,y)) + geom_point(size= 0.1)

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
Y_obs[,c('x','y')] <- Y_obs[,c('x','y')]/4300*76

plot(S)
points(Y_obs[,c('x','y')])

Lim <- 25.75

c <- rbind(c(3278, 4073), c(3609, 631))/4300*76
colnames(c) <- c('x', 'y')
e <- c(1.4, 1)
theta <- c(pi/2, 0)

prior_np <- list(b0 = c(log(0.05), 0.4), 
                 N = data.frame(a = c(log(150), 0), b = c(0.25, 150)), 
                 R_eff = data.frame(a = c(log(2.5), 0.001), b = c(0.25, 15)),
                 n = data.frame(a = log(c(0.5, 1)), b = c(0.5, 0.75)),
                 mu = data.frame(a = rep(23, 3), b = rep(27, 3)),
                 sigma = data.frame(a = rep(0.5, 3), b = rep(1.9, 3)))

Theta_np <- list(b0 = 0.05, 
                 N = matrix(c(log(150), logit(10/150)), nrow = 1), 
                 R_eff = matrix(c(2.5, 5), nrow = 1), 
                 n = matrix(log(rep(1, 2)), nrow = 1),
                 mu = matrix(rep(26.2, 3), nrow = 1), 
                 sigma = matrix(rep(1, 3), nrow = 1))

prior_cp <- list(b0 = c(log(0.05), 0.4), 
                 N = data.frame(a = c(log(150), 0), b = c(0.25, 150)), 
                 R_eff = data.frame(a = c(log(2.5), 0.001), b = c(0.25, 15)),
                 n = data.frame(a = log(c(0.5, 1)), b = c(0.5, 0.75)),
                 mu = data.frame(a = 23, b = 27),
                 sigma = data.frame(a = 0.5, b = 1.9))

Theta_cp <- list(b0 = 0.05, 
                 N = matrix(c(log(150), logit(10/150)), nrow = 1), 
                 R_eff = matrix(c(2.5, 5), nrow = 1), 
                 n = matrix(log(rep(1, 2)), nrow = 1),
                 mu = 26, 
                 sigma = 1)

res_np <- readRDS('data/res_np_v10acs.RDS')
res_cp <- readRDS('data/res_cp_v10acs.RDS')

# log-posterior for no-pooling model
log_post_mod_np_gal_BF <- function(sample, data){
  Y_obs <- data$Y_obs
  c <- data$c
  S <- data$S
  grid <- data$grid
  Ng <- data$Ng
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
  return(log_post_mod_np_gal(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y, p, Lim, prior, Ng))
}

# log-posterior for complete-pooling model
log_post_mod_cp_gal_BF <- function(sample, data){
  Y_obs <- data$Y_obs
  c <- data$c
  S <- data$S
  grid <- data$grid
  Ng <- data$Ng
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
  return(log_post_mod_cp_gal(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y, p, Lim, prior, Ng))
}

data_np <- list(Y_obs = Y_obs, c = c, S = S, grid = grid, Ng = 1, e = e, theta = theta, p = 1, Lim = 25.75, prior = prior_np)
samp_np <- res_np
colnames(samp_np[[1]]) <- colnames(samp_np[[2]]) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11', 'V12', 'V13')
lb_np <- c(0, 0, 1e-6, 0, 0.001, rep(0,2), rep(23+1e-6, 3), rep(0.5+1e-6, 3))
ub_np <- c(Inf, Inf, 150, Inf, 15, rep(Inf, 2), rep(27-1e-6, 3), rep(1.9-1e-6, 3))
names(lb_np) <- names(ub_np) <- colnames(samp_np[[1]])

data_cp <- list(Y_obs = Y_obs, c = c, S = S, grid = grid, Ng = 1, e = e, theta = theta, p = 1, Lim = 25.75, prior = prior_cp)
samp_cp <- res_cp
colnames(samp_cp[[1]]) <- colnames(samp_cp[[2]]) <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9')

lb_cp <- c(0, 0, 1e-6, 0, 0.001, 0, 0, 23+1e-6, 0.5+1e-6)
ub_cp <- c(Inf, Inf, 150, Inf, 15, Inf, Inf, 27-1e-6, 1.9-1e-6)
names(lb_cp) <- names(ub_cp) <- colnames(samp_cp[[1]])

set.seed(378091)
# Marginal likelihood for no-pooling model
M_np <- bridge_sampler(samples = samp_np, 
                       log_posterior = log_post_mod_np_gal_BF,
                       data = data_np,
                       lb = lb_np, ub = ub_np, verbose = T, silent = F, 
                       cores = 4, method = 'warp3', repetitions = 30)
set.seed(378091)
# Marginal likelihood for complete-pooling model
M_cp <- bridge_sampler(samples = samp_cp, 
                       log_posterior = log_post_mod_cp_gal_BF,
                       data = data_cp,
                       lb = lb_cp, ub = ub_cp, verbose = T, silent = F, 
                       cores = 4, method = 'warp3', repetitions = 30)
 
bf(M_np, M_cp)

saveRDS(list(M_np, M_cp), 'data/v10acs_BF.RDS')

# posterior diagnostics

res_np <- readRDS('data/res_prob_v10acs.RDS')
res_cp <- readRDS('data/res_nprob_v10acs.RDS')

res_np_df1 <- as_draws_df(res_np[[1]])
res_np_df2 <- as_draws_df(res_np[[2]])
res_np_df <- bind_draws(res_np_df1, res_np_df2, along = 'chain')

res_cp_df1 <- as_draws_df(res_cp[[1]])
res_cp_df2 <- as_draws_df(res_cp[[2]])
res_cp_df <- bind_draws(res_cp_df1, res_cp_df2, along = 'chain')

print(summarise_draws(res_np_df), n = 23)
summarise_draws(res_cp_df)

mean(unlist(as.vector(res_cp_df[,3])))




