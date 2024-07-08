library(tidyverse)
library(sf)
library(sp)
library(raster)
library(Rcpp)
library(RcppArmadillo)
library(spatstat)
library(VGAM)

args <- commandArgs(TRUE)
sim_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
seed <- c(7904, 5138, 6022)
sourceCpp('/project/6055726/dli346/GC_count_PP/source/cpp_help_func.cpp')
source('/project/6055726/dli346/GC_count_PP/source/fit_mod_MCMC.R')
Y_obs <- read_csv('/project/6055726/dli346/GC_count_PP/data/prob_GC_data/V12ACS_pGC_Jans.csv')

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

c <- rbind(c(3108, 3377), c(2514, -10), c(1305, 3353), c(686, 1629), c(2162, 1272), c(3692, 3717))/4300*76
colnames(c) <- c('x', 'y')
e <- c(1.14, 2, 1, 0.83, 1.3, 0.81)
theta <- c(pi/36, pi/6, 0, pi/12, 0, -pi/4)
UDG_ID <- c('W1','W2','W7')

prior <- list(b0 = c(log(0.035), 0.4), 
                 N = data.frame(a = c(log(50), log(100), log(20), 0, 0, 0), b = c(0.25, 0.25, 0.25, rep(50, 3))), 
                 R_eff = data.frame(a = c(log(1.5), log(4.5), log(1.5), log(1.8), log(3.7), log(1.2)), b = c(rep(0.25, 3), rep(0.5, 3))),
                 n = data.frame(a = log(c(rep(0.5, 3), rep(1, 3))), b = c(rep(0.5, 3), rep(0.75, 3))),
                 mu = data.frame(a = rep(26.3, 7), b = rep(0.5, 7)),
                 sigma = data.frame(a = rep(log(1.3), 7), b = rep(0.25, 7)))

tune <- list(b0 = 0.01, Ng = 0.25, Nu = 0.1, R = 0.1, n = 0.1, mu = 0.1, sigma = 0.1)

cf_error <- list(alpha = 6.7, m50 = 26.7, beta0 = 0.07477, beta1 = 0.75094, m1 = 25.5)

set.seed(seed[sim_id])
Theta <- list(b0 = runif(1, 0.025, 0.045), 
              N = matrix(c(rnorm(1, log(50), 0.25), rnorm(1, log(100), 0.25), rnorm(1, log(20), 0.5), rnorm(3, log(10), 0.5)), nrow = 1), 
              R_eff = matrix(c(rnorm(1, log(1.5), 0.25), rnorm(1, log(4.5), 0.25), rnorm(1, log(1.5), 0.25), rnorm(3, log(2), 0.5)), nrow = 1), 
              n = matrix(c(rnorm(3, log(0.5), 0.5), rnorm(3, log(1), 0.75)), nrow = 1),
              mu = matrix(rnorm(7, 26.2, 0.5), nrow = 1), 
              sigma = matrix(rnorm(7, log(1.2), 0.1), nrow = 1))

set.seed(seed[sim_id])
res_prob <- fit_mod_np_gal(S, grid, Theta, Ng = 3, UDG_ID, c, e, theta, Y_obs, p = 1, prior, M = 400000, tune, cf_error = cf_error)

# 56min - 50k
saveRDS(res_prob, paste0('/project/6055726/dli346/GC_count_PP/data/v12acs/res_prob_v12acs_Jans_seed_', seed[sim_id], '.RDS'))


