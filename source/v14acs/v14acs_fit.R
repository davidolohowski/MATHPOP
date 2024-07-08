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
Y_obs <- read_csv('/project/6055726/dli346/GC_count_PP/data/prob_GC_data/v14acs_pGC.csv')

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

c <- rbind(c(3368, 2716), c(39, 725), c(93, 1403), c(2485, 1499), c(3754, 442))/4300*76
colnames(c) <- c('x', 'y')
e <- c(1.25, 2, 0.81, 0.88, 0.76)
theta <- c(-pi/18, -pi/3, pi/4, pi/4, 0)
UDG_ID <- c('W7', 'W6', 'W4')

prior <- list(b0 = c(log(0.05), 0.4), 
                 N = data.frame(a = c(log(300), log(40), 0, 0, 0), b = c(0.25, 0.25, rep(50, 3))), 
                 R_eff = data.frame(a = c(log(15), log(3), log(1.2), log(1), log(3.7)), b = c(0.25, 0.25, rep(0.5, 3))),
                 n = data.frame(a = log(c(rep(0.5, 2), rep(1, 3))), b = c(rep(0.5, 2), rep(0.75, 3))),
                 mu = data.frame(a = rep(26.5, 6), b = rep(0.5, 6)),
                 sigma = data.frame(a = rep(log(1.3), 6), b = rep(0.5, 6)))

tune <- list(b0 = 0.001, Ng = 0.25, Nu = 0.1, R = 0.1, n = 0.1, mu = 0.1, sigma = 0.1)

cf_error <- list(alpha = 1.5, m50 = 25.75, beta0 = 0.08836, beta1 = 0.645, m1 = 25.5)

set.seed(seed[sim_id])
Theta <- list(b0 = runif(1, 0.04, 0.06), 
              N = matrix(c(rnorm(1, log(300), 0.25), rnorm(1, log(40), 0.25), rnorm(3, log(10), 0.5)), nrow = 1), 
              R_eff = matrix(c(rnorm(1, log(15), 0.25), rnorm(1, log(3), 0.25), rnorm(3, log(2), 0.5)), nrow = 1), 
              n = matrix(c(rnorm(2, log(0.5), 0.5), rnorm(3, log(1), 0.75)), nrow = 1),
              mu = matrix(rnorm(6, 26.2, 0.5), nrow = 1), 
              sigma = matrix(rnorm(6, log(1.2), 0.1), nrow = 1))

set.seed(seed[sim_id])
res_prob <- fit_mod_np_gal(S, grid, Theta, Ng = 2, UDG_ID, c, e, theta, Y_obs, p = 1, prior, M = 400000, tune, cf_error = cf_error)

# 50min - 50k
saveRDS(res_prob, paste0('/project/6055726/dli346/GC_count_PP/data/v14acs/res_prob_v14acs_seed_', seed[sim_id], '.RDS'))

