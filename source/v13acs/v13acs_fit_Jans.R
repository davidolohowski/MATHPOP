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
Y_obs <- read_csv('/project/6055726/dli346/GC_count_PP/data/prob_GC_data/V13ACS_pGC_Jans.csv')

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

c <- rbind(c(3088, 4132), c(1570, 1392), c(3776, 294), c(1199, 2062))/4300*76
colnames(c) <- c('x', 'y')
e <- c(1.36, 0.74, 1.7, 0.64)
theta <- c(pi/4, 0, 0, -pi/4)
UDG_ID <- c('W22', 'W25', 'R117')

prior <- list(b0 = c(log(0.06), 0.4), 
                 N = data.frame(a = c(log(50), 0, 0, 0), b = c(0.25, rep(50, 3))), 
                 R_eff = data.frame(a = c(log(1.75), log(2.1), log(1.1), log(1)), b = c(0.25, rep(0.5, 3))),
                 n = data.frame(a = log(c(rep(0.5, 1), rep(1, 3))), b = c(rep(0.5, 1), rep(0.75, 3))),
                 mu = data.frame(a = rep(26.5, 5), b = rep(0.5, 5)),
                 sigma = data.frame(a = rep(log(1.3), 5), b = rep(0.25, 5)))

tune <- list(b0 = 0.01, Ng = 0.25, Nu = 0.1, R = 0.1, n = 0.1, mu = 0.1, sigma = 0.1)

cf_error <- list(alpha = 6.9, m50 = 26.8, beta0 = 0.07477, beta1 = 0.75094, m1 = 25.5)

set.seed(seed[sim_id])
Theta <- list(b0 = runif(1, 0.05, 0.07), 
              N = matrix(c(rnorm(1, log(50), 0.25), rnorm(3, log(10), 0.5)), nrow = 1), 
              R_eff = matrix(c(rnorm(1, log(1.75), 0.25), rnorm(3, log(2), 0.5)), nrow = 1), 
              n = matrix(c(rnorm(1, log(0.5), 0.5), rnorm(3, log(1), 0.75)), nrow = 1),
              mu = matrix(rnorm(5, 26.2, 0.5), nrow = 1), 
              sigma = matrix(rnorm(5, log(1.2), 0.1), nrow = 1))

set.seed(seed[sim_id])
res_prob <- fit_mod_np_gal(S, grid, Theta, Ng = 1, UDG_ID, c, e, theta, Y_obs, p = 1, prior, M = 300000, tune, cf_error = cf_error)

# 32min - 50k
saveRDS(res_prob, paste0('/project/6055726/dli346/GC_count_PP/data/v13acs/res_prob_v13acs_Jans_seed_', seed[sim_id], '.RDS'))

