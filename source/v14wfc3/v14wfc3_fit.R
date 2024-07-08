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
Y_obs <- read_csv('/project/6055726/dli346/GC_count_PP/data/prob_GC_data/v14wfc3_pGC.csv')

X <- c(0, 62, 62, 0)
Y <- c(0, 0, 62, 62)
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
Y_obs[,c('x','y')] <- 62*Y_obs[,c('x','y')]/4400

c <- rbind(c(1979, 2957), c(2185, 983), c(358, 390))/4400*62
colnames(c) <- c('x', 'y')
e <- c(2, 1/0.52, 0.83)
theta <- c(4*pi/9, 0, 0)
UDG_ID <- c('W18', 'W19')

prior <- list(b0 = c(log(0.03), 0.4), 
                 N = data.frame(a = c(log(40), 0, 0), b = c(0.25, rep(50, 2))), 
                 R_eff = data.frame(a = c(log(2), log(3.1), log(1.84)), b = c(0.25, rep(0.5, 2))),
                 n = data.frame(a = log(c(rep(0.5, 1), rep(1, 2))), b = c(rep(0.5, 1), rep(0.75, 2))),
                 mu = data.frame(a = rep(26.3, 4), b = rep(0.5, 4)),
                 sigma = data.frame(a = rep(log(1.3), 4), b = rep(0.25, 4)))

tune <- list(b0 = 0.01, Ng = 0.25, Nu = 0.1, R = 0.1, n = 0.1, mu = 0.1, sigma = 0.1)

cf_error <- list(alpha = 1.57, m50 = 26.52, beta0 = 0.0977, beta1 = 0.613, m1 = 26)

set.seed(seed[sim_id])
Theta <- list(b0 = runif(1, 0.02, 0.05), 
              N = matrix(c(rnorm(1, log(40), 0.25), rnorm(2, log(10), 0.5)), nrow = 1), 
              R_eff = matrix(rnorm(3, log(2), 0.5), nrow = 1), 
              n = matrix(c(rnorm(1, log(0.5), 0.5), rnorm(2, log(1), 0.75)), nrow = 1),
              mu = matrix(rnorm(4, 26.2, .5), nrow = 1), 
              sigma = matrix(rnorm(4, log(1.2), 0.1), nrow = 1))

set.seed(seed[sim_id])
res_prob <- fit_mod_np_gal(S, grid, Theta, Ng = 1, UDG_ID, c, e, theta, Y_obs, p = 1, prior, M = 200000, tune, cf_error = cf_error)

# 31min - 50k
saveRDS(res_prob, paste0('/project/6055726/dli346/GC_count_PP/data/v14wfc3/res_prob_v14wfc3_seed_', seed[sim_id], '.RDS'))


