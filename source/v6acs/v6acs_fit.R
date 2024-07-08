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
Y_obs <- read_csv('/project/6055726/dli346/GC_count_PP/data/prob_GC_data/v6acs_pGC.csv')

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

c <- 76*rbind(c(900, 2905), c(1999, 2947), c(3385, 1726))/4300
e <- c(0.5, 0.65, 1)
theta <- c(0, 0, 0)
UDG_ID <- c('W13', 'W14', 'R116')

prior <- list(b0 = c(log(0.05), 0.4), 
                 N = data.frame(a = rep(0, 3), b = rep(50, 3)), 
                 R_eff = data.frame(a = rep(log(1.5), 3), b = rep(0.5, 3)),
                 n = data.frame(a = log(rep(1, 3)), b = rep(0.75, 3)),
                 mu = data.frame(a = rep(26.3, 4), b = rep(0.5, 4)),
                 sigma = data.frame(a = rep(log(1.3), 4), b = rep(0.25, 4)))

tune <- list(b0 = 0.001, N = 0.25, R = 0.2, n = 0.02, mu = 0.02, sigma = 0.02)

cf_error <- list(alpha = 1.5, m50 = 25.75, beta0 = 0.08836, beta1 = 0.645, m1 = 25.5)

set.seed(seed[sim_id])
Theta <- list(b0 = runif(1, 0.045, 0.055), 
              N = matrix(rnorm(3, log(10), 0.5), nrow = 1), 
              R_eff = matrix(rnorm(3, log(2), 0.5), nrow = 1), 
              n = matrix(rnorm(3, log(1), 0.75), nrow = 1),
              mu = matrix(rnorm(4, 26.2, 0.5), nrow = 1), 
              sigma = matrix(rnorm(4, log(1.2), 0.1), nrow = 1))

res_prob <- fit_mod_np(S, grid, Theta, UDG_ID, c, e, theta, Y_obs, p = 1, prior, M = 200000, tune, cf_error = cf_error)

saveRDS(res_prob, paste0('/project/6055726/dli346/GC_count_PP/data/v6acs/res_prob_v6acs_seed_', seed[sim_id], '.RDS'))


