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
library(priorsense)
library(ggpubr)
library(wesanderson)

sourceCpp('source/cpp_help_func.cpp')
source('source/fit_mod_MCMC.R')
Y_obs <- read_csv('data/prob_GC_data/v11acs_pGC.csv')

Y_obs <- as.data.frame(Y_obs)
Y_obs[,c('x','y')] <- 76*Y_obs[,c('x','y')]/4300

c <- rbind(c(1125, 1780), c(2300, 5800), c(2628, 1532), c(2517, 3289))/4300*76

W88_harris_mag <- filter(Y_obs, (x - c[3, 1])^2 + (y - c[3,2])^2 < 6^2)
W88_harris_mag$p <- rowMeans(W88_harris_mag[,8:507])
W88_harris_mag <- W88_harris_mag %>%
  dplyr::select(x, y, M, C, p)


log_lik_mag <- function(par, data = W88_harris_mag) {
  mu <- par[1]
  sig <- par[2]
  
  M <- data$M
  p <- data$p
  
  log_lik <- -sum(p*(dnorm(M, mu, sig, log = T) + log(f_cpp(M, 25.75)) - log(Phi_f_cpp(25.75, mu, sig))))
  
  return(log_lik)
}

res <- optim(c(25, 1), log_lik_mag, method = "L-BFGS-B", lower = c(24, 0.5), upper = c(26.7, 1.9), hessian = T)

pnorm((26.3 - res$par[1])/sqrt(solve(res$hessian)[1,1]), 0, 1, lower.tail = F)


ggplot(W88_harris_mag, aes(C,M, color = p)) + geom_point() + scale_y_reverse() +scale_color_viridis_c()


W88_harris_mag_fix <- filter(Y_obs, (x - c[3, 1])^2 + (y - c[3,2])^2 < 6^2)
W88_harris_mag_fix <- W88_harris_mag_fix %>%
  dplyr::select(x, y, M, C) %>%
  filter(C > 1.0 & C < 2.4)


log_lik_mag_fix <- function(par, data = W88_harris_mag_fix) {
  mu <- par[1]
  sig <- par[2]
  
  M <- data$M
  p <- data$p
  
  log_lik <- -sum((dnorm(M, mu, sig, log = T) + log(f_cpp(M, 25.75)) - log(Phi_f_cpp(25.75, mu, sig))))
  
  return(log_lik)
}

res <- optim(c(25, 1), log_lik_mag_fix, method = "L-BFGS-B", lower = c(24, 0.5), upper = c(26.7, 1.9), hessian = T)

pnorm((26.3 - res$par[1])/sqrt(solve(res$hessian)[1,1]), 0, 1, lower.tail = F)




