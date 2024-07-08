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

sourceCpp('source/cpp_help_func.cpp')
source('source/fit_mod_MCMC.R')
Y_obs <- read_csv('data/prob_GC_data/v7acs_pGC.csv')
Y_obs_J <- read_csv('data/prob_GC_data/V7ACS_pGC_Jans.csv')

Y_obs <- as.data.frame(Y_obs)
Y_obs[,c('x','y')] <- 76*Y_obs[,c('x','y')]/4300
Y_obs_J <- as.data.frame(Y_obs_J)
Y_obs_J[,c('x','y')] <- 76*Y_obs_J[,c('x','y')]/4300

c <- 76*rbind(c(661, 2987), c(1756, 752))/4300

R84_dat_harris <- filter(Y_obs, (x - c[1, 1])^2 + (y - c[1,2])^2 < 7.5^2)
R84_dat_harris$p <- rowMeans(R84_dat_harris[,8:507])
R84_dat_harris <- R84_dat_harris %>%
  filter(p > 0.05) %>%
  dplyr::select(x, y, M, C, p) %>%
  mutate(DATA = 'L24 ($p(\\mathrm{GC}) > 0.05$)')

sum(R84_dat_harris$p)

R84_dat_Jans <- filter(Y_obs_J, (x - c[1, 1])^2 + (y - c[1,2])^2 < 7.5^2)
R84_dat_Jans$p <- rowMeans(R84_dat_Jans[,8:507])
R84_dat_Jans <- R84_dat_Jans %>%
  filter(p > 0.05) %>%
  dplyr::select(x, y, M, C, p) %>%
  mutate(DATA = 'Prob J24 ($p(\\mathrm{GC}) > 0.05$)')


R84_dat_Jans_bin <- filter(Y_obs_J, (x - c[1, 1])^2 + (y - c[1,2])^2 < 7.5^2)
R84_dat_Jans_bin$p <- 1
R84_dat_Jans_bin <- R84_dat_Jans_bin %>%
  filter(p > 0 & M < 26.3 & C > 0.8 & C < 2.4) %>%
  dplyr::select(x, y, M, C, p) %>%
  mutate(DATA = 'Binary J24')

nrow(R84_dat_Jans_bin)

R84_dat <- bind_rows(R84_dat_harris, R84_dat_Jans, R84_dat_Jans_bin) %>%
  filter(C > 0.8)

ann_text <- data.frame(x = rep(15.5, 5), y = c(rep(c(59.5, 58.7),2), 59.1), 
                       lab = c("$N_{\\mathrm{GC}}: 19\\pm 9$ (Ours)",
                               '$\\mu_{\\mathrm{TO}}: 25.96_{-0.46}^{+0.43}$ mag (Ours)',
                               "$N_{\\mathrm{GC}}: 32\\pm 13$ (Ours)",
                               '$\\mu_{\\mathrm{TO}}: 26.12_{-0.39}^{+0.43}$ mag (Ours)',
                               '$N_{\\mathrm{GC}}: 43\\pm 6$ (J24)'),
                       DATA = c(rep(c('Prob J24 ($p(\\mathrm{GC}) > 0.05$)', 'L24 ($p(\\mathrm{GC}) > 0.05$)'), each = 2), 'Binary J24'))

R84_dat$DATA <- factor(R84_dat$DATA, levels = c("L24 ($p(\\mathrm{GC}) > 0.05$)", "Binary J24", "Prob J24 ($p(\\mathrm{GC}) > 0.05$)"))
ann_text$DATA <-  factor(ann_text$DATA, levels = c("L24 ($p(\\mathrm{GC}) > 0.05$)", "Binary J24", "Prob J24 ($p(\\mathrm{GC}) > 0.05$)"))

tikz(file = "R84_GC_F814W_plot.tex", standAlone=T, width = 8, height = 3)
ggplot(R84_dat, aes(x, y)) + geom_point(aes(color = M, size = p*1.5)) + facet_grid(~DATA) +
  geom_text(aes(x = x, y = y, label = lab), data = ann_text, size = 2.5) +
  scale_size_identity() + 
  xlim(c(11.68279 - 7.5, 11.68279 + 7.5)) + ylim(c(52.79349 - 7.5, 52.79349 + 7.5)) + 
  scale_color_viridis_c(direction = -1, name = 'F814W') + coord_fixed() +
  theme_minimal() + 
  theme(axis.title=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(),
        legend.justification = "left",
        legend.box.margin = margin(l = 0.1, unit = "cm")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
dev.off()
system('pdflatex R84_GC_F814W_plot.tex')

# tikz(file = "R84_harris_GC_plot.tex", standAlone=T, width = 4, height = 4)
# ggplot(R84_dat_harris, aes(x,y, color = M, size = p*2)) + geom_point() + 
#   scale_size_identity() +
#   xlim(c(11.68279 - 6, 11.68279 + 6)) + ylim(c(52.79349 - 6, 52.79349 + 6)) + 
#   scale_color_viridis_c(direction = -1, name = 'F814W') + 
#   theme_minimal() + theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) + 
#   xlab('X (kpc)') + ylab('Y (kpc)') + 
#   annotate('text', x = 14.5, y = 58.5, label = '$N_{\\mathrm{GC}}: 32\\pm 13$ (Ours)', size = 3) +
#   annotate('text', x = 14.5, y = 58, label = '$N_{\\mathrm{GC}}: 43\\pm 6$ (J24)', size = 3) +
#   annotate('text', x = 14.5, y = 57.5, label = '$\\mu_{\\mathrm{TO}}: 26.12_{-0.39}^{+0.43}$ mag (Ours)', size = 3)
# dev.off()
# system('pdflatex R84_harris_GC_plot.tex')
# 
# tikz(file = "R84_Jans_GC_plot.tex", standAlone=T, width = 4, height = 4)
# ggplot(R84_dat_Jans, aes(x,y, color = M, size = p*2)) + geom_point() + 
#   scale_size_identity() +
#   xlim(c(11.68279 - 6, 11.68279 + 6)) + ylim(c(52.79349 - 6, 52.79349 + 6)) + 
#   scale_color_viridis_c(direction = -1, name = 'F814W') +
#   theme_minimal() + theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) + 
#   xlab('X (kpc)') + ylab('Y (kpc)') + 
#   annotate('text', x = 14.5, y = 58.5, label = '$N_{\\mathrm{GC}}: 19\\pm 9$ (Ours)', size = 3) +
#   annotate('text', x = 14.5, y = 58, label = '$N_{\\mathrm{GC}}: 43\\pm 6$ (J24)', size = 3) +
#   annotate('text', x = 14.5, y = 57.5, label = '$\\mu_{\\mathrm{TO}}: 25.96_{-0.46}^{+0.43}$ mag (Ours)', size = 3)
# dev.off()
# system('pdflatex R84_Jans_GC_plot.tex')

# set up neceesary parameters and data for data analysis

# X <- c(0, 76, 76, 0)
# Y <- c(0, 0, 76, 76)
# S <- Polygon(cbind(X,Y))
# S <- SpatialPolygons(list(Polygons(list(S),'region')))
# S <- SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))
# 
# grid <- sp::makegrid(S, n = 100000)
# grid <- sp::SpatialPoints(grid, proj4string = CRS(proj4string(S)))
# grid <- raster::crop(grid, S)
# gridded(grid) <- T
# 
# dx1 <- unname(grid@grid@cellsize[1])
# dx2 <- unname(grid@grid@cellsize[2])
# 
# grid <- as.data.frame(grid)
# names(grid) <- c('x', 'y')
# grid <- as.matrix(grid)
# 
# grid <- list(grid, c(dx1, dx2))
# 
# Y_obs <- as.data.frame(Y_obs)
# Y_obs[,c('x','y')] <- 76*Y_obs[,c('x','y')]/4300
# 
# Lim <- 25.75
# 
# c <- 76*rbind(c(661, 2987), c(1756, 752))/4300
# e <- c(1,1)
# theta <- c(0,0)
# 
# prior_np <- list(b0 = c(log(0.035), 0.4), 
#                  N = data.frame(a = rep(0, 2), b = rep(150, 2)), 
#                  R_eff = data.frame(a = rep(0.001, 2), b = c(8, 5)),
#                  n = data.frame(a = log(rep(1, 2)), b = rep(0.75, 2)),
#                  mu = data.frame(a = rep(23, 3), b = rep(27, 3)),
#                  sigma = data.frame(a = rep(0.5, 3), b = rep(1.9, 3)))
# 
# Theta_np <- list(b0 = 0.03, 
#                  N = matrix(rep(logit(10/150), 2),nrow = 1), 
#                  R_eff = matrix(rep(1.5, 2), nrow = 1), 
#                  n = matrix(log(rep(1, 2)), nrow = 1),
#                  mu = matrix(rep(26.2, 3), nrow = 1), 
#                  sigma = matrix(rep(1, 3), nrow = 1))
# 
# prior_cp <- list(b0 = c(log(0.035), 0.4), 
#                  N = data.frame(a = rep(0, 2), b = rep(150, 2)), 
#                  R_eff = data.frame(a = rep(0.001, 2), b = c(8, 5)),
#                  n = data.frame(a = log(rep(1, 2)), b = rep(0.75, 2)),
#                  mu = data.frame(a = 23, b = 27),
#                  sigma = data.frame(a = 0.5, b = 1.9))
# 
# Theta_cp <- list(b0 = 0.03, 
#                  N = matrix(rep(logit(10/150), 2),nrow = 1), 
#                  R_eff = matrix(rep(1.5, 2), nrow = 1), 
#                  n = matrix(log(rep(1, 2)), nrow = 1),
#                  mu = 26, 
#                  sigma = 1)
# 
# # read in MCMC chains
# res_np <- readRDS('data/res_nprob_v7acs.RDS')
# res_cp <- readRDS('data/res_prob_v10acs.RDS')
# 
# # posterior diagnostics
# 
# res_np_df1 <- as_draws_df(res_np[[1]])
# res_np_df2 <- as_draws_df(res_np[[2]])
# res_np_df <- bind_draws(res_np_df1, res_np_df2, along = 'chain')
# 
# plot(acf(res_np[[1]][,2], lag = 200))
# 
# res_cp_df1 <- as_draws_df(res_nprob[[1]])
# res_cp_df2 <- as_draws_df(res_nprob[[2]])
# res_cp_df <- bind_draws(res_cp_df1, res_cp_df2, along = 'chain')
# 
# summarise_draws(res_np_df)
# summarise_draws(res_cp_df)
# 
# dat <- unlist(as.matrix(res_np_df))
# 
# cov(log(dat[,1:5]))
# 
# mlv(unlist(as.vector(res_cp_df[,3])))
# 
# # sensitivity check with power-scaling
# b_ratio <- function(b0, logmean = log(0.035), alpha){
#   return(dlnorm(b0, logmean, 0.4*alpha^(-0.5), log = T) - dlnorm(b0, log(0.035), 0.4, log = T))
# }
# 
# n_ratio <- function(n, logmean = log(1), alpha){
#   return(dlnorm(n, logmean, 0.75*alpha^(-0.5), log = T) - dlnorm(n, log(1), 0.75, log = T))
# }
# 
# #no-pooling model
# log_br_1_np <- psis(b_ratio(c(res_np[[1]][,1], res_np[[2]][,1]), alpha = 1.01), r_eff = 1737/75001)
# log_br_2_np <- psis(b_ratio(c(res_np[[1]][,1], res_np[[2]][,1]), alpha = 1/1.01), r_eff = 1737/75001)
# 
# log_n1r_1_np <- psis(n_ratio(c(res_np[[1]][,6], res_np[[2]][,6]), alpha = 1.01), r_eff = 1362/75001)
# log_n1r_2_np <- psis(n_ratio(c(res_np[[1]][,6], res_np[[2]][,6]), alpha = 1/1.01), r_eff = 1362/75001)
# 
# log_n2r_1_np <- psis(n_ratio(c(res_np[[1]][,7], res_np[[2]][,7]), alpha = 1.01), r_eff = 1493/75001)
# log_n2r_2_np <- psis(n_ratio(c(res_np[[1]][,7], res_np[[2]][,7]), alpha = 1/1.01), r_eff = 1493/75001)
# 
# log_r_1_np <- log_br_1_np$log_weights + log_n1r_1_np$log_weights + log_n2r_1_np$log_weights
# log_r_2_np <- log_br_2_np$log_weights + log_n1r_2_np$log_weights + log_n2r_2_np$log_weights
# 
# tf_dat_np <- rbind(res_np[[1]], res_np[[2]])
# 
# cjs_1_np <- apply(tf_dat_np, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 150002), y_weights = exp(log_r_1_np)))
# cjs_2_np <- apply(tf_dat_np, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 150002), y_weights = exp(log_r_2_np)))
# 
# #complete-pooling model
# log_br_1_cp <- psis(b_ratio(c(res_cp[[1]][,1], res_cp[[2]][,1]), alpha = 1.01), r_eff = 2026/75001)
# log_br_2_cp <- psis(b_ratio(c(res_cp[[1]][,1], res_cp[[2]][,1]), alpha = 1/1.01), r_eff = 2026/75001)
# 
# log_n1r_1_cp <- psis(n_ratio(c(res_cp[[1]][,6], res_cp[[2]][,6]), alpha = 1.01), r_eff = 1445/75001)
# log_n1r_2_cp <- psis(n_ratio(c(res_cp[[1]][,6], res_cp[[2]][,6]), alpha = 1/1.01), r_eff = 1445/75001)
# 
# log_n2r_1_cp <- psis(n_ratio(c(res_cp[[1]][,7], res_cp[[2]][,7]), alpha = 1.01), r_eff = 1329/75001)
# log_n2r_2_cp <- psis(n_ratio(c(res_cp[[1]][,7], res_cp[[2]][,7]), alpha = 1/1.01), r_eff = 1329/75001)
# 
# log_r_1_cp <- log_br_1_cp$log_weights + log_n1r_1_cp$log_weights + log_n2r_1_cp$log_weights
# log_r_2_cp <- log_br_2_cp$log_weights + log_n1r_2_cp$log_weights + log_n2r_2_cp$log_weights
# 
# tf_dat_cp <- rbind(res_cp[[1]], res_cp[[2]])
# 
# cjs_1_cp <- apply(tf_dat_cp, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 150002), y_weights = exp(log_r_1_cp)))
# cjs_2_cp <- apply(tf_dat_cp, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 150002), y_weights = exp(log_r_2_cp)))
# 
# # prior sensitivity
# ps_np <- (cjs_1_np + cjs_2_np)/2/log(1.01, base = 2)
# ps_cp <- (cjs_1_cp + cjs_2_cp)/2/log(1.01, base = 2)
# 
# log_lik_mod_np_PS <- function(sample, data, delta){
#   Y_obs <- data$Y_obs
#   c <- data$c
#   S <- data$S
#   grid <- data$grid
#   e <- data$e
#   theta <- data$theta
#   Lim <- data$Lim
#   p <- data$p
#   
#   K <- nrow(c)
#   b0 <- sample[1]
#   N <- sample[2:(K+1)]
#   R_eff <- sample[(K+2):(2*K+1)]
#   n <- sample[(2*K+2):(3*K+1)]
#   mu <- sample[(3*K+2):(4*K+2)]
#   sigma <- sample[(4*K+3):(5*K+3)]
#   
#   Y <- Y_obs[rbernoulli(nrow(Y_obs), p = Y_obs$GC_prob),]
#   return(delta*log_lik_mod_np(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y, p, Lim))
# }
# 
# # log-posterior for complete-pooling model
# log_lik_mod_cp_PS <- function(sample, data, delta){
#   Y_obs <- data$Y_obs
#   c <- data$c
#   S <- data$S
#   grid <- data$grid
#   e <- data$e
#   theta <- data$theta
#   Lim <- data$Lim
#   p <- data$p
#   
#   K <- nrow(c)
#   b0 <- sample[1]
#   N <- sample[2:(K+1)]
#   R_eff <- sample[(K+2):(2*K+1)]
#   n <- sample[(2*K+2):(3*K+1)]
#   mu <- sample[(3*K+2)]
#   sigma <- sample[(3*K+3)]
#   
#   Y <- Y_obs[rbernoulli(nrow(Y_obs), p = Y_obs$GC_prob),]
#   return(delta*log_lik_mod_cp(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y, p, Lim))
# }
# 
# # no-pooling likelihood sensitivity
# data_np_PS <- list(Y_obs = Y_obs, c = c, S = S, grid = grid, e = e, theta = theta, p = 1, Lim = 25.75)
# lik_ratio_np <- apply(tf_dat_np, 1, function(x){log_lik_mod_np_PS(x, data_np_PS, 0.01)})
# psis_lik_r_1_np <- psis(lik_ratio_np, r_eff = 771/75001)
# psis_lik_r_2_np <- psis(lik_ratio_np*100*(1/1.01 - 1), r_eff = 771/75001)
# 
# cjs_lik_1_np <- apply(tf_dat_np, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 150002), y_weights = exp(psis_lik_r_1_np$log_weights)))
# cjs_lik_2_np <- apply(tf_dat_np, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 150002), y_weights = exp(psis_lik_r_2_np$log_weights)))
# 
# lik_sens_np <- (cjs_lik_1_np + cjs_lik_2_np)/2/log(1.01, base = 2)
# 
# sense_table_np <- data.frame(Parameters = c('$\\beta_0$', '$N_1$', '$N_2$', '$R_\\mathrm{eff}^1$', '$R_\\mathrm{eff}^2$',
#                                          '$\\n_1$', '$\\n_2$', '$\\mu_0$', '$\\mu_1$', '$\\mu_2$', '$\\sigma_0$', '$\\sigma_1$', '$\\sigma_2$'),
#                           `Prior Sensitivity` = ps_np, `Likelihood Sensitivity` = lik_sens_np)
# 
# print(xtable(sense_table_np, digits=rep(3, 4), caption = "Prior and Likelihood Sensitivity under No-Pooling Model for V7-ACS Data", align = rep("c",4)), 
#       caption.placement = "top", include.rownames = FALSE, type = "latex", 
#       sanitize.text.function = function(x) {x})
# 
# # complete-pooling likelihood sensitivity
# data_cp_PS <- list(Y_obs = Y_obs, c = c, S = S, grid = grid, e = e, theta = theta, p = 1, Lim = 25.75)
# lik_ratio_cp <- apply(tf_dat_cp, 1, function(x){log_lik_mod_cp_PS(x, data_cp_PS, 0.01)})
# psis_lik_r_1_cp <- psis(lik_ratio_cp, r_eff = 800/75001)
# psis_lik_r_2_cp <- psis(lik_ratio_cp*100*(1/1.01 - 1), r_eff = 800/75001)
# 
# cjs_lik_1_cp <- apply(tf_dat_cp, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 150002), y_weights = exp(psis_lik_r_1_cp$log_weights)))
# cjs_lik_2_cp <- apply(tf_dat_cp, 2, function(x) cjs_dist(x, x, x_weights = rep(1, 150002), y_weights = exp(psis_lik_r_2_cp$log_weights)))
# 
# lik_sens_cp <- (cjs_lik_1_cp + cjs_lik_2_cp)/2/log(1.01, base = 2)
# 
# sense_table_cp <- data.frame(Parameters = c('$\\beta_0$', '$N_1$', '$N_2$', '$R_\\mathrm{eff}^1$', '$R_\\mathrm{eff}^2$',
#                                            '$n_1$', '$n_2$', '$\\mu_0$', '$\\sigma_0$'),
#                             `Prior Sensitivity` = ps_cp, `Likelihood Sensitivity` = lik_sens_cp)
# 
# print(xtable(sense_table_cp, digits=rep(3, 4), caption = "Prior and Likelihood Sensitivity under Complete-Pooling Model for V7-ACS Data", align = rep("c",4)), 
#       caption.placement = "top", include.rownames = FALSE, type = "latex", 
#       sanitize.text.function = function(x) {x})




