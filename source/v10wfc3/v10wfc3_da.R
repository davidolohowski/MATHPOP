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
library(reshape2)

sourceCpp('source/cpp_help_func.cpp')
source('source/fit_mod_MCMC.R')
Y_obs <- read_csv('data/prob_GC_data/v10wfc3_pGC.csv')
Y_obs_J <- read_csv('data/prob_GC_data/V10WFC3_pGC_Jans.csv')

Y_obs <- as.data.frame(Y_obs)
Y_obs[,c('x','y')] <- 62*Y_obs[,c('x','y')]/4400
Y_obs_J <- as.data.frame(Y_obs_J)
Y_obs_J[,c('x','y')] <- 62*Y_obs_J[,c('x','y')]/4400

c <- 62*rbind(c(808, 744), c(1930, 2653), c(2695, 2132))/4400

R27_dat_harris <- filter(Y_obs, (x - c[1, 1])^2 + (y - c[1,2])^2 < 7.5^2)
R27_dat_harris$p <- rowMeans(R27_dat_harris[,8:507])
R27_dat_harris <- R27_dat_harris %>%
  filter(p > 0.05) %>%
  dplyr::select(x, y, M, C, p) %>%
  mutate(DATA = 'L24 ($p(\\mathrm{GC}) > 0.05$)')

R27_dat_Jans <- filter(Y_obs_J, (x - c[1, 1])^2 + (y - c[1,2])^2 < 7.5^2)
R27_dat_Jans$p <- rowMeans(R27_dat_Jans[,8:507])
R27_dat_Jans <- R27_dat_Jans %>%
  filter(p > 0.05) %>%
  dplyr::select(x, y, M, C, p) %>%
  mutate(DATA = 'Prob J24 ($p(\\mathrm{GC}) > 0.05$)')

sum(R27_dat_Jans$p)

R27_dat_Jans_bin <- filter(Y_obs_J, (x - c[1, 1])^2 + (y - c[1,2])^2 < 7.5^2)
R27_dat_Jans_bin$p <- 1
R27_dat_Jans_bin <- R27_dat_Jans_bin %>%
  filter(p > 0 & M < 26.3 & C > 0.8 & C < 2.4) %>%
  dplyr::select(x, y, M, C, p) %>%
  mutate(DATA = 'Binary J24')

mean(R27_dat_Jans_bin$M)

R27_dat <- bind_rows(R27_dat_harris, R27_dat_Jans, R27_dat_Jans_bin)

ann_text_R27 <- data.frame(x = c(9, 9, 16.5, 9, 9, 16, 8.5, 16.5), y = c(17, 16, 5.5, 17, 16, 5.5, 16.5, 5.5), 
                           lab = c("$N_{\\mathrm{GC}}: 26_{-9}^{+11}$ (Ours)",
                                   '$\\mu_{\\mathrm{TO}}: 25.35_{-0.28}^{+0.38}$ mag (Ours)',
                                   '$N_{\\mathrm{cand}} = 26.3$',
                                   "$N_{\\mathrm{GC}}: 37_{-13}^{+10}$ (Ours)",
                                   '$\\mu_{\\mathrm{TO}}: 25.75_{-0.25}^{+0.32}$ mag (Ours)',
                                   '$N_{\\mathrm{cand}} = 27.6$',
                                   '$N_{\\mathrm{GC}}: 52\\pm 8$ (J24)',
                                   '$N_{\\mathrm{cand}} = 37$'),
                           DATA = c(rep(c('Prob J24 ($p(\\mathrm{GC}) > 0.05$)', 'L24 ($p(\\mathrm{GC}) > 0.05$)'), each = 3), rep('Binary J24', 2)))

R27_dat$DATA <- factor(R27_dat$DATA, levels = c("L24 ($p(\\mathrm{GC}) > 0.05$)", "Binary J24", "Prob J24 ($p(\\mathrm{GC}) > 0.05$)"))
ann_text_R27$DATA <-  factor(ann_text_R27$DATA, levels = c("L24 ($p(\\mathrm{GC}) > 0.05$)", "Binary J24", "Prob J24 ($p(\\mathrm{GC}) > 0.05$)"))
ann_text_R27$ID <- 'R27'

R27_dat$ID <- 'R27'

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

sum(R84_dat_Jans$p)

R84_dat_Jans_bin <- filter(Y_obs_J, (x - c[1, 1])^2 + (y - c[1,2])^2 < 7.5^2)
R84_dat_Jans_bin$p <- 1
R84_dat_Jans_bin <- R84_dat_Jans_bin %>%
  filter(p > 0 & M < 26.3 & C > 0.8 & C < 2.4) %>%
  dplyr::select(x, y, M, C, p) %>%
  mutate(DATA = 'Binary J24')

nrow(R84_dat_Jans_bin)

R84_dat <- bind_rows(R84_dat_harris, R84_dat_Jans, R84_dat_Jans_bin) %>%
  filter(C > 0.8)

ann_text_R84 <- data.frame(x = c(14.7, 14.7, 7.5, 14.4, 14.4, 7.5, 15.5, 7.5), y = c(59.5, 58.7, 47.5, 59.5, 58.7, 47.5, 59.1, 47.5), 
                       lab = c("$N_{\\mathrm{GC}}: 19\\pm 9$ (Ours)",
                               '$\\mu_{\\mathrm{TO}}: 25.96_{-0.46}^{+0.43}$ mag (Ours)',
                               '$N_{\\mathrm{cand}} = 18.0$',
                               "$N_{\\mathrm{GC}}: 32\\pm 13$ (Ours)",
                               '$\\mu_{\\mathrm{TO}}: 26.12_{-0.39}^{+0.43}$ mag (Ours)',
                               '$N_{\\mathrm{cand}} = 17.4$',
                               '$N_{\\mathrm{GC}}: 43\\pm 6$ (J24)',
                               '$N_{\\mathrm{cand}} = 26$'),
                       DATA = c(rep(c('Prob J24 ($p(\\mathrm{GC}) > 0.05$)', 'L24 ($p(\\mathrm{GC}) > 0.05$)'), each = 3), rep('Binary J24',2)))

R84_dat$DATA <- factor(R84_dat$DATA, levels = c("L24 ($p(\\mathrm{GC}) > 0.05$)", "Binary J24", "Prob J24 ($p(\\mathrm{GC}) > 0.05$)"))
ann_text_R84$DATA <-  factor(ann_text_R84$DATA, levels = c("L24 ($p(\\mathrm{GC}) > 0.05$)", "Binary J24", "Prob J24 ($p(\\mathrm{GC}) > 0.05$)"))
ann_text_R84$ID <- 'R84'

ann_text <- bind_rows(ann_text_R27, ann_text_R84)

R84_dat$ID <- 'R84'

LSBG_dat <- bind_rows(R27_dat, R84_dat)

tikz(file = "R27_R84_GC_F814W_plot.tex", standAlone=T, width = 8, height = 5)
ggplot(LSBG_dat, aes(x, y)) + geom_point(aes(color = M, size = p*1.5)) + facet_grid(ID~DATA, scales = 'free', switch = 'y') +
  scale_alpha_identity() + scale_size_identity() +
  geom_text(aes(x = x, y = y, label = lab), data = ann_text, size = 2.5) +
  scale_color_viridis_c(direction = -1, name = 'F814W') + 
  theme_bw() + 
  theme(axis.title=element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.justification = "left",
        legend.box.margin = margin(l = 0.1, unit = "cm"),
        aspect.ratio = 1) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
dev.off()
system('pdflatex R27_R84_GC_F814W_plot.tex')

res_prob_n <- readRDS('data/v10wfc3/res_prob_v10wfc3_seed_7904_n.RDS')
res_prob <- readRDS('data/v10wfc3/res_prob_v10wfc3_seed_7904.RDS')

res_prob_n_thinned <- res_prob_n[-c(1:20000),] %>%
  .[seq(1, 130000, by = 100),]

res_prob_thinned <- res_prob[-c(1:20000),] %>%
  .[seq(1, 130000, by = 100),]

e1 <- c(0.88, 1.35, 0.75)
theta1 <- c(0, -pi/12, 0)

e2 <- c(1.3, 1.35, 0.75)
theta2 <- c(pi/4, -pi/12, 0)

loc <- as.matrix(expand.grid(x = seq(11.38545 - 7.5, 11.38545 + 7.5, by = 0.24), y = seq(10.48364 - 7.5, 10.48364 + 7.5, by = 0.24)))

compute_mem_prob <- function(loc, res_prob, e, theta){
  prob <- numeric(nrow(res_prob))
  loc <- matrix(loc, ncol = 2)
  for (i in 1:nrow(res_prob)) {
    tot <- Sersic_ints(loc, c[1, ], res_prob[i, 'N_R27'], res_prob[i, 'R_R27'],
                       e[1], res_prob[i, 'n_R27'], theta[1]) + 
      Sersic_ints(loc, c[2, ], res_prob[i, 'N_W84'], res_prob[i, 'R_W84'],
                  e[2], res_prob[i, 'n_W84'], theta[2]) +
      Sersic_ints(loc, c[3, ], res_prob[i, 'N_W83'], res_prob[i, 'R_W83'],
                  e[3], res_prob[i, 'n_W83'], theta[3]) + res_prob[i, 'b0']
    
    prob[i] <- Sersic_ints(loc, c[1, ], res_prob[i, 'N_R27'], res_prob[i, 'R_R27'],
                        e[1], res_prob[i, 'n_R27'], theta[1])/tot 
  }
  return(mean(prob))
}


prob_dat_1 <- numeric(nrow(loc))
for (i in 1:nrow(loc)) {
  prob_dat_1[i] <- compute_mem_prob(loc[i,], res_prob_thinned, e1, theta1)
}

prob_dat_2 <- numeric(nrow(loc))
for (i in 1:nrow(loc)) {
  prob_dat_2[i] <- compute_mem_prob(loc[i,], res_prob_n_thinned, e2, theta2)
}

loc <- data.frame(loc)
p <- c(prob_dat_1)

R27_dat_loc <- bind_rows(loc[,c('x', 'y')]) %>%
  mutate(p = p)

tikz(file = "R27_GC_prob_plot.tex", standAlone=T, width = 4, height = 3)
ggplot(R27_dat_loc, aes(x,y)) + geom_contour_filled(aes(z = p)) +
  scale_fill_manual(values = wes_palette("Zissou1", n = 18, type = "continuous")[seq(1,18, by = 2)], 
                    name = '$\\hat{\\pi}_{\\mathrm{R27}}(s_i)$') +
  geom_point(data = R27_dat_harris, aes(x, y), size = 0.25) + coord_fixed() + 
  theme_minimal() + 
  theme(axis.title=element_blank(), 
        axis.text=element_blank(), 
        axis.ticks=element_blank(),
        strip.background = element_blank(),
        legend.justification = "left",
        legend.box.margin = margin(l = 0.1, unit = "cm")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, reverse = TRUE))
dev.off()
system('pdflatex R27_GC_prob_plot.tex')


#===============================================================================================#

R27_harris_mag <- filter(Y_obs, (x - c[1, 1])^2 + (y - c[1,2])^2 < 6^2)
R27_harris_mag$p <- rowMeans(R27_harris_mag[,8:507])
R27_harris_mag <- R27_dat_harris %>%
  dplyr::select(x, y, M, C, p)

R27_dat_Jans_mag <- filter(Y_obs_J, (x - c[1, 1])^2 + (y - c[1,2])^2 < 6^2)
R27_dat_Jans_mag <- R27_dat_Jans_mag %>%
  filter( M < 26.3 & C > 0.8 & C < 2.4) %>%
  dplyr::select(x, y, M, C)

hist(R27_dat_Jans_mag$M, breaks = 20)

log_lik_mag <- function(par, data = R27_harris_mag) {
  mu <- par[1]
  sig <- par[2]
  
  M <- data$M
  p <- data$p
  
  log_lik <- -sum(p*(dnorm(M, mu, sig, log = T) + log(f_cpp(M, 25.75)) - log(Phi_f_cpp(25.75, mu, sig))))
  
  return(log_lik)
}

res <- optim(c(25, 1), log_lik_mag, method = "L-BFGS-B", lower = c(24, 0.5), upper = c(26.7, 1.9), hessian = T)

pnorm((26.3 - res$par[1])/sqrt(solve(res$hessian)[1,1]), 0, 1, lower.tail = F)


R27_harris_mag_fix <- filter(Y_obs, (x - c[1, 1])^2 + (y - c[1,2])^2 < 6^2)
R27_harris_mag_fix <- R27_harris_mag_fix %>%
  dplyr::select(x, y, M, C) %>%
  filter(C > 1.0 & C < 2.4)


log_lik_mag_fix <- function(par, data = R27_harris_mag_fix) {
  mu <- par[1]
  sig <- par[2]
  
  M <- data$M
  p <- data$p
  
  log_lik <- -sum((dnorm(M, mu, sig, log = T) + log(f_cpp(M, 25.75)) - log(Phi_f_cpp(25.75, mu, sig))))
  
  return(log_lik)
}

res <- optim(c(25, 1), log_lik_mag, data = R27_harris_mag, method = "L-BFGS-B", lower = c(24, 0.5), upper = c(26.7, 1.9), hessian = T)

pnorm((26.3 - res$par[1])/sqrt(solve(res$hessian)[1,1]), 0, 1, lower.tail = F)

dat_bind <- bind_rows(R27_harris_mag_fix, W88_harris_mag_fix)


