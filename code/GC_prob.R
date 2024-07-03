library(tidyverse)
library(tikzDevice)
library(mixtools)

fnames <- list.files('data/point_source_data/')
ps_harris_fnames <- fnames[!grepl('Jans', fnames)]
dat_harris <- data.frame()
for (i in ps_harris_fnames) {
  df <- read_csv(paste0('data/point_source_data/',i))
  df <- bind_cols(df, field = rep(i, nrow(df)))
  dat_harris <- bind_rows(dat_harris, df)
}

tikz(file = "CMD_box_GC_harris.tex", standAlone=T, width = 3, height = 3)
ggplot(dat_harris, aes(F475W - F814W, F814W)) + geom_point(size = 0.05) + scale_y_reverse() +
  annotate('segment', x = 1.0, xend = 2.4, y = 26.5, yend = 26.5, colour = 'red') +
  annotate('segment', x = 1.0, xend = 1.0, y = 22, yend = 26.5, colour = 'red') +
  annotate('segment', x = 2.4, xend = 2.4, y = 22, yend = 26.5, colour = 'red') +
  annotate('segment', x = 1.0, xend = 2.4, y = 22, yend = 22, colour = 'red') +
  theme_minimal()
dev.off()
system('pdflatex CMD_box_GC_harris.tex')


CMD_harris <- mutate(dat_harris, C = F475W - F814W, M = F814W) %>%
  dplyr::select(x, y, RA, DEC, C, M, F814W, F475W, field) %>%
  filter(M > 22 & C > 0)

fname_GC_prob <- list.files('data/GC_prob/')
Harris_fnames <- fname_GC_prob[!grepl('Jans', fname_GC_prob)]
GC_prob_ls_Harris <- list()
for(i in 1:length(Harris_fnames)){
  ls <- readRDS(paste0('data/GC_prob/',Harris_fnames[i]))
  GC_prob_ls_Harris[[i]] <- ls
}

GC_probability_Harris <- do.call(cbind,lapply(GC_prob_ls_Harris, function(x) x$posterior[, which.max(x$lambdahat)]))

colnames(GC_probability_Harris) <- c(paste0('p', 1:500))

CMD_harris <- bind_cols(CMD_harris, GC_probability_Harris)

CMD_harris$p <- rowMeans(GC_probability_Harris)
CMD_harris$p_sd <- apply(GC_probability_Harris, 1, sd)

CMD_exclude <- mutate(dat_harris, C = F475W - F814W, M = F814W) %>%
  dplyr::select(x, y, RA, DEC, C, M, F814W, F475W, field) %>%
  filter(M < 22 | C < 0)

tikz(file = "GC_prob_harris_cropped.tex", standAlone=T, width = 3, height = 3)
ggplot(CMD_harris, aes(C, M, color = p)) + geom_point(aes(size = p_sd*2.5)) +
  #geom_point(data = CMD_exclude, aes(C, M), size = 0.25, color = 'grey') +
  scale_size_identity() +
  scale_y_reverse() + theme_minimal() +
  xlab('F475W - F814W') + ylab('F814W') +
  scale_colour_viridis_c(name = '$p(\\mathrm{GC})$', limits = c(0,1)) +
  annotate('segment', x = 1.0, xend = 2.4, y = 26.5, yend = 26.5, colour = 'red') +
  annotate('segment', x = 1.0, xend = 1.0, y = 22, yend = 26.5, colour = 'red') +
  annotate('segment', x = 2.4, xend = 2.4, y = 22, yend = 26.5, colour = 'red') +
  annotate('segment', x = 1.0, xend = 2.4, y = 22, yend = 22, colour = 'red') + xlim(c(0,3))+
  theme(legend.position="none")
dev.off()
system('pdflatex GC_prob_harris_cropped.tex')

for (i in fnames) {
  write_csv(filter(CMD_harris, field == i), paste0('data/prob_GC_data/',unlist(strsplit(i, split = '.', fixed = T))[1],'_pGC.csv'))
}

#dat <- mutate(dat, C = F475W - F814W, M = F814W) %>%
#  dplyr::select(x, y, RA, DEC, C, M, field) %>%
#  filter(M > 22 & C < 2.4 & C > 1)
#for (i in fnames) {
#  write_csv(filter(dat, field == i), paste0('data/',unlist(strsplit(i, split = '.', fixed = T))[1],'_npGC.csv'))
#}

# Janssens' fake star test
ps_Jans_fnames <- fnames[grepl('Jans', fnames)]
dat_Jans <- read_csv(paste0('data/point_source_data/',ps_Jans_fnames))

tikz(file = "CMD_box_GC_Jans.tex", standAlone=T, width = 3, height = 3)
ggplot(dat_Jans, aes(F475W - F814W, F814W)) + geom_point(size = 0.05) + scale_y_reverse() + 
  annotate('segment', x = 0.8, xend = 2.4, y = 26.3, yend = 26.3, colour = 'red') +
  annotate('segment', x = 0.8, xend = 0.8, y = 22, yend = 26.3, colour = 'red') +
  annotate('segment', x = 2.4, xend = 2.4, y = 22, yend = 26.3, colour = 'red') +
  annotate('segment', x = 0.8, xend = 2.4, y = 22, yend = 22, colour = 'red') +
  xlim(c(0,3)) +
  theme_minimal()
dev.off()
system('pdflatex CMD_box_GC_Jans.tex')

 
single_func <- function(par, dat){
  GCLF_TO <- par[1]
  GCLF_sig <- par[2]

  M <- dat$F814W

  nll <- -sum(log(truncnorm::dtruncnorm(M, -Inf, 26.3, GCLF_TO, GCLF_sig)))

  return(nll)
}

dat_trunc <- filter(dat_Jans, F814W < 26.3 & C > 0.8 & C < 2.4)

res_single <- optim(c(26.3, 1.2), single_func, dat = dat_trunc, 
                                        method = "L-BFGS-B", lower = c(24, 1), upper = c(50, 10), hessian = T)
# 
# two_func <- function(par, dat){
#   w <- par[1]
#   GCLF_TO <- par[2]
#   GCLF_sig <- par[3]
#   mu <- par[4]
#   sigma <- par[5]
#   
#   M <- dat$F814W
#   
#   nll <- -sum(log(w*dnorm(M, GCLF_TO, GCLF_sig)*f_cpp(M, 26.69, 6.56)/
#                     Phi_f_cpp(26.69, GCLF_TO, GCLF_sig, a = 6.56) + 
#                     (1-w)*dnorm(M, mu, sigma)))
#   
#   return(nll)
# }
# 
# 
# res_single <- optim(c(26.3, 1.2), single_func, dat = dat, 
#                     method = "L-BFGS-B", lower = c(24, 1), upper = c(30, 3))
# 
# res_two <- optim(c(0.6, 26.3, 1.2, 26.3, 0.2), two_func, dat = dat, 
#                  method = "L-BFGS-B", lower = c(0, 26, 1, 25, 0.1), upper = c(1, 27, 2, 27, 0.3))
# 
# # compute BIC = -422.5; two-component mixture is clearly better than a single GCLF component
# 2*res_two$value + log(nrow(dat))*5 - 2*res_single$value + log(nrow(dat))*2


mix_func <- function(par, dat){
  w <- par[1]
  wr <- par[2]
  GCLF_TO <- 26.3
  GCLF_sig <- par[3]
  GC_color_mu_r <- par[4]
  GC_color_sig_r <- par[5]
  GC_color_mu_b <- par[6]
  GC_color_sig_b <- par[7]
  mu <- par[8]
  sigma <- par[9]
  
  M <- dat$F814W
  C <- dat$C
  
  nll <- -sum(log(w*dnorm(M, GCLF_TO, GCLF_sig)*f_cpp(M, 26.69, 6.56)/Phi_f_cpp(26.69, GCLF_TO, GCLF_sig, a = 6.56)*(wr*dnorm(C, GC_color_mu_r, GC_color_sig_r) + (1-wr)*dnorm(C, GC_color_mu_b, GC_color_sig_b)) + 
                    (1-w)*dnorm(M, mu, sigma)*dunif(C, 0.8, 2.4)))
  
  return(nll)
}

set.seed(123456)
p_mat <- matrix(0, ncol = 500, nrow = nrow(dat_Jans))
par_mat <- matrix(0, ncol = 9, nrow = 500)
sim_dat <- data.frame()
for(i in 1:500){
  sim_CM <- as.data.frame(t(apply(dat_Jans[,c('C', 'F814W', 'M_err', 'C_err')], 1, 
                                  function(x){MASS::mvrnorm(n = 1, x[c(1,2)], Sigma = matrix(c(x[4]^2, x[3]^2, x[3]^2, x[3]^2), 2))})))
  
  dat <- filter(sim_CM, C > 0.8 & C < 2.4)[,c('C', 'F814W')]
  sim_dat <- bind_rows(sim_dat, dat)
  M <- dat$F814W
  C <- dat$C
  idx <- which(sim_CM$C > 0.8 & sim_CM$C < 2.4)
  res <- optim(c(0.6, 0.5, 1.2, 1.5, 0.16, 1.6, 0.1, 26.3, 0.3), mix_func, dat = dat, 
               method = "L-BFGS-B", lower = c(0, 0, 0.9, 1.25, 0.1, 1.5, 0.1, 26, 0.1), upper = c(1, 1, 1.5, 1.5, 0.5, 1.8, 0.5, 27, 0.6))$par
  
  par_mat[i,] <- res
  
  w <- res[1]
  wr <- res[2]
  GCLF_TO <- 26.3
  GCLF_sig <- res[3]
  GC_color_mu_r <- res[4]
  GC_color_sig_r <- res[5]
  GC_color_mu_b <- res[6]
  GC_color_sig_b <- res[7]
  mu <- res[8]
  sigma <- res[9]
  
  GCLF_comp <- w*dnorm(M, GCLF_TO, GCLF_sig)*f_cpp(M, 26.69, 6.56)/Phi_f_cpp(26.69, GCLF_TO, GCLF_sig, a = 6.56)*
    (wr*dnorm(C, GC_color_mu_r, GC_color_sig_r) + (1-wr)*dnorm(C, GC_color_mu_b, GC_color_sig_b))
  cont_comp <- (1-w)*dnorm(M, mu, sigma)*dunif(C, 0.8, 2.4)
  
  prob <- GCLF_comp/(GCLF_comp + cont_comp)
    
  p_mat[idx,i] <- prob
}

dat_Jans$p <- rowMeans(p_mat)
dat_Jans$p_sd <- apply(p_mat, 1, sd)

res_sd <- apply(par_mat, 2, sd)

w <- res[1]
wr <- res[2]
GCLF_TO <- 26.3
GCLF_sig <- res[3]
GC_color_mu_r <- res[4]
GC_color_sig_r <- res[5]
GC_color_mu_b <- res[6]
GC_color_sig_b <- res[7]
mu <- res[8]
sigma <- res[9]

Phi <- Phi_f_cpp(26.69, GCLF_TO, GCLF_sig, a = 6.56)

tikz(file = "F814W_mixture_Jans.tex", standAlone=T, width = 4, height = 3)
ggplot(sim_dat, aes(F814W)) + geom_histogram(aes(y = after_stat(count)/sum(count)/.05253/2),breaks = seq(22, 27.2, length = 50), alpha = 0.5) +
  geom_function(fun = function(x) w*dnorm(x, GCLF_TO, GCLF_sig)*f_cpp(x, 26.69, 6.56)/Phi + (1-w)*dnorm(x, mu, sigma)) +
  geom_function(fun = function(x) w*dnorm(x, GCLF_TO, GCLF_sig)*f_cpp(x, 26.69, 6.56)/Phi, color = "#EBCC2A") +
  geom_function(fun = function(x) (1-w)*dnorm(x, mu, sigma), color = 'purple') + ylab('Density') + theme_minimal() 
dev.off()
system('pdflatex F814W_mixture_Jans.tex')

tikz(file = "color_mixture_Jans.tex", standAlone=T, width = 4, height = 3)
ggplot(sim_dat, aes(C)) + geom_histogram(aes(y = after_stat(count)/sum(count)/0.032),breaks = seq(0.8, 2.4, length = 50), alpha = 0.5) +
  geom_function(fun = function(x) w*(wr*dnorm(x, GC_color_mu_r, GC_color_sig_r) + (1-wr)*dnorm(x, GC_color_mu_b, GC_color_sig_b)) + (1-w)*dunif(x, 0.8, 2.4)) +
  geom_function(fun = function(x) w*(wr*dnorm(x, GC_color_mu_r, GC_color_sig_r) + (1-wr)*dnorm(x, GC_color_mu_b, GC_color_sig_b)), color = "#EBCC2A") +
  geom_function(fun = function(x) w*wr*dnorm(x, GC_color_mu_r, GC_color_sig_r), color = '#3B9AB2') +
  geom_function(fun = function(x) w*(1-wr)*dnorm(x, GC_color_mu_b, GC_color_sig_b), color = "#FF0000") +
  geom_function(fun = function(x) (1-w)*dunif(x, 0.8, 2.4), color = 'purple') + ylab('Density') + theme_minimal() + xlab('F475W - F814W')
dev.off()
system('pdflatex color_mixture_Jans.tex')

tikz(file = "GC_prob_Jans.tex", standAlone=T, width = 4, height = 3)
ggplot(dat_Jans, aes(C, F814W, color = p)) + geom_point(aes(size = p_sd*2.5)) +
  scale_size_identity() +
  scale_y_reverse() + theme_minimal() +
  xlab('F475W - F814W') + ylab('F814W') +
  scale_colour_viridis_c(name = '$p(\\mathrm{GC})$') +
  annotate('segment', x = 0.8, xend = 2.4, y = 26.3, yend = 26.3, colour = 'red') +
  annotate('segment', x = 0.8, xend = 0.8, y = 22, yend = 26.3, colour = 'red') +
  annotate('segment', x = 2.4, xend = 2.4, y = 22, yend = 26.3, colour = 'red') +
  annotate('segment', x = 0.8, xend = 2.4, y = 22, yend = 22, colour = 'red') + xlim(c(0,3))
dev.off()
system('pdflatex GC_prob_Jans.tex')

frac_CMD_harris <- CMD_harris %>%
  dplyr::select(C, M, p, p_sd) %>%
  mutate(data = 'DOLPHOT')

frac_CMD_Jans <- dat_Jans %>%
  dplyr::select(C, F814W, p, p_sd) %>%
  mutate(data = 'SExtractor (J24)', M = F814W) %>%
  dplyr::select(C, M, p, p_sd, data)

set.seed(123456)
frac_CMD <- bind_rows(frac_CMD_harris[sample(1:nrow(frac_CMD_harris), 0.5*nrow(frac_CMD_harris), replace = F),],
                      frac_CMD_Jans[sample(1:nrow(frac_CMD_Jans), 0.5*nrow(frac_CMD_Jans), replace = F),])

tikz(file = "GC_prob.tex", standAlone=T, width = 6, height = 3)
ggplot(frac_CMD, aes(C, M, color = p)) + geom_point(aes(size = p_sd*2.5)) +
  scale_size_identity() +
  scale_y_reverse() + theme_minimal() +
  xlab('F475W - F814W') + ylab('F814W') +
  scale_colour_viridis_c(name = '$p(\\mathrm{GC})$') + 
  geom_segment(data = data.frame(xmin = c(0.8, 0.8, 2.4, 0.8), 
                                               xmax = c(2.4, 0.8, 2.4, 2.4), 
                                               ymin = c(26.3, 22, 22, 22), 
                                               ymax = c(26.3, 26.3, 26.3, 22),
                                 data = 'SExtractor (J24)'), 
           aes(x = xmin, xend = xmax, y = ymin, yend = ymax), colour = 'red') +
  geom_segment(data = data.frame(xmin = c(1, 1, 2.4, 1), 
                                 xmax = c(2.4, 1, 2.4, 2.4), 
                                 ymin = c(26.5, 22, 22, 22), 
                                 ymax = c(26.5, 26.5, 26.5, 22),
                                 data = 'DOLPHOT'), 
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax), colour = 'red') +
  facet_wrap(.~data)
dev.off()
system('pdflatex GC_prob.tex')

dat_Jans <- dat_Jans %>%
  dplyr::select(x, y, RA, DEC, C, F814W, field) %>%
  mutate(M = F814W) %>%
  dplyr::select(x, y, RA, DEC, C, M, field)

colnames(p_mat) <- c(paste0('p', 1:500))
p_mat <- data.frame(p_mat)

dat_Jans <- bind_cols(dat_Jans, p_mat)

for (i in unique(dat_Jans$field)) {
  write_csv(filter(dat_Jans, field == i), paste0('data/prob_GC_data/', i,'_pGC_Jans.csv'))
}


J_err <- read.table('data/binary_GC_data/Janssens_error.dat')
colnames(J_err) <- J_err[1,]
J_err <- J_err[-1,]

J_err <- J_err %>%
  mutate_at(c(1, 2, 3), as.numeric)

summary(lm(log(colourerr) ~ I(magi - 25.5), data = J_err))
# beta_0 = exp(-2.5933) = 0.07477
# beta_1 = 0.75094
J_err$magi475

ggplot(J_err, aes(magi, log(colourerr))) + geom_point() + geom_abline(slope = 0.75094, intercept = -21.74)

fake_v6acs_F814W <- readxl::read_xlsx('data/fake_data/v6_acs_fakestars.xlsx', sheet = 'F814W', col_types = 'numeric')
fake_v12acs_F814W <- readxl::read_xlsx('data/fake_data/v12_acs_fakestars.xlsx', sheet = 'F814W', col_types = 'numeric')
fake_v15acs_F814W <- readxl::read_xlsx('data/fake_data/v15_acs_fakestars.xlsx', sheet = 'F814W', col_types = 'numeric')
fake_acs_F814W <- bind_rows(fake_v6acs_F814W, fake_v12acs_F814W, fake_v15acs_F814W) %>%
  mutate(field = c(rep('v6acs', 47), rep('v12acs', 48), rep('v15acs', 46))) %>%
  group_by(min) %>%
  mutate(mag = sum(`<mag>`*nin)/sum(nin) - 0.249, 
         m_f = sum(nout)/sum(nin), m_ef = sqrt(sum(ef^2)/3),
         m_sigma = sum(sigma)/3) %>%
  ungroup() %>%
  dplyr::select(mag, m_f, m_ef, m_sigma) %>%
  mutate(t_mag = mag - 25.5, lg_sigma = log(m_sigma)) %>%
  distinct()


ggplot(fake_acs_F814W, aes(mag, m_f)) + geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = m_f - m_ef, ymax = m_f + m_ef)) + geom_function(fun = function(x) 1/(1+exp(1.498*(x - 25.75)))) +
  geom_function(fun = function(x) 1/(1+exp(6.5*(x - 26.6))), color = 'red')

summary(nls(m_f ~ 1/(1 + exp(a*(mag - m0))), data = fake_acs_F814W, start = list(a = 1.5, m0 = 26)))
summary(lm(lg_sigma ~ t_mag, data = fake_acs_F814W)) 

fake_v6acs_F475W <- readxl::read_xlsx('data/fake_data/v6_acs_fakestars.xlsx', sheet = 'F475W', col_types = 'numeric')
fake_v12acs_F475W <- readxl::read_xlsx('data/fake_data/v12_acs_fakestars.xlsx', sheet = 'F475W', col_types = 'numeric')
fake_v15acs_F475W <- readxl::read_xlsx('data/fake_data/v15_acs_fakestars.xlsx', sheet = 'F475W', col_types = 'numeric')
fake_acs_F475W <- bind_rows(fake_v6acs_F475W, fake_v12acs_F475W, fake_v15acs_F475W) %>%
  mutate(field = c(rep('v6acs', 34), rep('v12acs', 36), rep('v15acs', 36))) %>%
  group_by(min) %>%
  mutate(mag = sum(`<mag>`*nin)/sum(nin) - 0.533, 
         m_f = sum(nout)/sum(nin), m_ef = sqrt(sum(ef^2)/3),
         m_sigma = sum(sigma)/3) %>%
  ungroup() %>%
  dplyr::select(mag, m_f, m_ef, m_sigma) %>%
  mutate(t_mag = mag - 26, lg_sigma = log(m_sigma)) %>%
  distinct()


ggplot(fake_acs_F475W, aes(mag, m_f)) + geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = m_f - m_ef, ymax = m_f + m_ef)) + geom_function(fun = function(x) 1/(1+exp(1.66*(x - 26.967))))

summary(nls(m_f ~ 1/(1 + exp(a*(mag - m0))), data = fake_acs_F475W, start = list(a = 1.62, m0 = 27.5)))
summary(lm(lg_sigma ~ t_mag, data = fake_acs_F475W))

acs_fake <- bind_rows(fake_acs_F475W, fake_acs_F814W) %>%
  mutate(Filter = c(rep('F475W', 36), rep('F814W', 48)))

tikz(file = "err_acs.tex", standAlone=T, width = 6, height = 4)
ggplot(acs_fake, aes(mag, lg_sigma, colour = Filter)) + geom_point(size=0.5) +
  geom_function(fun = function(x) -2.55165 + 0.699*(x- 26), color = '#F8766D') +
  geom_function(fun = function(x) -2.426 + 0.645*(x - 25.5)) +
  xlab('Magnitude') + ylab('$\\log(\\sigma)$') + theme_minimal()
dev.off()
system('pdflatex err_acs.tex')

fake_v6wfc_F814W <- readxl::read_xlsx('data/fake_data/v6_wfc3_fakestars.xlsx', sheet = 'F814W', col_types = 'numeric')
fake_v12wfc_F814W <- readxl::read_xlsx('data/fake_data/v12_wfc3_fakestars.xlsx', sheet = 'F814W', col_types = 'numeric')
fake_v15wfc_F814W <- readxl::read_xlsx('data/fake_data/v15_wfc3_fakestars.xlsx', sheet = 'F814W', col_types = 'numeric')
fake_wfc_F814W <- bind_rows(fake_v6wfc_F814W, fake_v12wfc_F814W, fake_v15wfc_F814W) %>%
  mutate(field = c(rep('v6wfc', 48), rep('v12wfc', 45), rep('v15wfc', 48))) %>%
  group_by(min) %>%
  mutate(mag = sum(`<mag>`*nin)/sum(nin) - 0.249, 
         m_f = sum(nout)/sum(nin), m_ef = sqrt(sum(ef^2)/3),
         m_sigma = sum(sigma)/3) %>%
  ungroup() %>%
  dplyr::select(mag, m_f, m_ef, m_sigma) %>%
  mutate(t_mag = mag - 25.5, lg_sigma = log(m_sigma)) %>%
  distinct()

fake_wfc_F814W <- fake_wfc_F814W[-48,]

summary(nls(m_f ~ 1/(1 + exp(a*(mag - m0))), data = fake_wfc_F814W, start = list(a = 1.5, m0 = 26)))
summary(lm(lg_sigma ~ t_mag, data = fake_wfc_F814W)) 

fake_v6wfc_F475X <- readxl::read_xlsx('data/fake_data/v6_wfc3_fakestars.xlsx', sheet = 'F475X', col_types = 'numeric')
fake_v12wfc_F475X <- readxl::read_xlsx('data/fake_data/v12_wfc3_fakestars.xlsx', sheet = 'F475X', col_types = 'numeric')
fake_v15wfc_F475X <- readxl::read_xlsx('data/fake_data/v15_wfc3_fakestars.xlsx', sheet = 'F475X', col_types = 'numeric')
fake_wfc_F475X <- bind_rows(fake_v6wfc_F475X, fake_v12wfc_F475X, fake_v15wfc_F475X) %>%
  mutate(field = c(rep('v6wfc', 36), rep('v12wfc', 36), rep('v15wfc', 36))) %>%
  group_by(min) %>%
  mutate(mag = sum(`<mag>`*nin)/sum(nin) - 0.51, 
         m_f = sum(nout)/sum(nin), m_ef = sqrt(sum(ef^2)/3),
         m_sigma = sum(sigma)/3) %>%
  mutate(mag = mag + 0.16*1.5) %>%
  ungroup() %>%
  dplyr::select(mag, m_f, m_ef, m_sigma) %>%
  mutate(t_mag = mag - 26, lg_sigma = log(m_sigma)) %>%
  distinct()

fake_wfc_F475X$lg_sigma[c(29:32)] <- NaN

summary(nls(m_f ~ 1/(1 + exp(a*(mag - m0))), data = fake_wfc_F475X, start = list(a = 1.5, m0 = 26)))
summary(lm(lg_sigma ~ t_mag, data = fake_wfc_F475X)) 

wfc_fake <- bind_rows(fake_wfc_F475X, fake_wfc_F814W) %>%
  mutate(Filter = c(rep('F475W', 36), rep('F814W', 47)))

fake_test <- bind_rows(acs_fake, wfc_fake) %>%
  mutate(Camera = c(rep('ACS', 84), rep('WFC3', 83)))

m <- seq(20, 31.5, length = 114)
df_fitted <- data.frame(mag = rep(m, 4), 
                        m_f = c(1/(1+exp(1.57*(m - 26.52))),
                                1/(1+exp(1.5*(m - 25.75))),
                                1/(1+exp(1.23*(m - 28.02))),
                                1/(1+exp(1.66*(m - 26.93)))),
                        Filter = c(rep('F814W', 228), rep('F475W', 228)),
                        Camera = rep(c('WFC3', 'ACS', 'WFC3', 'ACS'), each = 114))


tikz(file = "fm.tex", standAlone=T, width = 6, height = 3)
ggplot(fake_test, aes(mag, m_f, colour = Filter)) + geom_point(size= 0.5) + 
  geom_errorbar(aes(ymin = m_f - m_ef, ymax = m_f + m_ef)) + 
  geom_line(data = df_fitted, aes(mag, m_f)) + 
  xlab('Magnitude') + ylab('$f$') + theme_minimal() + facet_wrap(~Camera) 
dev.off()
system('pdflatex fm.tex')

m <- seq(19.5, 31, length = 114)
err_fitted <- data.frame(mag = rep(m, 4), 
                         lg_sigma = c(-2.426 + 0.645*(m - 25.5),
                                -2.552 + 0.699*(m - 26),
                                -2.326 + 0.613*(m - 25.5),
                                -2.912 + 0.652*(m - 26)),
                        Camera = c(rep('ACS', 228), rep('WFC3', 228)),
                        Filter = rep(c('F814W', 'F475W', 'F814W', 'F475W'), each = 114))

tikz(file = "err_m.tex", standAlone=T, width = 6, height = 3)
ggplot(fake_test, aes(mag, lg_sigma, colour = Filter)) + geom_point(size=0.5) +
  geom_line(data = err_fitted, aes(mag, lg_sigma)) +
  xlab('Magnitude') + ylab('$\\log(\\sigma)$') + theme_minimal() +facet_wrap(~ Camera)
dev.off()
system('pdflatex err_m.tex')

