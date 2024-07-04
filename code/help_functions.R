# Sersic intensity
Sersic_ints <- function(X, #point pattern
                        c, # center
                        N, # average number of points
                        R_eff, # effective radius
                        e = 1, # ellipticity
                        n = 0.5, theta = 0 # Sersic index and inclination angle
){
  bn <- zipfR::Rgamma.inv(2*n, 0.5)
  r <- (((X[,1] - c[1])*cos(theta) - (X[,2]-c[2])*sin(theta))^2/R_eff^2 + ((X[,1] - c[1])*sin(theta) + (X[,2]-c[2])*cos(theta))^2/R_eff^2/e^2)^0.5
  return(exp(log(N) + (2*n)*log(bn) - bn*(r)^(1/n) - log(n) - lgamma(2*n) - log(2*pi) - 2*log(R_eff)- log(e)))
}

expit <- function(x, a, b){
  v <- 1/(1+exp(-x))
  return(a + (b-a)*v)
}

logit <- function(u, a, b){
  p <- pmin(1, pmax(0, (u-a)/(b-a))) 
  return(pmin(999, pmax(-999,log(p/(1-p)))))
}
#integrate Sersic intensity
integrate_Sersic <- function(grid, c, N, R_eff, e = 1, n = 0.5, theta = 0){
  dx1 <- grid[[2]][1]
  dx2 <- grid[[2]][2]
  
  return(unname(sum(Sersic_ints(grid[[1]], c, N, R_eff, e, n, theta))*dx1*dx2))
}

f <- function(M, Lim, alpha = 1.5){
  1/(1 + exp(alpha*(M-Lim)))
}

#simulate point pattern from Sersic
simulate_Sersic <- function(N, c, R_eff, e, n, theta){
  bn <- zipfR::Rgamma.inv(2*n, 0.5)
  z <- runif(N)
  x <- zipfR::Rgamma.inv(2*n, z)
  R <- (x/bn)^n*R_eff
  t <- runif(N, 0, 2*pi)
  A <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),2,2)
  X <- solve(A, t(cbind(R*cos(t), R*sin(t)*e)))
  return(data.frame(x = c[1] + X[1,], y = c[2] + X[2,]))
} 

# Simulate GCLF
simulate_GCLF <- function(N, mu, sigma){
  return(rnorm(N, mu, sigma))
}

simulate_Y <- function(S, b0, c, N, R_eff, e, n, theta, mu, sigma){
  A <- sf::st_area(st_as_sf(S))
  K <- nrow(c)
  N_y <- rpois(K, N)
  N_b <- rpois(1, b0*A)
  if(N_b == 0){
    Y <- data.frame()
  }
  else{
    Y0 <- as.data.frame(sp::spsample(S, N_b, type = 'random'))
    Y <- data.frame(x = Y0[,1], y = Y0[,2], M = rnorm(N_b, mu[1], sigma[1]), id = rep(0,N_b)) 
  }
  for (i in 1:K) {
    if(N_y[i] == 0)
    {
      Y <- Y
    }
    else{
      Y <- bind_rows(Y, data.frame(as.data.frame(simulate_Sersic(N_y[i], c[i,], R_eff[i], e[i], n[i], theta[i])), M = rnorm(N_y[i], mu[i+1], sigma[i+1]), id = rep(i, N_y[i])))
    }
  }
  Y <- Y[unlist(sf::st_intersects(st_as_sf(S), st_as_sf(Y, coords = c('x','y')))),]
  return(Y)
}


err_M <- function(M, m0 = 25.5, a = 0.0884, b = 0.645){
  a*exp(b*(M - m0))
}

simulate_Y_noisy <- function(S, b0, c, N, R_eff, e, n, theta, mu, sigma, a = 0.0884, b = 0.645, m0 = 25.5){
  A <- sf::st_area(st_as_sf(S))
  K <- nrow(c)
  N_b <- floor(b0*A)
  if(N_b == 0){
    Y <- data.frame()
  }
  else{
    Y0 <- as.data.frame(sp::spsample(S, N_b, type = 'random'))
    M <- rnorm(N_b, mu[1], sigma[1])
    M_err <- rnorm(N_b, M, err_M(M))
    Y <- data.frame(x = Y0[,1], y = Y0[,2], M = M_err, Mt = M, id = rep(0,N_b), e_M = err_M(M, m0, a, b)) 
  }
  for (i in 1:K) {
    if(N[i] == 0)
    {
      Y <- Y
    }
    else{
      M <- rnorm(N[i], mu[i+1], sigma[i+1])
      M_err <- rnorm(N[i], M, err_M(M))
      Y <- bind_rows(Y, data.frame(as.data.frame(simulate_Sersic(N[i], c[i,], R_eff[i], e[i], n[i], theta[i])), M = M_err, Mt = M, id = rep(i, N[i]), e_M = err_M(M)))
    }
  }
  Y <- Y[unlist(sf::st_intersects(st_as_sf(S), st_as_sf(Y, coords = c('x','y')))),]
  return(Y)
}

simulate_Yf <- function(Y, Lim){
  return(Y[rbernoulli(nrow(Y), p = f(Y$Mt, Lim)),])
}

simulate_Yf_noise <- function(Y, Lim, alpha = 1.5){
  return(Y[rbernoulli(nrow(Y), p = f(Y$M, Lim, alpha)),])
}

simulate_Yobs <- function(YT, ps){
  Z <- rbernoulli(nrow(YT), ps)
  return(YT[Z,])
}

phi_f <- function(M, mu, sigma, Lim){
  dnorm(M, mu, sigma)*f(M, Lim)
}

psi_f <- function(M, mu, sigma, Lim, m0, b0, b1, a){
  phi_eM_cpp(M, mu, sigma, m0, b0, b1)*f(M, Lim, a)
}

log_lik_MATHPOP <- function(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, cf_error){
  
  alpha <- cf_error$alpha
  m50 <- cf_error$m50
  beta0 <- cf_error$beta0
  beta1 <- cf_error$beta1
  m1 <- cf_error$m1
  
  if(b0 <= 0 || any(N < 0) || any(R_eff < 0) || any(sigma < 0)){
    return(-Inf)
  }
  else{
    A <- sf::st_area(sf::st_as_sf(S))
    K <- nrow(c)
    Theta <- list()
    for (i in 1:K) {
      Theta[[i]] <- list(c = c[i,], N = N[i], R_eff = R_eff[i], e = e[i], n = n[i], theta = theta[i], mu = mu[i+1], sigma = sigma[i+1])
    }

    norm_const <- sum(unlist(lapply(Theta, function(x){integrate_Sersic(grid, x$c, x$N, x$R_eff, x$e, x$n, x$theta)*p*p_eM_cpp(Lim = m50, x$mu, x$sigma, m0 = m1, b0 = beta0, b1 = beta1, a = alpha)})))
    L <- - norm_const - b0*A*p*p_eM_cpp(Lim = m50, mu[1], sigma[1], m0 = m1, b0 = beta0, b1 = beta1, a = alpha)
    loc <- b0*p*psi_f(Y_obs$M, mu[1], sigma[1], Lim = m50, m0 = m1, b0 = beta0, b1 = beta1, a = alpha)
    loc <- loc + vapply(data.table::transpose(lapply(Theta, function(x){Sersic_ints(Y_obs[,c('x','y')], x$c, x$N, x$R_eff, x$e, x$n, x$theta)*p*psi_f(Y_obs$M, x$mu, x$sigma, Lim = m50, m0 = m1, b0 = beta0, b1 = beta1, a = alpha)})), sum, 0)
    L <- L + sum(log(loc)) 
    return(L)
  }
}

# log_post_MATHPOP <- function(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, prior, cf_error){
#     lp <- log_lik_mod_np(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, cf_error) +
#       dlnorm(b0, prior$b0[1], prior$b0[2], log = T) +
#       sum(dnorm(mu, prior$mu[,1], prior$mu[,2], log = T)) +
#       sum(dlnorm(sigma, prior$sigma[,1], prior$sigma[,2], log = T)) + sum(log(sigma)) +
#       sum(dfoldnorm(N, mean = prior$N[,1], sd = prior$N[,2], log = T)) + sum(log(N)) +
#       sum(dlnorm(R_eff, prior$R_eff[,1], prior$R_eff[,2], log = T)) + sum(log(R_eff)) +
#       sum(dlnorm(n, prior$n[,1], prior$n[,2], log = T)) + sum(log(n))
#     return(lp)
# }

log_post_MATHPOP <- function(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, prior, Ng, cf_error){
  if(Ng == 0)
  {
    lp <- log_lik_MATHPOP(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, cf_error) +
      dlnorm(b0, prior$IGM$b0[1], prior$IGM$b0[2], log = T) +
      dnorm(mu[1], prior$IGM$mu[1], prior$IGM$mu[2], log = T) +
      dlnorm(sigma[1], prior$IGM$sigma[1], prior$IGM$sigma[2], log = T) + 
      sum(dnorm(mu[-1], prior$UDG$mu[,1], prior$UDG$mu[,2], log = T)) +
      sum(dlnorm(sigma[-1], prior$UDG$sigma[,1], prior$UDG$sigma[,2], log = T)) + sum(log(sigma)) +
      sum(dfoldnorm(N, mean = prior$UDG$N[,1], sd = prior$UDG$N[,2], log = T)) + sum(log(N)) +
      sum(dlnorm(R_eff, prior$UDG$R_eff[,1], prior$UDG$R_eff[,2], log = T)) + sum(log(R_eff)) +
      sum(dlnorm(n, prior$UDG$n[,1], prior$UDG$n[,2], log = T)) + sum(log(n))
    return(lp)
  }
  else{
    K <- nrow(c)
    lp <- log_lik_MATHPOP(S, grid, b0, c, N, R_eff, e, n, theta, mu, sigma, Y_obs, p, cf_error) +
      dlnorm(b0, prior$IGM$b0[1], prior$IGM$b0[2], log = T) +
      dnorm(mu[1], prior$IGM$mu[1], prior$IGM$mu[2], log = T) +
      dlnorm(sigma[1], prior$IGM$sigma[1], prior$IGM$sigma[2], log = T) + 
      sum(dnorm(mu[-1], c(prior$gal$mu[,1], prior$UDG$mu[,1]), c(prior$gal$mu[,2], prior$UDG$mu[,2]), log = T)) + 
      sum(dlnorm(sigma[-1], c(prior$gal$sigma[,1], prior$UDG$sigma[,1]), c(prior$gal$sigma[,2], prior$UDG$sigma[,2]), log = T)) + sum(log(sigma)) +
      sum(dlnorm(N[1:Ng], prior$gal$N[,1], prior$gal$N[,2], log = T)) + 
      sum(dfoldnorm(N[(Ng+1):K], prior$UDG$N[,1], prior$UDG$N[,2], log = T)) + sum(log(N)) +
      sum(dlnorm(R_eff, c(prior$gal$R_eff[,1], prior$UDG$R_eff[,1]), c(prior$gal$R_eff[,2], prior$UDG$R_eff[,2]), log = T)) + sum(log(R_eff)) +
      sum(dlnorm(n, c(prior$gal$n[,1], prior$UDG$n[,1]), c(prior$gal$n[,2], prior$UDG$n[,2]), log = T)) + sum(log(n))
    return(lp)
  }
}

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


meas_uncertain_mix_func <- function(dat, n_iter, seed = 12345){
  set.seed(seed)
  p_mat <- matrix(0, ncol = n_iter, nrow = nrow(dat))
  par_mat <- matrix(0, ncol = 9, nrow = n_iter)
  sim_dat <- data.frame()
  for(i in 1:n_iter){
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
  return(list(prob = p_mat, par = par_mat, sim = sim_dat))
}
