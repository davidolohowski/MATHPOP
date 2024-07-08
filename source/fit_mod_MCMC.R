source('source/help_functions.R')
library(progress)

fit_mod_np <- function(S, grid, Theta, UDG_ID, c, e, theta, Y, p, prior, M, tune, n = 1000, gamma = 0.1, prob_model = TRUE, cf_error){
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta || Current MCMC iteration: :current]",
                         total = M,
                         complete = "=",
                         incomplete = "-",
                         current = ">",
                         clear = FALSE,
                         width = 200) 
  K <- nrow(c)
  res <- matrix(0, nrow = M+1, ncol = 5*K+3)
  col_names <- c('b0', paste0('N_', UDG_ID), paste0('R_', UDG_ID), paste0('n_', UDG_ID), 'mu_0', paste0('mu_', UDG_ID), 'sigma_0', paste0('sigma_', UDG_ID))
  colnames(res) <- col_names
  res[1,] <- c(Theta$b0, Theta$N, Theta$R_eff, Theta$n, Theta$mu, Theta$sigma)
  t_b <- tune$b0
  t_N <- tune$N
  t_R <- tune$R
  t_n <- tune$n
  t_m <- tune$mu
  t_s <- tune$sigma
  gamma <- 2.38^2/(5*K+3)
  idx <- sample(1:500, M, replace = T)
  if(prob_model){
    Y_dat <- Y[, 1:7]
    Y_prob <- as.matrix(Y[, 8:507]) 
  }
  for (i in 1:M) {
    if(i <= n){
      if(prob_model){
        Y_obs <- Y_dat[rbernoulli(nrow(Y_dat), p = Y_prob[,idx[i]]),]
      }
      else{
        Y_obs <- Y
      }
      b0_old <- res[i, 1]
      N_old <- exp(res[i, 2:(K+1)])
      R_eff_old <- exp(res[i,(K+2):(2*K+1)])
      n_old <- exp(res[i,(2*K+2):(3*K+1)])
      mu_old <- res[i, (3*K+2):(4*K+2)]
      sigma_old <- exp(res[i, (4*K+3):(5*K+3)])
      
      b0_new <- rnorm(1, b0_old, t_b)
      N_new <- exp(rnorm(K, log(N_old), t_N))
      R_eff_new <- exp(rnorm(K, log(R_eff_old), t_R))
      n_new <- exp(rnorm(K, log(n_old), t_n))
      mu_new <- rnorm(K+1, mu_old, t_m)
      sigma_new <- exp(rnorm(K+1, log(sigma_old), t_s))
      
      r <- log_post_mod_np(S, grid, b0_new, c, N_new, R_eff_new, e, n_new, theta, mu_new, sigma_new, Y_obs, p, prior, cf_error = cf_error) - 
        log_post_mod_np(S, grid, b0_old, c, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, cf_error = cf_error)
      
      C <- cov(res[c(1:n),])
      
      Mu <- colMeans(res[c(1:n),])
      
      if(log(runif(1)) < r){
        res[i+1, ] <- c(b0_new, log(N_new), 
                        log(R_eff_new), log(n_new), mu_new, log(sigma_new))
      }
      else{
        res[i+1,] <- res[i,]
      }
      pb$tick()
    }
    else{
      if(prob_model){
        Y_obs <- Y_dat[rbernoulli(nrow(Y_dat), p = Y_prob[,idx[i]]),]
      }
      else{
        Y_obs <- Y
      }
      b0_old <- res[i, 1]
      N_old <- exp(res[i, 2:(K+1)])
      R_eff_old <- exp(res[i,(K+2):(2*K+1)])
      n_old <- exp(res[i,(2*K+2):(3*K+1)])
      mu_old <- res[i, (3*K+2):(4*K+2)]
      sigma_old <- exp(res[i, (4*K+3):(5*K+3)])
      
      X <- res[i,]
      C <- (i-2)/(i-1)*C + (X-Mu)%*%t(X-Mu)/i
      Mu <- (Mu*(i-1) + X)/i
      
      Theta_new <- MASS::mvrnorm(1, mu = X, 
                                 Sigma = gamma*(C + diag(1e-7, nrow = 5*K+3)))
      
      b0_new <- Theta_new[1]
      N_new <- exp(Theta_new[2:(K+1)])
      R_eff_new <- exp(Theta_new[(K+2):(2*K+1)])
      n_new <- exp(Theta_new[(2*K+2):(3*K+1)])
      mu_new <- Theta_new[(3*K+2):(4*K+2)]
      sigma_new <- exp(Theta_new[(4*K+3):(5*K+3)])
      
      r <- log_post_mod_np(S, grid, b0_new, c, N_new, R_eff_new, e, n_new, theta, mu_new, sigma_new, Y_obs, p, prior, cf_error = cf_error) -
        log_post_mod_np(S, grid, b0_old, c, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, cf_error = cf_error)
      
      if(log(runif(1)) < r){
        res[i+1, ] <- c(b0_new, log(N_new), 
                        log(R_eff_new), log(n_new), mu_new, log(sigma_new))
      }
      else{
        res[i+1,] <- res[i,]
      }
      pb$tick()
    }
  }
  res[,2:(K+1)] <- exp(res[,2:(K+1)])
  res[,(K+2):(2*K+1)] <- exp(res[,(K+2):(2*K+1)])
  res[,(2*K+2):(3*K+1)] <- exp(res[,(2*K+2):(3*K+1)])
  res[,(4*K+3):(5*K+3)] <- exp(res[,(4*K+3):(5*K+3)])
  return(res)
}

fit_mod_np_gal <- function(S, grid, Theta, Ng, UDG_ID, c, e, theta, Y, p, prior, M, tune, n = 1000, gamma = 0.05, prob_model = TRUE, cf_error){
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta || Current MCMC iteration: :current]",
                         total = M,
                         complete = "=",
                         current = ">",
                         clear = FALSE,
                         width = 200) 
  K <- nrow(c)
  res <- matrix(0, nrow = M+1, ncol = 5*K+3)
  col_names <- c('b0', paste0('N_gal_', 1:Ng), paste0('N_', UDG_ID), paste0('R_gal_', 1:Ng), paste0('R_', UDG_ID), paste0('n_gal_', 1:Ng), paste0('n_', UDG_ID), 'mu_0', paste0('mu_gal_', 1:Ng), paste0('mu_', UDG_ID), 'sigma_0', paste0('sigma_gal_', 1:Ng), paste0('sigma_', UDG_ID))
  colnames(res) <- col_names
  res[1,] <- c(Theta$b0, Theta$N, Theta$R_eff, Theta$n, Theta$mu, Theta$sigma)
  t_b <- tune$b0
  t_Ng <- tune$Ng
  t_Nu <- tune$Nu
  t_R <- tune$R
  t_n <- tune$n
  t_m <- tune$mu
  t_s <- tune$sigma
  gamma <- 2.38^2/(5*K+3)
  if(prob_model) {
    Y_dat <- Y[, 1:7]
    Y_prob <- Y[, 8:507]
  }
  idx <- sample(1:500, M, replace = T)
  for (i in 1:M) {
    if(i <= n){
      if(prob_model){
        Y_obs <- Y_dat[rbernoulli(nrow(Y_dat), p = Y_prob[,idx[i]]),]
      }
      else{
        Y_obs <- Y
      }
      b0_old <- res[i, 1]
      N_old_g <- exp(res[i, 2:(Ng+1)])
      N_old_u <- exp(res[i, (Ng+2):(K+1)])
      N_old <- c(N_old_g, N_old_u)
      R_eff_old <- exp(res[i,(K+2):(2*K+1)])
      n_old <- exp(res[i,(2*K+2):(3*K+1)])
      mu_old <- res[i, (3*K+2):(4*K+2)]
      sigma_old <- exp(res[i, (4*K+3):(5*K+3)])
      
      b0_new <- rnorm(1, b0_old, t_b)
      N_new_g <- exp(rnorm(Ng, log(N_old_g), t_Ng))
      N_new_u <- exp(rnorm(K-Ng, res[i, (Ng+2):(K+1)], t_Nu))
      N_new <- c(N_new_g, N_new_u)
      R_eff_new <- exp(rnorm(K, log(R_eff_old), t_R))
      n_new <- exp(rnorm(K, log(n_old), t_n))
      mu_new <- rnorm(K+1, mu_old, t_m)
      sigma_new <- exp(rnorm(K+1, log(sigma_old), t_s))
      
      r <- log_post_mod_np_gal(S, grid, b0_new, c, N_new, R_eff_new, e, n_new, theta, mu_new, sigma_new, Y_obs, p, prior, Ng, cf_error = cf_error)-
        log_post_mod_np_gal(S, grid, b0_old, c, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, Ng, cf_error = cf_error)
      
      C <- cov(res[c(1:n),])
      
      Mu <- colMeans(res[c(1:n),])
      
      if(log(runif(1)) < r){
        res[i+1, ] <- c(b0_new, log(N_new_g), log(N_new_u), 
                        log(R_eff_new), log(n_new), mu_new, log(sigma_new))
      }
      else{
        res[i+1,] <- res[i,]
      }
      pb$tick()
    }
    else{
      if(prob_model){
        Y_obs <- Y_dat[rbernoulli(nrow(Y_dat), p = Y_prob[,idx[i]]),]
      }
      else{
        Y_obs <- Y
      }
      b0_old <- res[i, 1]
      N_old_g <- exp(res[i, 2:(Ng+1)])
      N_old_u <- exp(res[i, (Ng+2):(K+1)])
      N_old <- c(N_old_g, N_old_u)
      R_eff_old <- exp(res[i,(K+2):(2*K+1)])
      n_old <- exp(res[i,(2*K+2):(3*K+1)])
      mu_old <- res[i, (3*K+2):(4*K+2)]
      sigma_old <- exp(res[i, (4*K+3):(5*K+3)])
      
      X <- res[i,]
      C <- (i-2)/(i-1)*C + (X-Mu)%*%t(X-Mu)/i
      Mu <- (Mu*(i-1) + X)/i
      
      Theta_new <- MASS::mvrnorm(1, mu = X, 
                                 Sigma = gamma*(C + 1e-7*diag(1, nrow = 5*K+3)))
      
      b0_new <- Theta_new[1]
      N_new_g <- exp(Theta_new[2:(Ng+1)])
      N_new_u <- exp(Theta_new[(Ng+2):(K+1)])
      N_new <- c(N_new_g, N_new_u)
      R_eff_new <- exp(Theta_new[(K+2):(2*K+1)])
      n_new <- exp(Theta_new[(2*K+2):(3*K+1)])
      mu_new <- Theta_new[(3*K+2):(4*K+2)]
      sigma_new <- exp(Theta_new[(4*K+3):(5*K+3)])
      
      r <- log_post_mod_np_gal(S, grid, b0_new, c, N_new, R_eff_new, e, n_new, theta, mu_new, sigma_new, Y_obs, p, prior, Ng, cf_error = cf_error) -
        log_post_mod_np_gal(S, grid, b0_old, c, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, Ng, cf_error = cf_error)
      
      if(log(runif(1)) < r){
        res[i+1, ] <- c(b0_new, log(N_new_g), log(N_new_u), 
                        log(R_eff_new), log(n_new), mu_new, log(sigma_new))
      }
      else{
        res[i+1,] <- res[i,]
      }
      pb$tick()
    }
  }
  res[,2:(Ng+1)] <- exp(res[,2:(Ng+1)])
  res[,(Ng+2):(K+1)] <- exp(res[,(Ng+2):(K+1)])
  res[,(K+2):(2*K+1)] <- exp(res[,(K+2):(2*K+1)])
  res[,(2*K+2):(3*K+1)] <- exp(res[,(2*K+2):(3*K+1)])
  res[,(4*K+3):(5*K+3)] <- exp(res[,(4*K+3):(5*K+3)])
  return(res)
}

fit_mod_np_aMHwG <- function(S, grid, Theta, UDG_ID, c, e, theta, Y, p, prior, M, tune, n = 1000, prob_model = TRUE, cam = 'ACS'){
  K <- nrow(c)
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta || Current MCMC iteration: :current]",
                         total = M*(K+1)+n,
                         complete = "=",
                         incomplete = "-",
                         current = ">",
                         clear = FALSE,
                         width = 200) 
  
  res <- matrix(0, nrow = n + M*(K+1) + 1, ncol = 5*K+3)
  col_names <- c('b0', paste0('N_', UDG_ID), paste0('R_', UDG_ID), paste0('n_', UDG_ID), 'mu_0', paste0('mu_', UDG_ID), 'sigma_0', paste0('sigma_', UDG_ID))
  colnames(res) <- col_names
  res[1,] <- c(Theta$b0, Theta$N, Theta$R_eff, Theta$n, Theta$mu, Theta$sigma)
  t_b <- tune$b0
  t_N <- tune$N
  t_R <- tune$R
  t_n <- tune$n
  t_m <- tune$mu
  t_s <- tune$sigma
  i <- 1
  while(i <= M*(K+1)+n) {
    if(i <= n){
      if(prob_model){
        Y_obs <- Y[rbernoulli(nrow(Y), p = Y$GC_prob),]
      }
      else{
        Y_obs <- Y
      }
      b0_old <- res[i, 1]
      N_old <- exp(res[i, 2:(K+1)])
      R_eff_old <- exp(res[i,(K+2):(2*K+1)])
      n_old <- exp(res[i,(2*K+2):(3*K+1)])
      mu_old <- res[i, (3*K+2):(4*K+2)]
      sigma_old <- res[i, (4*K+3):(5*K+3)]
      
      b0_new <- rnorm(1, b0_old, t_b)
      N_new <- exp(rnorm(K, log(N_old), t_N))
      R_eff_new <- exp(rnorm(K, log(R_eff_old), t_R))
      n_new <- exp(rnorm(K, log(n_old), t_n))
      mu_new <- rnorm(K+1, mu_old, t_m)
      sigma_new <- rnorm(K+1, sigma_old, t_s)
      
      r <- log_post_mod_np(S, grid, b0_new, c, N_new, R_eff_new, e, n_new, theta, mu_new, sigma_new, Y_obs, p, prior, cam = cam) - 
        log_post_mod_np(S, grid, b0_old, c, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, cam = cam)
      
      C <- cov(res[c(1:n),])
      
      Mu <- colMeans(res[c(1:n),])
      
      if(log(runif(1)) < r){
        res[i+1, ] <- c(b0_new, log(N_new), 
                        log(R_eff_new), log(n_new), mu_new, sigma_new)
        i <- i+1
      }
      else{
        res[i+1,] <- res[i,]
        i <- i+1
      }
      pb$tick()
    }
    else{
      if(prob_model){
        Y_obs <- Y[rbernoulli(nrow(Y), p = Y$GC_prob),]
      }
      else{
        Y_obs <- Y
      }
      for(j in 1:(K+1)){
        if(j == 1){
          b0_old <- res[i, 1]
          N_old <- exp(res[i, 2:(K+1)])
          R_eff_old <- exp(res[i, (K+2):(2*K+1)])
          n_old <- exp(res[i, (2*K+2):(3*K+1)])
          mu_old <- res[i, (3*K+2):(4*K+2)]
          sigma_old <- res[i, (4*K+3):(5*K+3)]
          
          X <- res[i,]
          C <- (i-2)/(i-1)*C + (X-Mu)%*%t(X-Mu)/(i)
          Mu <- (Mu*(i-1) + X)/(i)
          
          gamma <- 2.38^2/3
          Theta_new <- MASS::mvrnorm(1, mu = X, 
                                     Sigma = gamma*(C + diag(1e-7, nrow = 5*K+3)))
          
          b0_new <- Theta_new[1]
          mu_new <- c(Theta_new[3*K+2], mu_old[-1])
          sigma_new <- c(Theta_new[4*K+3], sigma_old[-1])
          
          r <- log_post_mod_np(S, grid, b0_new, c, N_old, R_eff_old, e, n_old, theta, mu_new, sigma_new, Y_obs, p, prior, cam = cam) -
            log_post_mod_np(S, grid, b0_old, c, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, cam = cam)
          
          if(log(runif(1)) < r){
            res[i+1, ] <- c(b0_new, log(N_old), 
                            log(R_eff_old), log(n_old), mu_new, sigma_new)
            i <- i+1
          }
          else{
            res[i+1,] <- res[i,]
            i <- i+1
          } 
          pb$tick()
        }
        else{
          b0_old <- res[i, 1]
          N_old <- N_new <- exp(res[i, 2:(K+1)])
          R_eff_old <- R_eff_new <- exp(res[i, (K+2):(2*K+1)])
          n_old <- n_new <- exp(res[i, (2*K+2):(3*K+1)])
          mu_old <- mu_new <- res[i, (3*K+2):(4*K+2)]
          sigma_old <- sigma_new <- res[i, (4*K+3):(5*K+3)]
          
          X <- res[i,]
          C <- (i-2)/(i-1)*C + (X-Mu)%*%t(X-Mu)/(i)
          Mu <- (Mu*(i-1) + X)/(i)
          
          gamma <- 2.38^2/5
          Theta_new <- MASS::mvrnorm(1, mu = X, 
                                     Sigma = gamma*(C + diag(1e-7, nrow = 5*K+3)))
          
          N_new[j-1] <- exp(Theta_new[j])
          R_eff_new[j-1] <- exp(Theta_new[K+j])
          n_new[j-1] <- exp(Theta_new[2*K+j])
          mu_new[j] <- Theta_new[3*K+j+1]
          sigma_new[j] <- Theta_new[4*K+j+2]
          
          r <- log_post_mod_np(S, grid, b0_old, c, N_new, R_eff_new, e, n_new, theta, mu_new, sigma_new, Y_obs, p, prior, cam = cam) -
            log_post_mod_np(S, grid, b0_old, c, N_old, R_eff_old, e, n_old, theta, mu_old, sigma_old, Y_obs, p, prior, cam = cam)
          
          if(log(runif(1)) < r){
            res[i+1, ] <- c(b0_old, log(N_new), 
                            log(R_eff_new), log(n_new), mu_new, sigma_new)
            i <- i+1
          }
          else{
            res[i+1,] <- res[i,]
            i <- i+1
          } 
          pb$tick()
        }
      }
    }
  }
  res[,2:(K+1)] <- exp(res[,2:(K+1)])
  res[,(K+2):(2*K+1)] <- exp(res[,(K+2):(2*K+1)])
  res[,(2*K+2):(3*K+1)] <- exp(res[,(2*K+2):(3*K+1)])
  return(res)
}








