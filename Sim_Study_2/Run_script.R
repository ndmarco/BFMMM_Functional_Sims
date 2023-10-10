library(devtools)
# Install BFMMM package from github
# If already installed, do not run the following line
install_github('ndmarco/BayesFMMM')

library(BayesFMMM)
library(MASS)
library(DirichletReg)
library(future.apply)

run_sim <- function(iter){
  set.seed(iter)
  ## Load sample data
  Y <- readRDS(system.file("test-data", "Sim_data.RDS", package = "BayesFMMM"))
  time <- readRDS(system.file("test-data", "time.RDS", package = "BayesFMMM"))

  ## Set Hyperparameters
  tot_mcmc_iters <- 150
  n_try <- 1
  k <- 3
  n_funct <- 40
  basis_degree <- 3
  n_eigen <- 3
  boundary_knots <- c(0, 1000)
  internal_knots <- c(200, 400, 600, 800)

  ## Run function
  x <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                               basis_degree, n_eigen, boundary_knots,
                               internal_knots)
  B <- x$B[[1]]
  n_obs = 200
  nu <- matrix(0,nrow =3, ncol = 8)
  p <- matrix(0, 8, 8)
  for(i in 1:8){
    p[1,1] = 1
    if(i > 1){
      p[i,i] = 2
      p[i-1,i] = -1
      p[i,i-1] = -1
    }
  }
  nu[1,] <- mvrnorm(n=1, mu = seq(6,-8, -2), Sigma = 4*p)
  nu[2,] <- mvrnorm(n=1, mu = seq(-8, 6, 2), Sigma = 4*p)
  nu[3,] <- mvrnorm(n=1, mu = rep(0,8), Sigma = 4*p)
  Phi_1 <- rnorm(24, 0, 1)
  Phi_2 <- rnorm(24, 0, 0.5)
  Phi_3 <- rnorm(24, 0, 0.2)
  Phi <- array(0, dim = c(3, 8, 3))
  Phi[,,1] <- matrix(Phi_1, nrow = 3)
  Phi[,,2] <- matrix(Phi_2, nrow = 3)
  Phi[,,3] <- matrix(Phi_3, nrow = 3)


  chi <- matrix(rnorm(n_obs *3, 0, 1), ncol = 3, nrow=n_obs)

  Z <- matrix(0, nrow = n_obs, ncol = 3)
  alpha <- c(30, 1, 1)
  for(i in 1:(n_obs * 0.2)){
    Z[i,] <- rdirichlet(1, alpha)
  }
  alpha <- c(1, 30, 1)
  for(i in (n_obs * 0.2 + 1):(n_obs * 0.4)){
    Z[i,] <- rdirichlet(1, alpha)
  }
  alpha <- c(1, 1, 30)
  for(i in (n_obs * 0.4 + 1):(n_obs * 0.6)){
    Z[i,] <- rdirichlet(1, alpha)
  }
  alpha <- c(1, 1, 1)
  for(i in (n_obs * 0.6 + 1):n_obs){
    Z[i,] <- rdirichlet(1, alpha)
  }

  y <- rep(0,100)
  y <- rep(list(y), n_obs)
  time <- time[[1]]
  time <- rep(list(time), n_obs)
  for(i in 1:n_obs){
    mean = rep(0,10)
    for(j in 1:3){
      mean = mean + Z[i,j] * B %*% nu[j,]
      for(m in 1:3){
        mean = mean + Z[i,j] * chi[i,m] * B %*% Phi[j, ,m]
      }
    }
    y[[i]] = mvrnorm(n = 1, mean, diag(0.001, 100))
  }
  dir.create("data")
  x <- list("y" = y, "nu" = nu, "Z" = Z, "Phi" = Phi, "Chi" = chi)
  saveRDS(x, paste("./data/data", iter, ".RDS", sep = ""))
  Y <- y


  ####### 2 Clusters #########
  ############################
  ## Set Hyperparameters
  tot_mcmc_iters <- 2000
  n_try <- 50
  k <- 2
  n_funct <- 200
  basis_degree <- 3
  n_eigen <- 3
  boundary_knots <- c(0, 1000)
  internal_knots <- c(200, 400, 600, 800)
  dir.create(paste0("2_clusters/trace", iter))

  ## Get Estimates of Z and nu
  est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots)
  tot_mcmc_iters <- 4000
  n_try <- 5
  ## Get estimates of other parameters
  est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                          basis_degree, n_eigen, boundary_knots,
                          internal_knots, est1)
  dir_i <- paste( "./2_clusters/trace", iter, "/", sep="")
  tot_mcmc_iters <- 100000
  MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                basis_degree, n_eigen, boundary_knots,
                                internal_knots, est1, est2, dir = dir_i,
                                thinning_num = 10, r_stored_iters = 10000)

  ####### 3 Clusters #########
  ############################
  ## Set Hyperparameters
  tot_mcmc_iters <- 2000
  n_try <- 50
  k <- 3
  n_funct <- 200
  basis_degree <- 3
  n_eigen <- 3
  boundary_knots <- c(0, 1000)
  internal_knots <- c(200, 400, 600, 800)
  dir.create(paste0("3_clusters/trace", iter))

  ## Get Estimates of Z and nu
  est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots)
  tot_mcmc_iters <- 4000
  n_try <- 5
  ## Get estimates of other parameters
  est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                          basis_degree, n_eigen, boundary_knots,
                          internal_knots, est1)
  dir_i <- paste("./3_clusters/trace", iter, "/", sep="")
  tot_mcmc_iters <- 100000
  MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                basis_degree, n_eigen, boundary_knots,
                                internal_knots, est1, est2, dir = dir_i,
                                thinning_num = 10, r_stored_iters = 10000)

  ####### 4 Clusters #########
  ############################
  ## Set Hyperparameters
  tot_mcmc_iters <- 2000
  n_try <- 50
  k <- 4
  n_funct <- 200
  basis_degree <- 3
  n_eigen <- 3
  boundary_knots <- c(0, 1000)
  internal_knots <- c(200, 400, 600, 800)
  dir.create(paste0("4_clusters/trace", iter))

  ## Get Estimates of Z and nu
  est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots)
  tot_mcmc_iters <- 4000
  n_try <- 5
  ## Get estimates of other parameters
  est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                          basis_degree, n_eigen, boundary_knots,
                          internal_knots, est1)
  dir_i <- paste("./4_clusters/trace", iter, "/", sep="")
  tot_mcmc_iters <- 100000
  MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                basis_degree, n_eigen, boundary_knots,
                                internal_knots, est1, est2, dir = dir_i,
                                thinning_num = 10, r_stored_iters = 10000)

  ####### 5 Clusters #########
  ############################
  ## Set Hyperparameters
  tot_mcmc_iters <- 2000
  n_try <- 50
  k <- 5
  n_funct <- 200
  basis_degree <- 3
  n_eigen <- 3
  boundary_knots <- c(0, 1000)
  internal_knots <- c(200, 400, 600, 800)
  dir.create(paste0("5_clusters/trace", iter))

  ## Get Estimates of Z and nu
  est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots)
  tot_mcmc_iters <- 4000
  n_try <- 5
  ## Get estimates of other parameters
  est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                          basis_degree, n_eigen, boundary_knots,
                          internal_knots, est1)
  dir_i <- paste("./5_clusters/trace", iter, "/", sep="")
  tot_mcmc_iters <- 100000
  MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                basis_degree, n_eigen, boundary_knots,
                                internal_knots, est1, est2, dir = dir_i,
                                thinning_num = 10, r_stored_iters = 10000)
}

############################################
## Before running, change the directories ##
############################################

setwd()

##File Structure should be as follows:

## 2_clusters
######## trace1
######## ...
######## trace10
## 3_clusters
######## trace1
######## ...
######## trace10
## 4_clusters
######## trace1
######## ...
######## trace10
## 5_clusters
######## trace1
######## ...
######## trace10
## data
dir.create("2_clusters")
dir.create("3_clusters")
dir.create("4_clusters")
dir.create("5_clusters")
dir.create("data")

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)

already_ran <- dir(paste0(getwd(), "/5_clusters/"))
to_run <- which(!paste0("trace", 1:10) %in% already_ran)
seeds <- to_run
start_time <- Sys.time()
future_lapply(seeds, function(this_seed) run_sim(this_seed))
end_time <- Sys.time()

# get elapsed time of the simulation
total_time <- end_time - start_time
