#####################################################
### Extra Simulations for K = 3            ##########
### Found in the Appendix  Section B.1     ##########
#####################################################

## Load relevant packages
library(MASS)
library(DirichletReg)
library(devtools)
library(future.apply)

# Install BFMMM package from github
# If already installed, do not run the following line
install_github('ndmarco/BayesFMMM')

library(BayesFMMM)


run_sim3 <- function(iter){
  set.seed(iter)
  n_obs_vec <- c(50, 100, 200)
  for(n in 1:length(n_obs_vec)){
    dir.create(paste0(n_obs_vec[n],"_obs"))
    dir.create(paste0(n_obs_vec[n],"_obs/sim",iter))
    n_obs <- n_obs_vec[n]
    Y <- readRDS(system.file("test-data", "Sim_data.RDS", package = "BayesFMMM"))
    time <- readRDS(system.file("test-data", "time.RDS", package = "BayesFMMM"))

    ## Set Hyperparameters
    tot_mcmc_iters <- 150
    n_try <- 1
    k <- 2
    n_funct <- 40
    basis_degree <- 3
    n_eigen <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)

    ## Run function to get B splines
    x <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                 basis_degree, n_eigen, boundary_knots,
                                 internal_knots)
    B <- x$B[[1]]

    # Specify parameters for generating functions
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

    nu[1,] <- mvrnorm(n=1, mu = seq(6,-8, -2), Sigma = 2*p)
    nu[2,] <- mvrnorm(n=1, mu = seq(-8, 6, 2), Sigma = 2*p)
    nu[3,] <- mvrnorm(n=1, mu = rep(4, 8), Sigma = 2*p)
    decomp <- svd(nu, nv = 8)
    Phi_i <- matrix(0, nrow = 3, ncol = 8)
    Phi_i[1,] <- t(rnorm(5, 0, 1)) %*% t(decomp$v[,4:8])
    Phi_i[2,] <- t(rnorm(5, 0, 1)) %*% t(decomp$v[,4:8])
    Phi_i[3,] <- t(rnorm(5, 0, 1)) %*% t(decomp$v[,4:8])
    Phi_1 <- Phi_i
    Phi_i[1,] <- t(rnorm(5, 0, .9)) %*% t(decomp$v[,4:8])
    Phi_i[2,] <- t(rnorm(5, 0, .9)) %*% t(decomp$v[,4:8])
    Phi_i[3,] <- t(rnorm(5, 0, .9)) %*% t(decomp$v[,4:8])
    Phi_2 <- Phi_i
    Phi_i[1,] <- t(rnorm(5, 0, 0.7)) %*% t(decomp$v[,4:8])
    Phi_i[2,] <- t(rnorm(5, 0, 0.7)) %*% t(decomp$v[,4:8])
    Phi_i[3,] <- t(rnorm(5, 0, 0.7)) %*% t(decomp$v[,4:8])
    Phi_3 <- Phi_i
    Phi <- array(0, dim = c(3, 8, 3))
    Phi[,,1] <- matrix(Phi_1, nrow = 2)
    Phi[,,2] <- matrix(Phi_2, nrow = 2)
    Phi[,,3] <- matrix(Phi_3, nrow = 2)

    chi <- matrix(rnorm(n_obs *3, 0, 1), ncol = 3, nrow=n_obs)

    Z <- matrix(0, nrow = n_obs, ncol = 3)
    alpha <- c(100, 1, 1)
    for(i in 1:(n_obs * 0.2)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 100, 1)
    for(i in (n_obs * 0.2 + 1):(n_obs * 0.4)){
      Z[i,] <- rdirichlet(1, alpha)
    }
    alpha <- c(1, 1, 100)
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

    # Generate the observed sample paths
    for(i in 1:n_obs){
      mean = rep(0,100)
      for(j in 1:3){
        mean = mean + Z[i,j] * B %*% nu[j,]
        for(m in 1:3){
          mean = mean + Z[i,j] * chi[i,m] * B %*% Phi[j, ,m]
        }
      }
      y[[i]] = mvrnorm(n = 1, mean, diag(0.001, 100))
    }

    x <- list("y" = y, "nu" = nu, "Z" = Z, "Phi" = Phi, "Chi" = chi)

    saveRDS(x, paste("./", n_obs_vec[n],"_obs/sim",iter, "/truth.RDS", sep = ""))


    ## Set Hyperparameters
    tot_mcmc_iters <- 2000
    n_try <- 50
    k <- 3
    n_funct <- n_obs
    basis_degree <- 3
    n_eigen <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    Y <-y

    ## Get Estimates of Z and nu
    est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                    basis_degree, n_eigen, boundary_knots,
                                    internal_knots, b = 10000)
    tot_mcmc_iters <- 20000
    n_try <- 10
    ## Get estimates of other parameters
    est2 <- BFMMM_Theta_est(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                            basis_degree, n_eigen, boundary_knots,
                            internal_knots, est1)
    dir_i <- paste("./", n_obs_vec[n],"_obs/sim",iter, "/", sep="")
    tot_mcmc_iters <- 500000
    MCMC.chain <-BFMMM_warm_start(tot_mcmc_iters, k, Y, time, n_funct,
                                  basis_degree, n_eigen, boundary_knots,
                                  internal_knots, est1, est2, dir = dir_i,
                                  thinning_num = 100, r_stored_iters = 10000)
  }
}


##### Run Simulation

### Set working dir (different than the main simulation study)
setwd("")

ncpu <- min(5, availableCores())
#
plan(multisession, workers = ncpu)

already_ran <- dir(paste0(getwd(), "/200_obs"))
to_run <- which(!paste0("sim", 1:10) %in% already_ran)
seeds <- to_run
start_time <- Sys.time()
future_lapply(seeds, function(this_seed) run_sim3(this_seed))
end_time <- Sys.time()

# get elapsed time of the simulation
total_time <- end_time - start_time
