library(BayesFMMM)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(pracma)
library(reticulate)

###############################################################################
######## 3 Features (Found in Appendix)         ###############################
###############################################################################

### Set working dir for 3 features
setwd()

Y <- readRDS(system.file("test-data", "Sim_data.RDS", package = "BayesFMMM"))
time <- readRDS(system.file("test-data", "time.RDS", package = "BayesFMMM"))
tot_mcmc_iters <- 150
n_try <- 1
k <- 2
n_funct <- 40
basis_degree <- 3
n_eigen <- 3
boundary_knots <- c(0, 1000)
internal_knots <- c(200, 400, 600, 800)

est1 <- BFMMM_Nu_Z_multiple_try(tot_mcmc_iters, n_try, k, Y, time, n_funct,
                                basis_degree, n_eigen, boundary_knots,
                                internal_knots)

## Get basis functions to reconstruct the true functions
## Basis functions are B-splines of degree 3 with the boundary knots and internal knots specified above
B <- est1$B

# Initialize placeholders for estimates of the mean square integrated error of mean and covariance functions
err_Z <- matrix(0, 10, 3)
int_err_mean1 <- matrix(0, 10, 3)
int_err_mean2 <- matrix(0, 10, 3)
int_err_mean3 <- matrix(0, 10, 3)
int_err_cov12 <- matrix(0, 10, 3)
int_err_cov1 <- matrix(0, 10, 3)
int_err_cov2 <- matrix(0, 10, 3)
int_err_cov3 <- matrix(0, 10, 3)
int_err_cov13 <- matrix(0, 10, 3)
int_err_cov23 <- matrix(0, 10, 3)

# Initialize placeholders for estimates of the squared norm of mean and covariance functions
norm_mu1 <- matrix(1, nrow = 10, ncol = 3)
norm_mu2 <- matrix(1, nrow = 10, ncol = 3)
norm_mu3 <- matrix(1, nrow = 10, ncol = 3)
norm_C1 <- matrix(1, nrow = 10, ncol = 3)
norm_C2 <- matrix(1, nrow = 10, ncol = 3)
norm_C3 <- matrix(1, nrow = 10, ncol = 3)
norm_C12 <- matrix(1, nrow = 10, ncol = 3)
norm_C13 <- matrix(1, nrow = 10, ncol = 3)
norm_C23 <- matrix(1, nrow = 10, ncol = 3)

tern <- function(Z) {
  X <- matrix(0, nrow = nrow(Z), ncol = 2)
  X[,1] <- (Z[,2] + 0.5* Z[,3])
  X[,2] <- sqrt(0.75) * Z[,3]
  return(X)
}

anti_tern <- function(tZ){
  X <- matrix(0, nrow = nrow(tZ), ncol = 3)
  X[,3] <- tZ[,2] / sqrt(0.75)
  X[,2] <- tZ[,1] - 0.5 * X[,3]
  X[,1] <- 1 - (X[,3] + X[,2])
  return(X)
}

### May have to change path to python3
Sys.setenv(RETICULATE_PYTHON = "/usr/local/bin/python3")
use_python("/usr/local/bin/python3")
cv2 <- import("cv2")
np <- import("numpy")

rescale_function <- function(Z){
  # Reduce data to two dimensions
  tZ <- tern(Z)


  #### Python functions to find triangle that encloses the points with the smallest area
  tz <- np_array(tZ)
  tz = tz$astype(np$float32)
  tri_py <- cv2$minEnclosingTriangle(tz)

  # Define points for triangle in two dimensions
  triangle_t <- matrix(0, ncol = 2, nrow = 3)
  triangle_t[1,] <- tri_py[[2]][1,1,]
  triangle_t[2,] <- tri_py[[2]][2,1,]
  triangle_t[3,] <- tri_py[[2]][3,1,]

  # Define points for triangle in original space
  triangle <- anti_tern(triangle_t)

  #
  #
  # # Find distance from remaining points to the midpoint between the first two points of the triangle
  # # ph1 <- dist(rbind((tZ[hpts[ind[1,1]],] + tZ[hpts[ind[1,1]],])/2, tZ[hpts[-c(ind[1,1], ind[1,2])],]))
  # # ph1 <- as.matrix(ph1)
  # ind1 <- which.max(ph[ind[1,1],] + ph[ind[1,2],] )
  #
  # triangle[3,] <- Z[hpts[ind1],]
  # triangle_t[3,] <- tZ[hpts[ind1],]



  ### Plot graphically in 2 dimensions
  #######################################################
  # plot(tZ, cex = 0.5)
  # lines(rbind(triangle_t, triangle_t[1,]), col="green")
  #######################################################


  # Reorient data to the original triangle
  inv_trans <- solve(triangle, tol = 1e-16)
  trans_Z <- Z %*% inv_trans

  # New data may not lie on the simplex
  # Shrink data to fit in triangle
  diag_mat <- diag(1,3)
  l <- 1

  while(sum(rowSums(trans_Z < 0) > 0) > 0){
    if(sum(rowSums(trans_Z < 0) > 0) > 1){
      shrink_factor <- colSums(trans_Z[rowSums(trans_Z < 0) > 0,])
    }else{
      shrink_factor <- trans_Z[rowSums(trans_Z < 0) > 0,]
    }

    diag_mat <- diag_mat %*% matrix(c(1 - (shrink_factor[1] * 0.001), 0.5 * (shrink_factor[1] * 0.001), 0.5 *
                                        (shrink_factor[1] * 0.001), 0.5*  (shrink_factor[2] * 0.001),
                                      1- (shrink_factor[2] * 0.001), 0.5 * (shrink_factor[2] * 0.001),
                                      0.5 * (shrink_factor[3] * 0.001), 0.5 * (shrink_factor[3] * 0.001),
                                      1-(shrink_factor[3] * 0.001)), ncol = 3, nrow=3, byrow = T)
    inv_trans <- solve(triangle, diag_mat)
    trans_Z <- Z %*% inv_trans
    l <- l + 1
  }

  trans_mat <- solve(inv_trans)

  return(list("t_mat" = trans_mat, "rescaled_Z" = trans_Z))
}

# Estimate the R-MISE for each of the 3 sample sizes (N = 40, N = 80, N = 160)
for(j in 1:3){
  if(j == 1){
    dir <- "./50_obs/"
  }
  if(j == 2){
    dir <- "./100_obs/"
  }
  if(j == 3){
    dir <- "./200_obs/"
  }
  for(i in 1:10){
    # Read in true parameters
    x <- readRDS(paste(dir, "sim", i, "/truth.RDS", sep=""))
    x$Phi_true <- x$Phi
    x$nu_true <- x$nu
    Z_true <- x$Z

    # Get true means
    nu_1_true <-  B[[1]] %*% t(t(x$nu_true[1,]))
    nu_2_true <-  B[[1]] %*% t(t(x$nu_true[2,]))
    nu_3_true <- B[[1]] %*% t(t(x$nu_true[3,]))

    # Get true covariance and cross-covariance functions
    cov1_true <- matrix(0, 100, 100)
    for(q in 1:100){
      for(k in 1:100){
        for(m in 1:3){
          cov1_true[q,k] <- cov1_true[q,k] + x$Phi_true[1, ,m] %*% t(t(B[[1]][q,])) %*% B[[1]][k,] %*% t(t(x$Phi_true[1, ,m]))
        }
      }
    }

    cov2_true <- matrix(0, 100, 100)
    for(q in 1:100){
      for(k in 1:100){
        for(m in 1:3){
          cov2_true[q,k] <- cov2_true[q,k] + x$Phi_true[2, ,m] %*% t(t(B[[1]][q,])) %*% B[[1]][k,] %*% t(t(x$Phi_true[2, ,m]))
        }
      }
    }

    cov3_true <- matrix(0, 100, 100)
    for(q in 1:100){
      for(k in 1:100){
        for(m in 1:3){
          cov3_true[q,k] <- cov3_true[q,k] + x$Phi_true[3, ,m] %*% t(t(B[[1]][q,])) %*% B[[1]][k,] %*% t(t(x$Phi_true[3, ,m]))
        }
      }
    }

    cov12_true <- matrix(0, 100, 100)
    for(q in 1:100){
      for(k in 1:100){
        for(m in 1:3){
          cov12_true[q,k] <- cov12_true[q,k] + x$Phi_true[1, ,m] %*% t(t(B[[1]][q,])) %*% B[[1]][k,] %*% t(t(x$Phi_true[2, ,m]))
        }
      }
    }
    cov13_true <- matrix(0, 100, 100)
    for(q in 1:100){
      for(k in 1:100){
        for(m in 1:3){
          cov13_true[q,k] <- cov13_true[q,k] + x$Phi_true[1, ,m] %*% t(t(B[[1]][q,])) %*% B[[1]][k,] %*% t(t(x$Phi_true[3, ,m]))
        }
      }
    }
    cov23_true <- matrix(0, 100, 100)
    for(q in 1:100){
      for(k in 1:100){
        for(m in 1:3){
          cov23_true[q,k] <- cov23_true[q,k] + x$Phi_true[2, ,m] %*% t(t(B[[1]][q,])) %*% B[[1]][k,] %*% t(t(x$Phi_true[3, ,m]))
        }
      }
    }

    # Hyperparameters for the sampling already conducted
    n_files <- 50
    basis_degree <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    dir_i <- paste(dir, "sim", i, "/", sep = "")

    Z_est <- ZCI(dir_i, n_files)

    Z_rescaled <- Z_est$Z_trace
    trans_mat <- matrix(NA, nrow = 3 * length(Z_est$Z_trace[1,1,]), ncol = 3)
    for(iter in 1:4500){
      obj <- try(rescale_function(Z_est$Z_trace[,,iter]), silent = T)
      if(!inherits(obj, "try-error")){
        Z_rescaled[,,iter] <- obj$rescaled_Z
        trans_mat[((iter-1)*3 + 1):(iter*3), ] = obj$t_mat
        # potential label switching
        if(iter > 1){
          mat <- solve(t(Z_rescaled[,,iter]) %*% Z_rescaled[,,iter]) %*% t(Z_rescaled[,,iter]) %*% Z_rescaled[,,iter-1]
          mat <- round(mat)
          Z_rescaled[,,iter] <- Z_rescaled[,,iter] %*% mat
          trans_mat[((iter-1)*3 + 1):(iter*3), ] = solve(solve(obj$t_mat) %*% mat)
        }

      }else{
        Z_rescaled[,,iter] <- Z_rescaled[,,iter-1]
        trans_mat[((iter-1)*3 + 1):(iter*3), ] = trans_mat[((iter-2)*3 + 1):((iter-1)*3), ]
      }
    }

    Z_median <- matrix(0, nrow = nrow(Z_rescaled), ncol = ncol(Z_rescaled))
    for(iter in 1:nrow(Z_rescaled)){
      for(j_iter in 1:ncol(Z_rescaled)){
        Z_median[iter,j_iter] = median(Z_rescaled[iter,j_iter,])
      }
    }

    ### Account for label Switching
    mat <- solve(t(Z_median) %*% Z_median) %*% t(Z_median) %*% Z_true
    mat <- round(mat)
    indices <- rep(NA, 3)
    for(r in 1:3){
      indices[r] = which(mat[,r] > 0)
    }

    Z_median <- Z_median %*% mat


    # Get estimated mean, covariance, and cross-covariance functions
    nu_1 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, indices[1], burnin_prop = 0.7, trans_mats = trans_mat)
    nu_2 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, indices[2], burnin_prop = 0.7, trans_mats = trans_mat)
    nu_3 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, indices[3], burnin_prop = 0.7, trans_mats = trans_mat)
    cov1 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, indices[1], indices[1], burnin_prop = 0.7, trans_mats = trans_mat)
    cov2 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, indices[2], indices[2], burnin_prop = 0.7, trans_mats = trans_mat)
    cov3 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, indices[3], indices[3], burnin_prop = 0.7, trans_mats = trans_mat)
    cov12 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                    boundary_knots, internal_knots, indices[1], indices[2], burnin_prop = 0.7, trans_mats = trans_mat)
    cov13 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                    boundary_knots, internal_knots, indices[1], indices[3], burnin_prop = 0.7, trans_mats = trans_mat)
    cov23 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                    boundary_knots, internal_knots, indices[2], indices[3], burnin_prop = 0.7, trans_mats = trans_mat)


    # Calculate root mean square error in allocation matrix
    err_Z[i,j] <- sqrt(norm(Z_median - Z_true, "F") / (2 * nrow(Z_median)))

    # Estimate mean integrated square error for mean and covariance functions
    int_err_mean1[i,j] <- sum((nu_1$CI_50 - nu_1_true)^2) * 10
    int_err_mean2[i,j] <- sum((nu_2$CI_50 - nu_2_true)^2) * 10
    int_err_mean3[i,j] <- sum((nu_3$CI_50 - nu_3_true)^2) * 10
    int_err_cov1[i,j] <- sum((cov1$CI_50 - cov1_true)^2) * 100
    int_err_cov2[i,j] <- sum((cov2$CI_50 - cov2_true)^2) * 100
    int_err_cov3[i,j] <- sum((cov3$CI_50 - cov3_true)^2) * 100
    int_err_cov12[i,j] <- sum((cov12$CI_50 - cov12_true)^2) * 100
    int_err_cov13[i,j] <- sum((cov13$CI_50 - cov13_true)^2) * 100
    int_err_cov23[i,j] <- sum((cov23$CI_50 - cov23_true)^2) * 100

    # Estimates of the squared norm of mean and covariance functions
    norm_mu1[i,j] <- sum((nu_1_true)^2) * 10
    norm_mu2[i,j] <- sum((nu_2_true)^2) * 10
    norm_mu3[i,j] <- sum((nu_3_true)^2) * 10
    norm_C1[i,j] <- sum((cov1_true)^2) * 100
    norm_C2[i,j] <- sum((cov2_true)^2) * 100
    norm_C3[i,j] <- sum((cov3_true)^2) * 100
    norm_C12[i,j] <- sum((cov12_true)^2) * 100
    norm_C13[i,j] <- sum((cov13_true)^2) * 100
    norm_C23[i,j] <- sum((cov23_true)^2) * 100

    print(i)
    print(j)

  }
}

int_err_mean1 / norm_mu1
int_err_mean2 / norm_mu2
int_err_mean3 / norm_mu3

# Plot R-MISE for covariance function of feature 1
C1_RMSE <- matrix(0, 30, 2)
C1_RMSE[1:10,1] <- (int_err_cov1[,1] / norm_C1[,1])
C1_RMSE[1:10,2] <- 50
C1_RMSE[11:20,1] <- (int_err_cov1[,2] / norm_C1[,2])
C1_RMSE[11:20,2] <- 100
C1_RMSE[21:30,1] <- (int_err_cov1[,3] / norm_C1[,3])
C1_RMSE[21:30,2] <- 200
C1_RMSE <- as.data.frame(C1_RMSE)
colnames(C1_RMSE) <- c("R-MISE", "N")
C1_RMSE$N <- as.factor(C1_RMSE$N)
p1 <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5,  0.1, 0.02)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))
# Plot R-MISE for covariance function of feature 2
C2_RMSE <- matrix(0, 30, 2)
C2_RMSE[1:10,1] <- (int_err_cov2[,1] / norm_C2[,1])
C2_RMSE[1:10,2] <- 50
C2_RMSE[11:20,1] <- (int_err_cov2[,2] / norm_C2[,2])
C2_RMSE[11:20,2] <- 100
C2_RMSE[21:30,1] <- (int_err_cov2[,3] / norm_C2[,3])
C2_RMSE[21:30,2] <- 200
C2_RMSE <- as.data.frame(C2_RMSE)
colnames(C2_RMSE) <- c("R-MISE", "N")
C2_RMSE$N <- as.factor(C2_RMSE$N)
p2 <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5, 0.1, 0.02)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for covariance function of feature 3
C3_RMSE <- matrix(0, 30, 2)
C3_RMSE[1:10,1] <- (int_err_cov3[,1] / norm_C3[,1])
C3_RMSE[1:10,2] <- 50
C3_RMSE[11:20,1] <- (int_err_cov3[,2] / norm_C3[,2])
C3_RMSE[11:20,2] <- 100
C3_RMSE[21:30,1] <- (int_err_cov3[,3] / norm_C3[,3])
C3_RMSE[21:30,2] <- 200
C3_RMSE <- as.data.frame(C3_RMSE)
colnames(C3_RMSE) <- c("R-MISE", "N")
C3_RMSE$N <- as.factor(C3_RMSE$N)
p3 <- ggplot(C3_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5, 0.1, 0.02)) + ggtitle(TeX("$C^{(3,3)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for cross-covariance function of features 1 and 2
C12_RMSE <- matrix(0, 30, 2)
C12_RMSE[1:10,1] <- (int_err_cov12[,1] / norm_C12[,1])
C12_RMSE[1:10,2] <- 50
C12_RMSE[11:20,1] <- (int_err_cov12[,2] / norm_C12[,2])
C12_RMSE[11:20,2] <- 100
C12_RMSE[21:30,1] <- (int_err_cov12[,3] / norm_C12[,3])
C12_RMSE[21:30,2] <- 200
C12_RMSE <- as.data.frame(C12_RMSE)
colnames(C12_RMSE) <- c("R-MISE", "N")
C12_RMSE$N <- as.factor(C12_RMSE$N)
p4 <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5,  0.1, 0.02)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for cross-covariance function of features 2 and 3
C23_RMSE <- matrix(0, 30, 2)
C23_RMSE[1:10,1] <- (int_err_cov23[,1] / norm_C23[,1])
C23_RMSE[1:10,2] <- 50
C23_RMSE[11:20,1] <- (int_err_cov23[,2] / norm_C23[,2])
C23_RMSE[11:20,2] <- 100
C23_RMSE[21:30,1] <- (int_err_cov23[,3] / norm_C23[,3])
C23_RMSE[21:30,2] <- 200
C23_RMSE <- as.data.frame(C23_RMSE)
colnames(C23_RMSE) <- c("R-MISE", "N")
C23_RMSE$N <- as.factor(C23_RMSE$N)
p5 <- ggplot(C23_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5, 0.1, 0.02)) + ggtitle(TeX("$C^{(2,3)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for cross-covariance function of features 1 and 3
C13_RMSE <- matrix(0, 30, 2)
C13_RMSE[1:10,1] <- (int_err_cov13[,1] / norm_C13[,1])
C13_RMSE[1:10,2] <- 50
C13_RMSE[11:20,1] <- (int_err_cov13[,2] / norm_C13[,2])
C13_RMSE[11:20,2] <- 100
C13_RMSE[21:30,1] <- (int_err_cov13[,3] / norm_C13[,3])
C13_RMSE[21:30,2] <- 200
C13_RMSE <- as.data.frame(C13_RMSE)
colnames(C13_RMSE) <- c("R-MISE", "N")
C13_RMSE$N <- as.factor(C13_RMSE$N)
p6 <- ggplot(C13_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5, 0.1, 0.02)) + ggtitle(TeX("$C^{(1,3)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))



# Plot R-MISE for the mean function of feature 1
mu1_RMSE <- matrix(0, 30, 2)
mu1_RMSE[1:10,1] <- (int_err_mean1[,1] / norm_mu1[,1])
mu1_RMSE[1:10,2] <- 50
mu1_RMSE[11:20,1] <- (int_err_mean1[,2] / norm_mu1[,2])
mu1_RMSE[11:20,2] <- 100
mu1_RMSE[21:30,1] <- (int_err_mean1[,3] / norm_mu1[,3])
mu1_RMSE[21:30,2] <- 200
mu1_RMSE <- as.data.frame(mu1_RMSE)
colnames(mu1_RMSE) <- c("R-MISE", "N")
mu1_RMSE$N <- as.factor(mu1_RMSE$N)
p7 <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.1, 0.02, 0.005, 0.001))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for the mean function of feature 2
mu2_RMSE <- matrix(0, 30, 2)
mu2_RMSE[1:10,1] <- (int_err_mean2[,1] / norm_mu2[,1])
mu2_RMSE[1:10,2] <- 50
mu2_RMSE[11:20,1] <- (int_err_mean2[,2] / norm_mu2[,2])
mu2_RMSE[11:20,2] <- 100
mu2_RMSE[21:30,1] <- (int_err_mean2[,3] / norm_mu2[,3])
mu2_RMSE[21:30,2] <- 200
mu2_RMSE <- as.data.frame(mu2_RMSE)
colnames(mu2_RMSE) <- c("R-MISE", "N")
mu2_RMSE$N <- as.factor(mu2_RMSE$N)
p8 <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.1, 0.02, 0.005, 0.001)) + ggtitle(TeX("$\\mu_2$"))+
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for the mean function of feature 1
mu3_RMSE <- matrix(0, 30, 2)
mu3_RMSE[1:10,1] <- (int_err_mean3[,1] / norm_mu3[,1])
mu3_RMSE[1:10,2] <- 50
mu3_RMSE[11:20,1] <- (int_err_mean3[,2] / norm_mu3[,2])
mu3_RMSE[11:20,2] <- 100
mu3_RMSE[21:30,1] <- (int_err_mean3[,3] / norm_mu3[,3])
mu3_RMSE[21:30,2] <- 200
mu3_RMSE <- as.data.frame(mu3_RMSE)
colnames(mu3_RMSE) <- c("R-MISE", "N")
mu3_RMSE$N <- as.factor(mu3_RMSE$N)
p9 <- ggplot(mu3_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.1, 0.02, 0.005, 0.001))  + ggtitle(TeX("$\\mu_3$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot the RMSE for the allocation parameters
Z_RMSE <- matrix(0, 30, 2)
Z_RMSE[1:10,1] <- err_Z[,1]
Z_RMSE[1:10,2] <- 50
Z_RMSE[11:20,1] <- err_Z[,2]
Z_RMSE[11:20,2] <- 100
Z_RMSE[21:30,1] <- err_Z[,3]
Z_RMSE[21:30,2] <- 200
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N")
Z_RMSE$N <- as.factor(Z_RMSE$N)
p10 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") + scale_y_continuous( trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.08, 0.04, 0.02, 0.01)) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                       legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))


# Arrange all plots together
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2)
