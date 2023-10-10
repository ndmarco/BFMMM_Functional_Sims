library(BayesFMMM)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(pracma)

### Set working dir
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
err_Z <- matrix(0, 50, 3)
int_err_mean1 <- matrix(0, 50, 3)
int_err_mean2 <- matrix(0, 50, 3)
int_err_cov12 <- matrix(0, 50, 3)
int_err_cov1 <- matrix(0, 50, 3)
int_err_cov2 <- matrix(0, 50, 3)

# Initialize placeholders for estimates of the squared norm of mean and covariance functions
norm_mu1 <- matrix(1, nrow = 50, ncol = 3)
norm_mu2 <- matrix(1, nrow = 50, ncol = 3)
norm_C1 <- matrix(1, nrow = 50, ncol = 3)
norm_C2 <- matrix(1, nrow = 50, ncol = 3)
norm_C12 <- matrix(1, nrow = 50, ncol = 3)

# Estimate the R-MISE for each of the 3 sample sizes (N = 40, N = 80, N = 160)
for(j in 1:3){
  if(j == 1){
    dir <- "./40_obs/"
  }
  if(j == 2){
    dir <- "./80_obs/"
  }
  if(j == 3){
    dir <- "./160_obs/"
  }
  for(i in 1:50){
    # Read in true parameters
    x <- readRDS(paste(dir, "sim", i, "/truth.RDS", sep=""))
    x$Phi_true <- x$Phi
    x$nu_true <- x$nu
    Z_true <- x$Z

    # Get true means
    nu_1_true <-  B[[1]] %*% t(t(x$nu_true[1,]))
    nu_2_true <-  B[[1]] %*% t(t(x$nu_true[2,]))

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

    cov12_true <- matrix(0, 100, 100)
    for(q in 1:100){
      for(k in 1:100){
        for(m in 1:3){
          cov12_true[q,k] <- cov12_true[q,k] + x$Phi_true[1, ,m] %*% t(t(B[[1]][q,])) %*% B[[1]][k,] %*% t(t(x$Phi_true[2, ,m]))
        }
      }
    }

    # Hyperparameters for the sampling already conducted
    n_files <- 50
    basis_degree <- 3
    boundary_knots <- c(0, 1000)
    internal_knots <- c(200, 400, 600, 800)
    dir_i <- paste(dir, "sim", i, "/", sep = "")

    # Get estimated mean, covariance, and cross-covariance functions
    nu_1 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 1, burnin_prop = 0.7)
    nu_2 <- FMeanCI(dir_i, n_files, time[[1]], basis_degree, boundary_knots, internal_knots, 2, burnin_prop = 0.7)
    cov1 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 1, 1, burnin_prop = 0.7)
    cov2 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                   boundary_knots, internal_knots, 2, 2, burnin_prop = 0.7)
    cov12 <- FCovCI(dir_i, n_files, time[[1]], time[[1]], basis_degree,
                    boundary_knots, internal_knots, 1, 2, burnin_prop = 0.7)

    # Get estimated allocation parameters
    Z_est <- ZCI(dir_i, n_files)

    ## Deal with potential label switching issues
    if(sum(abs(nu_1$CI_50 - nu_1_true)) > sum(abs(nu_1$CI_50 - nu_2_true))){
      Z_est_i <- Z_est
      Z_est$CI_50[,1] <- Z_est_i$CI_50[,2]
      Z_est$CI_50[,2] <- Z_est_i$CI_50[,1]
      nu_i <- nu_1
      nu_1 <- nu_2
      nu_2 <- nu_i
      cov_i <- cov1
      cov1 <- cov2
      cov2 <- cov_i
      cov12$CI_50<- t(cov12$CI_50)
      cov12$CI_Upper<- t(cov12$CI_Upper)
      cov12$CI_Lower<- t(cov12$CI_Lower)
    }
    # Calculate root mean square error in allocation matrix
    err_Z[i,j] <- sqrt(norm(Z_est$CI_50 - Z_true, "F") / (2 * nrow(Z_est$CI_50)))

    # Estimate mean integrated square error for mean and covariance functions
    int_err_mean1[i,j] <- sum((nu_1$CI_50 - nu_1_true)^2) * 10
    int_err_mean2[i,j] <- sum((nu_2$CI_50 - nu_2_true)^2) * 10
    int_err_cov1[i,j] <- sum((cov1$CI_50 - cov1_true)^2) * 100
    int_err_cov2[i,j] <- sum((cov2$CI_50 - cov2_true)^2) * 100
    int_err_cov12[i,j] <- sum((cov12$CI_50 - cov12_true)^2) * 100

    # Estimates of the squared norm of mean and covariance functions
    norm_mu1[i,j] <- sum((nu_1_true)^2) * 10
    norm_mu2[i,j] <- sum((nu_2_true)^2) * 10
    norm_C1[i,j] <- sum((cov1_true)^2) * 100
    norm_C2[i,j] <- sum((cov2_true)^2) * 100
    norm_C12[i,j] <- sum((cov12_true)^2) * 100

    print(i)
    print(j)

  }
}



# Plot R-MISE for covariance function of feature 1
C1_RMSE <- matrix(0, 150, 2)
C1_RMSE[1:50,1] <- (int_err_cov1[,1] / norm_C1[,1])
C1_RMSE[1:50,2] <- 40
C1_RMSE[51:100,1] <- (int_err_cov1[,2] / norm_C1[,2])
C1_RMSE[51:100,2] <- 80
C1_RMSE[101:150,1] <- (int_err_cov1[,3] / norm_C1[,3])
C1_RMSE[101:150,2] <- 160
C1_RMSE <- as.data.frame(C1_RMSE)
colnames(C1_RMSE) <- c("R-MISE", "N")
C1_RMSE$N <- as.factor(C1_RMSE$N)
p1 <- ggplot(C1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5,  0.1, 0.02)) + ggtitle(TeX("$C^{(1,1)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))
# Plot R-MISE for covariance function of feature 2
C2_RMSE <- matrix(0, 150, 2)
C2_RMSE[1:50,1] <- (int_err_cov2[,1] / norm_C2[,1])
C2_RMSE[1:50,2] <- 40
C2_RMSE[51:100,1] <- (int_err_cov2[,2] / norm_C2[,2])
C2_RMSE[51:100,2] <- 80
C2_RMSE[101:150,1] <- (int_err_cov2[,3] / norm_C2[,3])
C2_RMSE[101:150,2] <- 160
C2_RMSE <- as.data.frame(C2_RMSE)
colnames(C2_RMSE) <- c("R-MISE", "N")
C2_RMSE$N <- as.factor(C2_RMSE$N)
p2 <- ggplot(C2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5, 0.1, 0.02)) + ggtitle(TeX("$C^{(2,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(),  axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for cross-covariance function
C12_RMSE <- matrix(0, 150, 2)
C12_RMSE[1:50,1] <- (int_err_cov12[,1] / norm_C12[,1])
C12_RMSE[1:50,2] <- 40
C12_RMSE[51:100,1] <- (int_err_cov12[,2] / norm_C12[,2])
C12_RMSE[51:100,2] <- 80
C12_RMSE[101:150,1] <- (int_err_cov12[,3] / norm_C12[,3])
C12_RMSE[101:150,2] <- 160
C12_RMSE <- as.data.frame(C12_RMSE)
colnames(C12_RMSE) <- c("R-MISE", "N")
C12_RMSE$N <- as.factor(C12_RMSE$N)
p3 <- ggplot(C12_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 100), breaks = c(10, 0.5, 0.1, 0.02)) + ggtitle(TeX("$C^{(1,2)}$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for the mean function of feature 1
mu1_RMSE <- matrix(0, 150, 2)
mu1_RMSE[1:50,1] <- (int_err_mean1[,1] / norm_mu1[,1])
mu1_RMSE[1:50,2] <- 40
mu1_RMSE[51:100,1] <- (int_err_mean1[,2] / norm_mu1[,2])
mu1_RMSE[51:100,2] <- 80
mu1_RMSE[101:150,1] <- (int_err_mean1[,3] / norm_mu1[,3])
mu1_RMSE[101:150,2] <- 160
mu1_RMSE <- as.data.frame(mu1_RMSE)
colnames(mu1_RMSE) <- c("R-MISE", "N")
mu1_RMSE$N <- as.factor(mu1_RMSE$N)
p4 <- ggplot(mu1_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.1, 0.02, 0.005, 0.001))  + ggtitle(TeX("$\\mu_1$")) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot R-MISE for the mean function of feature 2
mu2_RMSE <- matrix(0, 150, 2)
mu2_RMSE[1:50,1] <- (int_err_mean2[,1] / norm_mu2[,1])
mu2_RMSE[1:50,2] <- 40
mu2_RMSE[51:100,1] <- (int_err_mean2[,2] / norm_mu2[,2])
mu2_RMSE[51:100,2] <- 80
mu2_RMSE[101:150,1] <- (int_err_mean2[,3] / norm_mu2[,3])
mu2_RMSE[101:150,2] <- 160
mu2_RMSE <- as.data.frame(mu2_RMSE)
colnames(mu2_RMSE) <- c("R-MISE", "N")
mu2_RMSE$N <- as.factor(mu2_RMSE$N)
p5 <- ggplot(mu2_RMSE, aes(x=N, y=`R-MISE`)) + scale_y_continuous(labels = scales::percent, trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.1, 0.02, 0.005, 0.001)) + ggtitle(TeX("$\\mu_2$"))+
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Plot the RMSE for the allocation paramaeters
Z_RMSE <- matrix(0, 150, 2)
Z_RMSE[1:50,1] <- err_Z[,1]
Z_RMSE[1:50,2] <- 40
Z_RMSE[51:100,1] <- err_Z[,2]
Z_RMSE[51:100,2] <- 80
Z_RMSE[101:150,1] <- err_Z[,3]
Z_RMSE[101:150,2] <- 160
Z_RMSE <- as.data.frame(Z_RMSE)
colnames(Z_RMSE) <- c("RMSE", "N")
Z_RMSE$N <- as.factor(Z_RMSE$N)
p6 <- ggplot(Z_RMSE, aes(x=N, y=`RMSE`)) + ggtitle("Z") + scale_y_continuous( trans = 'log2', minor_breaks = scales::pretty_breaks(n = 10), breaks = c(0.08, 0.04, 0.02, 0.01)) +
  geom_boxplot() +  theme_bw() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"),
                                            legend.position = "none",plot.title = element_text(hjust = 0.5), text = element_text(size = 15))

# Arrange all plots together
grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                                                       c(4, 4, 5, 5, 6, 6)))
