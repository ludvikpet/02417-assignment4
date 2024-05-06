# Preliminaries:
rm(list=ls())
library(ggplot2)
library(FKF)
library(MASS)
library(hash)
library(numDeriv)
cwd <- setwd("C:/Users/margr/OneDrive - Danmarks Tekniske Universitet/Skrivebord/DTU/sem10/time-series-analysis/assignments/assignment4/02417-assignment4/")
data_path <- paste(cwd, "/data/", sep="")

################
#### PART 1 ####
################

# Load data:
retrieve_data <- function(data_path){
  data <- list()
  for (i in 1:4) {
    name <- paste("rain", i, ".csv", sep="")
    print(name)
    data[[i]] <- read.csv(paste(data_path, name, sep=""))
    
    # Add column for cummulative sum of u elements:
    data[[i]]$u_cumsum <- cumsum(data[[i]]$u)
  }
  return(data)
}
data <- retrieve_data(data_path)

# Test if data has been loaded correctly:
data[[1]]

# Plot data:
plot_data <- function(data) {
  file_path <- paste(cwd, "/water_over_time.png", sep="")
  png(filename = file_path, width = 1200, height = 700)
  par(mfrow = c(2, 2))
  for (i in 1:length(data)) {
    x <- data[[i]]$minutes
    u <- data[[i]]$u
    u_cumsum <- data[[i]]$u_cumsum
    y <- data[[i]]$y
    title <- paste("rain", i, " - Water level (y) and incomming rain (u) over time", sep="")
    par(mar = c(5, 4, 4, 8))
    plot(x, y, type = "l", col = "blue", ylim= c(0,max(data[[1]]$y)), xlim = c(0, max(data[[4]]$minutes)), xlab = "Time (minutes)", ylab = "Water (100 m^3)", main = title)
    par(new = TRUE)
    plot(x, u, type = "l", col = "red", ylim = c(0, max(data[[3]]$u)), ylab = "", xlab = "", axes = FALSE)
    axis(side = 4, las = 1)
    mtext("Output pr. minute (u, 100 m^3)", side = 4, line = 3)
    legend("topright", legend = c("y (water in bassin)", "u (output pr. minute)"), col = c("blue", "red"), lty = 1)
  }
  dev.off()
}
# plot_data(data) # Uncomment to create .png of the parameters over time for the four bassins 

################
#### PART 4 ####
################

# Read in data:
A <- function(a) {
  return(matrix(c(1 - a, 0, 0, 0,
                  a, 1 - a, 0, 0,
                  0, a, 1 - a, 0,
                  0, 0, a, 0.98), nrow = 4, byrow=TRUE))
}

B <- matrix(c(1, 0, 0, 0), ncol = 1)

C <- matrix(c(0, 0, 0, 1), ncol = 4)

G <- function(X_t1) {
  return(matrix(c(abs(X_t1[1]), 0, 0, 0, # 1st, 2nd, 3rd, 4th water bassin
                  0, abs(X_t1[2]), 0, 0,
                  0, 0, abs(X_t1[3]), 0,
                  0, 0, 0, abs(X_t1[4])), nrow = 4, byrow=TRUE)) # Only take the abs val, not sqrt
}

V1 <- function(sigma_sq_1) { # System noise
  return(diag(rep(sigma_sq_1), 4))
}

u <- function(idx, t) {
  return (c(data[[idx]]$u[t]))
}

##################
#### PART 4.1 ####
##################

# Calculate the negative log likelihood using the Kalman filter:
negloglik <- function(prm, plot=FALSE) {
  ## --------------------------------
  # Overwrite for the parameters given
  all[names(prm)] <- prm
  # Recreate the parameter values
  a_prm <- prm["a"]
  # A <- matrix(all[c("A1","A2","A3","A4",
  #                   "A5","A6","A7","A8",
  #                   "A9","A10","A11","A12",
  #                   "A13","A14","A15","A16")], nrow=4, byrow=TRUE)
  A_new <- A(a_prm)
  B <- matrix(all[c("B1","B2","B1","B2")], nrow=4)
  C <- matrix(all[c("C1","C2","C3","C4")], nrow=1)
  S2 <- matrix(all["Sigma2"])
  # Initial state vector and covariance
  X <- all[c("X01","X02", "X03", "X04")]
  SigmaX <- matrix(all[c("SigmaX01","SigmaX02","SigmaX03","SigmaX04",
                         "SigmaX05","SigmaX06","SigmaX07","SigmaX08",
                         "SigmaX09","SigmaX010","SigmaX011","SigmaX012",
                         "SigmaX013","SigmaX014","SigmaX015","SigmaX016")], nrow=4, byrow=TRUE)
  S1 <- prm["Sigma1"]
  
  ## --------------------------------
  # For keeping results
  n <- length(y)
  lik <- rep(NA,n)
  idx <- h[[paste0(n)]]
  df <- data.frame(
    time = data[[4]]$minutes[1:n],
    y_pred = rep(NA, n),
    y = y,
    XX_pred_noise = rep(NA, n),
    YY_pred_noise = rep(NA, n)
  )
  
  # Start on the iterations
  i <- 1
  while(i < n) {
    # ---- PREDICT ----
    SigmaX <- A_new %*% SigmaX %*% t(A_new) + G(X) %*% V1(S1) %*% t(G(X))
    SigmaY <- C %*% SigmaX %*% t(C) + S2 # Based on S1_{t+1}
    X <- A_new %*% X + B %*% u(idx, i) # Update after variance calculations!
    
    # Keep the prediction:
    df$y_pred[i] <- C %*% X
    df$XX_pred_noise[i] <- sqrt(SigmaX[4,4])
    df$YY_pred_noise[i] <- sqrt(SigmaY)
    
    # ---- Step ----
    i <- i + 1
    # The likelihood
    lik[i] <- dnorm(y[[i]], mean=C%*%X, sd=sqrt(SigmaY))
    
    K <- SigmaX %*% t(C) %*% solve(SigmaY)
    X <- X + K %*% (y[i] - C %*% X)
    SigmaX <- SigmaX - K %*% SigmaY %*% t(K)
  }
  
  # Plot the predicted and the reconstructed output with variance:
  if(plot == TRUE) {
    par(mar = c(5, 4, 4, 8))
    title <- paste("Predicted and ground truth output with variance bounds (rain", idx, ")", sep="")
    plot(df$time, df$y, type="l", col="blue", ylim=range(df$y,na.rm=TRUE), xlab="Time (minutes)", ylab="Water (100 m^3)", main=title)
    lines(df$time, df$y_pred, col="red")
    lines(df$time, df$y_pred + df$YY_pred_noise, col="green", lty=2)
    lines(df$time, df$y_pred - df$YY_pred_noise, col="green", lty=2)
    legend("topright", legend = c("Observed", "Predicted", "Pred. obs. noise"), col = c("blue", "red", "green"), lty = 1)
  } # Had to add it here, which is stupid (R...)
  
  # Calculate the negative log-likelihood
  return(-sum(log(lik[!is.na(lik)])))
}

# Define the initial parameter values and bounds used for estimation:
a_init <- 0.1  # Initial guess for the transition rate 'a'
S1_init <- 1000 # As long as uncertainty is large enough, choice of X is irrelevant!
S2_init <- 5
initial_params <- c(a = a_init, Sigma1 = S1_init, Sigma2 = S2_init)  # Including other parameters if they are to be optimized

# Define bounds
lower_bounds <- c(a = 0.01, Sigma1 = 0.001, Sigma2 = 0.1)  # Lower bounds for each parameter
upper_bounds <- c(a = 0.1, Sigma1 = 1, Sigma2 = 5)  # Upper bounds for each parameter

X0 <- c(0,0,0,0)
SigmaX0 <- V1(S1_init)

# Prm in vector for estimation:
all <- c(
  A = c(A(a_init)),
  B = c(B),
  C = c(C),
  Sigma2 = S2_init,
  X0 = X0,
  SigmaX0 = c(SigmaX0)
)
h <- hash()
h[["378"]] <- 1
h[["420"]] <- 2
h[["366"]] <- 3
h[["481"]] <- 4 # Find a better solution than this ...

##################
#### PART 4.2 ####
##################

negloglik(initial_params, plot=TRUE)
# Estimate the parameters using MLE:
res <- optim(initial_params, negloglik, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, control = list(maxit = 1000))
negloglik(res$par, plot=TRUE)
# What we see from the plot: Given a wrong starting guess, as the variance is large enough, we'll approximate the data
# after some small time-horizon. I.e. what do we learn? If we add an adequate uncertainty, then the filter will be able
# to model well to the data.

# However, why does the predicted values at entry increase when y is constant at 0?

#----------------------------------------------------------------
# What about the estimated uncertainty of the parameters with MLE?
#----------------------------------------------------------------
# Approximate the Hessian
H <- hessian(negloglik, res$par)
H
res
# Fisher's information matrix
I <- qr.solve(H)

# Covariance matrix of the estimates
cov2cor(I)

# Standard error (Standard deviation estimate of the theta parameters, so mu and exp(sigma^2))
sqrt(diag(I))

##################
#### PART 4.3 ####
##################
data[[4]]$y
res_lst <- list()
for(i in 1:length(data)) {
  y <- data[[i]]$y
  res <- optim(initial_params, negloglik, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, control = list(maxit = 1000))
  
  res_lst[[i]] = res
}

# Result presentation:
res_lst[[1]]$par
res_lst[[2]]$par
res_lst[[3]]$par
res_lst[[4]]$par

##################
#### PART 4.3 ####
##################

file_path <- paste(cwd, "/model_validation.png", sep="")
png(filename = file_path, width = 1200, height = 700)
par(mfrow = c(2, 2))
for(i in 1:length(data)) { 
  y <- data[[i]]$y
  negloglik(res_lst[[i]]$par, plot=TRUE) # Creates plot to be added to plot grid
}
dev.off()
