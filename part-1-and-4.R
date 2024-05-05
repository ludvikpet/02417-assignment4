rm(list=ls())
library(ggplot)
library("FKF")
cwd <- setwd("C:/Users/margr/OneDrive - Danmarks Tekniske Universitet/Skrivebord/DTU/sem10/time-series-analysis/assignments/assignment4/02417-assignment4/")
data_path <- paste(cwd, "/data/", sep="")


# Load data:
data <- list()
for (i in 1:4) {
  name <- paste("rain", i, ".csv", sep="")
  data[[i]] <- read.csv(paste(data_path, name, sep=""))
  
  # Add column for cummulative sum of u elements:
  data[[i]]$u_cumsum <- cumsum(data[[i]]$u)
}
data[1]
max(unlist(data[[1]]$u))
max(data[[i]]$minutes)

n <- length(data[[1]]$minutes)

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
plot_data(data)
print(data[[1]]$y)
################
#### PART 4 ####
################

# Read in data:
A <- function(a) {
  return(matrix(c(1 - a, 0, 0, 0,
                  a, 1 - a, 0, 0,
                  0, a, 1 - a, 0,
                  0, 0, a, 0.98), nrow = 4, ncol = 4, byrow=TRUE))
}

B <- matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1)

C <- matrix(c(0, 0, 0, 1), nrow = 4, ncol = 1)

G <- function(X_t1) {
  return(matrix(c(sqrt(X_t1[1]), 0, 0, 0, # 1st, 2nd, 3rd, 4th water bassin
                  0, sqrt(X_t1[2]), 0, 0,
                  0, 0, sqrt(X_t1[3]), 0,
                  0, 0, 0, sqrt(X_t1[4])), nrow = 4, ncol = 4, byrow=TRUE))
}

V1 <- function(sigma_sq_1) { # System noise
  return(matrix(c(sigma_sq_1, 0, 0, 0, # 1st, 2nd, 3rd, 4th water bassin
                  0, sigma_sq_1, 0, 0,
                  0, 0, sigma_sq_1, 0,
                  0, 0, 0, sigma_sq_1), nrow = 4, ncol = 4, byrow=TRUE))
}

u <- function(data, t) {
  return(matrix(c(data[4]$y[t], 0, 0, 0, # 1st, 2nd, 3rd, 4th water bassin
                  0, data[4]$y[t], 0, 0,
                  0, 0, data[4]$y[t], 0,
                  0, 0, 0, data[4]$y[t]), nrow = 4, ncol = 4, byrow=TRUE))
}

# Initial values:
X0 <- c(0, 0, 0, 0)
a <- 0.1 # transition rate, range: [0.01, 0.1]
Sigma1 <- 1  # variance range: [0.001; 1], hence we start with largest
Sigma2 <- 5 # variance range: [0.1; 5], hence we start with largest

##################
#### PART 4.1 ####
##################

y <- rep(NA, n) # Store output
lik <- rep(NA, n)
X <- X0
SigmaX0 <- V1(Sigma1)

for (i in 2:n){
    X <- A %*% X + B%*%g + chol(Sigma1) %*% matrix(rnorm(2),ncol=1)
    y[i] <- C %*% X + sqrt(Sigma2) %*% rnorm(1)
}

# Plot the simulated height
plot(y, type="l", xlab="time", ylab="Observed altitude [m]", ylim=range(y,na.rm=TRUE))

# Calculate the negative log likelihood using the Kalman filter:
negloglik <- function(prm){
  ## --------------------------------
  # Overwrite for the parameters given
  all[names(prm)] <- prm
  # Recreate the parameter values
  A <- matrix(all[c("A1","A2","A3","A4")], nrow=2)
  B <- matrix(all[c("B1","B2")], nrow=2)
  C <- matrix(all[c("C1","C2")], nrow=1)
  Sigma1 <- matrix(all[c("Sigma11","Sigma12","Sigma13","Sigma14")], nrow=2)
  Sigma2 <- matrix(all["Sigma2"])
  # Initial state vector and covariance
  X <- all[c("X01","X02")]
  SigmaX <- matrix(all[c("SigmaX01","SigmaX02","SigmaX03","SigmaX04")], nrow=2)
  ## --------------------------------
  # For keeping results
  n <- length(y)
  lik <- rep(NA,n)
  # Start on the iterations
  i <- 1
  while(i < n){
    # ---- PREDICT ----
    X <- A %*% X + B %*% g
    SigmaX <- A %*% SigmaX %*% t(A) + Sigma1
    SigmaY <- C %*% SigmaX %*% t(C) + Sigma2
    # ---- Step ----
    i <- i + 1
    # The likelihood
    lik[i] <- dnorm(y[i], mean=C%*%X, sd=sqrt(SigmaY))
    # ---- UPDATE ----
    K <- SigmaX %*% t(C) %*% solve(SigmaY)
    X <- X + K %*% (y[i] - C %*% X)
    SigmaX <- SigmaX - K %*% SigmaY %*% t(K)
  }
  # Calculate the negative log-likelihood
  return(-sum(log(lik[!is.na(lik)])))
}

# Prm in vector for estimation:
prm <- c(
  A = c(A(a)),
  B = c(B),
  C = c(C),
  G = c(G(X0)),
  Sigma1 = Sigma1,
  Sigma2 = Sigma2,
  X0 = X0,
  SigmaX0 = c(SigmaX0)
)

# Estimate the parameters using MLE: