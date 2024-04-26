rm(list=ls())
library(ggplot)
cwd <- setwd("C:/Users/margr/OneDrive - Danmarks Tekniske Universitet/Skrivebord/DTU/sem10/time-series-analysis/assignments/assignment4/")
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
max(unlist(data[1]$u))
max(data[[i]]$minutes)

# Plot data:
plot_data <- function(data) {
  file_path <- paste(cwd, "/test.png", sep="")
  png(filename = file_path, width = 1200, height = 700)
  par(mfrow = c(2, 2))
  # for (i in 1:length(data)) {
  #   x <- data[[i]]$minutes
  #   u <- data[[i]]$u*10
  #   u_cumsum <- data[[i]]$u_cumsum
  #   y <- data[[i]]$y
  #   title <- paste("rain", i, " - Water level (y) and incomming rain (u) over time", sep="")
  #   plot(x, y, type = "l", col = "blue", ylim= c(0,max(data[[1]]$y)), xlim = c(0, max(data[[4]]$minutes)), xlab = "Time (minutes)", ylab = "Water (100 m^3)", main = title)
  #   lines(x, u, col = "red")
  #   # lines(x, u_cumsum, col = "green")
  #   legend("topright", legend = c("y (water in bassin)", "u (output pr. minute, scaled x 10)"), col = c("blue", "red"), lty = 1)
  # }
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
    axis(side = 4, col = "red", col.axis = "red", las = 1)
    mtext("Output pr. minute (u, 100 m^3)", side = 4, line = 3, col = "red")
    legend("topright", legend = c("y (water in bassin)", "u (output pr. minute)"), col = c("blue", "red"), lty = 1)
  }
  dev.off()
}
plot_data(data)
