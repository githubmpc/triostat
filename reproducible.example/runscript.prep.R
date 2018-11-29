# directory where run files are located/ to be saved
#setwd("/gpfs/data/chaklab/home/choum02/batchrun27")
setwd("~/Desktop/Chakravarti_Lab/git/testing scripts/reproducible.example")

simsList <- vector("list", 12)

# Sim1
p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 50
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[1]] <- params

# Sim2
p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 100
params <- data.frame(cbind(p, theta, sigma2, N))

# plots
#truth.plot <- data.frame(truth$data$log_ratio, as.factor(truth$data$copy_number))
#names(truth.plot) <- c("log_ratio", "copy_number")
#ggplot(truth.plot, aes(log_ratio, fill = copy_number)) +
#  geom_histogram(bins=100)
simsList[[2]] <- params

# Sim3
p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 250
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[3]] <- params

# Sim4
p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[4]] <- params

# Sim 5
p <- c(0.09, 0.42, 0.49)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 50
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[5]] <- params

# Sim 6
p <- c(0.09, 0.42, 0.49)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 100
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[6]] <- params

# Sim 7
p <- c(0.09, 0.42, 0.49)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 250
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[7]] <- params

# Sim 8
p <- c(0.09, 0.42, 0.49)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[8]] <- params

# Sim 9
p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 50
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[9]] <- params

# Sim 10
p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 100
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[10]] <- params

# Sim 11
p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 250
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[11]] <- params

# Sim 12
p <- c(0.01, 0.18, 0.81)
theta <- c(-4,-2, 0)
sigma2 <- c(1, 1, 1)
N <- 1000
params <- data.frame(cbind(p, theta, sigma2, N))
simsList[[12]] <- params

params.input <- unlist(simsList)
sim.params <- data.frame(matrix(data=params.input, nrow=12, ncol=12, byrow=T))
sim.params <- sim.params[,-c(11:12)]
colnames(sim.params) <- c("p1", "p2", "p3", "theta1", "theta2", "theta3", "sigma2a", "sigma2b", "sigma2c", "N")
saveRDS(sim.params,"sim.params.overlap.rds")