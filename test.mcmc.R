library(devtools)
setwd("~/Desktop/Chakravarti_Lab/git")
load_all("CNPBayes_trios")
#setwd("~/Desktop/Chakravarti_Lab/git/triostat")
#load_all("triostat")
library(gtools)

p <- c(0.25, 0.5, 0.25)
theta <- c(-3,0.15, 1.2)
sigma2 <- c(0.1, 0.2, 0.2)
N <- 30
params <- data.frame(cbind(p, theta, sigma2))

truth.plot <- data.frame(truth$data$log_ratio, as.factor(truth$data$copy_number))
names(truth.plot) <- c("log_ratio", "copy_number")
ggplot(truth.plot, aes(log_ratio, fill = copy_number)) +
geom_histogram(bins=100)


seed <- 4990917
set.seed(seed)

# note maplabel and hp defined manually here for now
maplabel <- c(0,1,2)
k <- length(maplabel)
hp <- HyperparametersTrios(k = 3)
nbatch <- 1
mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
truth <- simulate_data_multi2(params, N=N,
                              batches = rep(c(1:nbatch),
                                            length.out = 3*N),
                              error=0, mprob, maplabel)

#mp <- McmcParams(iter=4000, burnin=4000, thin=1, nStarts=3)

#tbm.init <- TBM(triodata=truth$data,
#            mprob=mprob,
#            maplabel=maplabel)

true.stats <- component_stats(truth$data)
thetas <- true.stats$mean
sigmas <- true.stats$sd
pis <- true.stats$p

# check the two CN combinations we are interested in exists
# CN2,1,0 and CN 2,1,1
cn.df <- data.frame(matrix(ncol=3, data=truth$data$copy_number))
# row 5 and 6 has 2,1,1 and 1,2,1
# row 18 and 29 has 1,2,0 and 2,1,0 combinations
colnames(cn.df) <- c("cn.mot","cn.fat","cn.off")

# generate the large matrix
genotype.df <- data.frame(permutations(3,3,c(0:2), repeats.allowed=T))
colnames(genotype.df) <- c("cn.mot","cn.fat","cn.off")

# simplifying assumption - thetas, sigmas,pis fixed
genotype.df$mott <- thetas[(genotype.df$cn.mot+1)]
genotype.df$fatt <- thetas[(genotype.df$cn.fat+1)]
genotype.df$offt <- thetas[(genotype.df$cn.fat+1)]
genotype.df$mots <- sigmas[(genotype.df$cn.mot+1)]
genotype.df$fats <- sigmas[(genotype.df$cn.fat+1)]
genotype.df$offs <- sigmas[(genotype.df$cn.fat+1)]
genotype.df$motp <- pis[(genotype.df$cn.mot+1)]
genotype.df$fatp <- pis[(genotype.df$cn.fat+1)]
genotype.df$offp <- pis[(genotype.df$cn.fat+1)]
mendel <- rbinom(n=nrow(cn.df), 1, prob=0.9)
mendel.count <- sum(mendel)

# sort out the data
inten.mat <- data.frame(matrix(ncol=3, data=truth$data$log_ratio))
colnames(inten.mat) <- c("cmot", "cfat", "coff")

# compare row 5 vs 29
inten.mat1 <- inten.mat[5,]
inten.mat2 <- inten.mat[29,]

# test mcmc loop
genotype.df$cmot.dens <- dnorm(inten.mat2$cmot, mean=genotype.df$mott, sd=(genotype.df$mots)^0.5)
genotype.df$cfat.dens <- dnorm(inten.mat2$cfat, mean=genotype.df$fatt, sd=(genotype.df$fats)^0.5)
genotype.df$coff.dens <- dnorm(inten.mat2$coff, mean=genotype.df$offt, sd=(genotype.df$offs)^0.5)
eta <- rbeta(1, 1+mendel.count, 1+(N-mendel.count))
mendel <- rbinom(n=nrow(cn.df), 1, prob=eta)
mendel.count <- sum(mendel)
genotype.df$m0 <- genotype.df$motp * genotype.df$fatp * genotype.df$offp * (1-eta)
mprob.unlist <- mprob[,c(1:3)]
trio.prob <- c(t(mprob.unlist))
genotype.df$tp <- trio.prob
genotype.df$m1 <- genotype.df$motp * genotype.df$fatp * genotype.df$tp * eta
genotype.df$nump <- genotype.df$cmot.dens * genotype.df$cfat.dens * genotype.df$coff.dens * (genotype.df$m0 + genotype.df$m1)
geno.denom <- sum(genotype.df$nump)
genotype.df$prob <- genotype.df$nump / geno.denom
row.select <- sample(1:27, 1, prob=genotype.df$prob)
genotype.assign <- genotype.df[row.select,c(1:3)]

# the actual MCMC - here we test inten.mat2 (CN2,1,0)
genotype.mat <- matrix(nrow=20, ncol=3)
genotype.211 <- vector(length=20)
genotype.210 <- vector(length=20)

for (i in 1:20){
  
  genotype.df$cmot.dens <- dnorm(inten.mat2$cmot, mean=genotype.df$mott, sd=(genotype.df$mots)^0.5)
  genotype.df$cfat.dens <- dnorm(inten.mat2$cfat, mean=genotype.df$fatt, sd=(genotype.df$fats)^0.5)
  genotype.df$coff.dens <- dnorm(inten.mat2$coff, mean=genotype.df$offt, sd=(genotype.df$offs)^0.5)
  eta <- rbeta(1, 1+mendel.count, 1+(N-mendel.count))
  mendel <- rbinom(n=nrow(cn.df), 1, prob=eta)
  mendel.count <- sum(mendel)
  genotype.df$m0 <- genotype.df$motp * genotype.df$fatp * genotype.df$offp * (1-eta)
  mprob.unlist <- mprob[,c(1:3)]
  trio.prob <- c(t(mprob.unlist))
  genotype.df$tp <- trio.prob
  genotype.df$m1 <- genotype.df$motp * genotype.df$fatp * genotype.df$tp * eta
  genotype.df$nump <- genotype.df$cmot.dens * genotype.df$cfat.dens * genotype.df$coff.dens * (genotype.df$m0 + genotype.df$m1)
  geno.denom <- sum(genotype.df$nump)
  genotype.df$prob <- genotype.df$nump / geno.denom
  row.select <- sample(1:27, 1, prob=genotype.df$prob)
  genotype.assign <- as.numeric(genotype.df[row.select,c(1:3)])
  
  genotype.mat[i,] <- genotype.assign
  genotype.210[i] <- genotype.df$prob[22]
  genotype.211[i] <- genotype.df$prob[23]
}

mean(genotype.210)
mean(genotype.211)
