library(devtools)
setwd("~/Desktop/Chakravarti_Lab/git")
load_all("CNPBayes_trios")
#setwd("~/Desktop/Chakravarti_Lab/git/triostat")
#load_all("triostat")
library(gtools)
library(tidyverse)
library(ggplot2)

p <- c(0.25, 0.5, 0.25)
theta <- c(-3,0.15, 1.2)
sigma2 <- c(0.1, 0.2, 0.2)
N <- 30
params <- data.frame(cbind(p, theta, sigma2))

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

# for simplicity, we are using the true parameters of
# theta, sigma and pi in calculating our probabilities
true.stats <- component_stats(truth$data)
thetas <- true.stats$mean
sigmas <- true.stats$sd
pis <- true.stats$p

# check the two CN combinations we are interested in exists in our simulate data
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
genotype.df$offt <- thetas[(genotype.df$cn.off+1)]
genotype.df$mots <- sigmas[(genotype.df$cn.mot+1)]
genotype.df$fats <- sigmas[(genotype.df$cn.fat+1)]
genotype.df$offs <- sigmas[(genotype.df$cn.off+1)]
genotype.df$motp <- pis[(genotype.df$cn.mot+1)]
genotype.df$fatp <- pis[(genotype.df$cn.fat+1)]
genotype.df$offp <- pis[(genotype.df$cn.off+1)]
mendel <- rbinom(n=nrow(cn.df), 1, prob=0.9)
mendel.count <- sum(mendel)

# sort out the data
inten.mat <- data.frame(matrix(ncol=3, data=truth$data$log_ratio))
colnames(inten.mat) <- c("cmot", "cfat", "coff")

# compare row 5 vs 29
# We use simulated LRR data from above to generate our probability table
inten.mat1 <- inten.mat[5,]
inten.mat2 <- inten.mat[18,]

# we use inten.mat2 which is -0.155, 0.263 and -2.54 for
# mum, dad and child respectively
# substitute any other values into inten.mat3 to run other values
inten.mat3 <- inten.mat2
#inten.mat3 <- inten.mat1

# generate the table
genotype.df$cmot.dens <- dnorm(inten.mat3$cmot, mean=genotype.df$mott, sd=(genotype.df$mots)^0.5)
genotype.df$cfat.dens <- dnorm(inten.mat3$cfat, mean=genotype.df$fatt, sd=(genotype.df$fats)^0.5)
genotype.df$coff.dens <- dnorm(inten.mat3$coff, mean=genotype.df$offt, sd=(genotype.df$offs)^0.5)

# note this bit is how I am handling the Mendelian indicator
# the counts are based off the simulated dataset above
eta <- rbeta(1, 1+mendel.count, 1+(N-mendel.count))
mendel <- rbinom(n=nrow(cn.df), 1, prob=eta)
mendel.count <- sum(mendel)

#fix eta
eta <-1

mprob.unlist <- mprob[,c(1:3)]
trio.prob <- c(t(mprob.unlist))
genotype.df$tp <- trio.prob
genotype.df$mend <- genotype.df$coff.dens/trio.prob
genotype.df$mend[is.infinite(genotype.df$mend)]<-0
genotype.df$mend.denom <- (genotype.df$mend + genotype.df$cmot.dens + genotype.df$cfat.dens)
genotype.df$mend.num <- (genotype.df$mend * genotype.df$cmot.dens * genotype.df$cfat.dens) 
genotype.df$m <- genotype.df$mend.num / genotype.df$mend.denom

genotype.df$m0 <- genotype.df$motp * genotype.df$fatp * genotype.df$offp * (1-genotype.df$m)
genotype.df$m1 <- genotype.df$motp * genotype.df$fatp * genotype.df$tp * genotype.df$m
genotype.df$nump <- genotype.df$cmot.dens * genotype.df$cfat.dens * genotype.df$coff.dens * (genotype.df$m0 + genotype.df$m1)
geno.denom <- sum(genotype.df$nump)
genotype.df$prob <- genotype.df$nump / geno.denom

# update genotypes to the trio
row.select <- sample(1:27, 1, prob=genotype.df$prob)
genotype.assign <- as.numeric(genotype.df[row.select,c(1:3)])

# relabel the table to be consistent with whiteboard
genotype.probs <- genotype.df[,c(1:3, 13:15,21,22,24)]
colnames(genotype.probs) <- c("cn.mum", "cn.dad", "cn.child", "A", "B", "C", "D", "E", "G")

genotype.probs$abc <- genotype.probs$A*genotype.probs$B*genotype.probs$C
# the probability for genotype 2,1,1
genotype.probs[23,9]

# the probability for genotype 2,1,0
genotype.probs[22,9]

mat <- unite(genotype.probs, "m,f", c("cn.mum", "cn.dad"), sep=",") %>%
  unite("m,f,o", c("m,f", "cn.child"), sep=",") %>%
  as.tibble %>%
  mutate(`m,f,o`=factor(`m,f,o`))
ggplot(mat, aes(`m,f,o`, G)) +
  geom_point()
ggsave("joint_prob6.pdf", width=12, height=6)

###############################
###below is test gibbs sampler
#############################

# the actual MCMC - here we test inten.mat2 (CN2,1,0)
genotype.mat <- matrix(nrow=20, ncol=3)
genotype.true <- vector(length=20)
genotype.mend <- vector(length=20)

for (i in 1:20){
  
  genotype.df$cmot.dens <- dnorm(inten.mat2$cmot, mean=genotype.df$mott, sd=(genotype.df$mots)^0.5)
  genotype.df$cfat.dens <- dnorm(inten.mat2$cfat, mean=genotype.df$fatt, sd=(genotype.df$fats)^0.5)
  genotype.df$coff.dens <- dnorm(inten.mat2$coff, mean=genotype.df$offt, sd=(genotype.df$offs)^0.5)
  mprob.unlist <- mprob[,c(1:3)]
  trio.prob <- c(t(mprob.unlist))
  genotype.df$tp <- trio.prob
  genotype.df$mend <- genotype.df$coff.dens/trio.prob
  genotype.df$mend[is.infinite(genotype.df$mend)]<-0
  genotype.df$mend.denom <- (genotype.df$mend + genotype.df$cmot.dens + genotype.df$cfat.dens)
  genotype.df$mend.num <- (genotype.df$mend * genotype.df$cmot.dens * genotype.df$cfat.dens) 
  genotype.df$m <- genotype.df$mend.num / genotype.df$mend.denom
  
  genotype.df$m0 <- genotype.df$motp * genotype.df$fatp * genotype.df$offp * (1-genotype.df$m)
  genotype.df$m1 <- genotype.df$motp * genotype.df$fatp * genotype.df$tp * genotype.df$m
  genotype.df$nump <- genotype.df$cmot.dens * genotype.df$cfat.dens * genotype.df$coff.dens * (genotype.df$m0 + genotype.df$m1)
  geno.denom <- sum(genotype.df$nump)
  genotype.df$prob <- genotype.df$nump / geno.denom
  
  # update genotypes to the trio
  row.select <- sample(1:27, 1, prob=genotype.df$prob)
  genotype.assign <- as.numeric(genotype.df[row.select,c(1:3)])
  
  genotype.mat[i,] <- genotype.assign
  genotype.mend[i] <- genotype.df$prob[13]
  genotype.true[i] <- genotype.df$prob[16]
}

mean(genotype.mend)
mean(genotype.true)
