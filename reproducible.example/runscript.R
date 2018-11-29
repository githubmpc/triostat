#! /usr/bin/env Rscript

library(devtools)
#library(Rcpp)
library(coda)
library(tidyverse)
library(bindrcpp)

# directory where CNPBayes and triostats is located
setwd("/gpfs/data/chaklab/home/choum02")
load_all("CNPBayes_trios")
load_all("triostats")

# directory where parameter files are located
#setwd("/gpfs/data/chaklab/home/choum02/batchrun27")
setwd("~/Desktop/Chakravarti_Lab/git/testing scripts/reproducible.example")

# choose your parameter file
params.all <- readRDS("sim.params.rds")


# commandArgs specified as array in batch script
ab <- commandArgs(trailingOnly=TRUE)  %>%
  as.integer
ab <- ab

#seeds <- readRDS(file.path("./params_del", "params_seeds.rds"))
seeds <- readRDS("params_seeds.rds")
seed <- seeds[ab]
set.seed(seed)

##--------------------------------------------------
##
message("cnpbayes_trios")
##
##--------------------------------------------------

# note maplabel and hp defined manually here for now
maplabel <- c(0,1,2)
k <- length(maplabel)
hp <- HyperparametersTrios(k = 3)
nbatch <- 1

truth <- simulateData(params.all, ab, maplabel, nbatch)

mp <- McmcParams(iter=1500, burnin=1500, thin=1, nStarts=3)
model <- gibbs_trios(model="TBM", dat=as.tibble(truth$data),
                     batches=truth$data$batches,
                     mp=mp, mprob=mprob, maplabel=maplabel,
                     k_range=c(3, 3),
                     max_burnin=6000)

mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
                            batches=batch(model[[1]]),
                            mp=mp,
                            k_range=c(3, 3),
                            max_burnin=6000)

#mb2 <- gibbs(model="MB", dat=truth$data$log_ratio,
#             batches=rep(c(1:nbatch),
#                         length.out = 3*N),
#             mp=mp, k_range=c(k, k), max_burnin=2000)

trio.results <- sum.fit(model, truth, maplabel)
mb.results <- sum.fit(mb2, truth, maplabel)

# other random stats
true.stats <- component_stats(truth$data)
truth.parent.pi <- table(trio.results$cn.truth.parent)/length(trio.results$cn.truth.parent)
truth.child.pi <- table(trio.results$cn.truth.child)/length(trio.results$cn.truth.child)
child.pi <- table(trio.results$results.child)/length(trio.results$results.child)
parent.pi <- table(trio.results$results.parent)/length(trio.results$results.parent)
logrr <- truth$data$log_ratio

# input all the info
summaryResults <- list(params = params,
                       sim.params = truth$params,
                       N = length(logrr)/3,
                       SimLogRR = logrr,
                       SimTruth = true.stats,
                       SimPi = true.stats$p,
                       SimParPi = truth.parent.pi,
                       SimChildPi = truth.child.pi,
                       SimThetas = true.stats$mean,
                       SimVar = (true.stats$sd)^2,
                       TruthCNcall = truth$data$copy_number,
                       TruthCNpar = trio.results$cn.truth.parent,
                       TruthCNoff = trio.results$cn.truth.child,
                       TrioCNcall = trio.results$results@z,
                       TrioCNpar = trio.results$results.parent,
                       TrioCNoff = trio.results$results.child,
                       MBCNcall = mb.results$results@z,
                       MBCNpar = mb.results$results.parent,
                       MBCNoff = mb.results$results.child,
                       Accuracy = trio.results$prop.true.overall,
                       AccuracyParents = trio.results$prop.true.parent,
                       AccuracyChild = trio.results$prop.true.child,
                       AccuracyMB = mb.results$prop.true.overall,
                       AccuracyMBParents = mb.results$prop.true.parent,
                       AccuracyMBChild = mb.results$prop.true.child,
                       ModelPi = model[[1]]@modes$mixprob,
                       ModelParentPi = parent.pi,
                       ModelChildPi = child.pi,
                       ModelTheta = model[[1]]@modes$theta,
                       ModelSigma2 = model[[1]]@modes$sigma2
)

## save to some output directory
pathout <- file.path("/gpfs/data/chaklab/home/choum02/batchrun27/results")
#pathout <- file.path("~/Desktop/Chakravarti_Lab/git")
results.out <- paste0("params_", ab, ".rds")
#saveRDS(summary.results, file=file.path(pathout, params.rds))
setwd("/gpfs/data/chaklab/home/choum02/batchrun27/results")
saveRDS(summaryResults, file=file.path(pathout, results.out))
