simulateData <- function(params.file, ab, maplabel, nbatch) {
  paramN <- nrow(params.all)
  params.ab <- ab %% paramN

  params <- unlist(params.all[params.ab,c(1:9)])
  N <- params.all[params.ab,10]
  params <- data.frame(matrix(data=params, nrow=3,ncol=3))
  colnames(params) <- c("p", "theta", "sigma2")


  mprob <- mprob.matrix(tau=c(0.5, 0.5, 0.5), maplabel, error=0)
  truth <- simulate_data_multi2(params, N=N,
                                batches = rep(c(1:nbatch),
                                              length.out = 3*N),
                                error=0, mprob, maplabel)

truth

}

sum.fit <- function(model, truth, maplabel){
  is_offspring <- truth$data$family_member=="o"
  cn_offspring <- truth$data$copy_number[is_offspring]

  results <- z2cn(model[[1]], maplabel)
  is_offspring <- model[[1]]@triodata$family_member=="o"
  cn.truth <- truth$data$copy_number
  cn.truth.parent <- cn.truth[!is_offspring]
  cn.truth.child <- cn.truth[is_offspring]
  results.parent <- results@z[!is_offspring]
  results.child <- results@z[is_offspring]

  prop.true.overall <- sum(results@z == cn.truth, na.rm=T) / length(cn.truth)
  prop.true.parent <- sum(results.parent == cn.truth.parent, na.rm=T) / length(cn.truth.parent)
  prop.true.child <- sum(results.child == cn.truth.child, na.rm=T) / length(cn.truth.child)

  sumFit <- list(results = results,
                 cn.truth = cn.truth,
                 cn.truth.parent = cn.truth.parent,
                 cn.truth.child = cn.truth.child,
                 results.parent = results.parent,
                 results.child = results.child,
                 prop.true.overall = prop.true.overall,
                 prop.true.parent = prop.true.parent,
                 prop.true.child = prop.true.child)
  sumFit
}

simdata.read <- function(results.file){
  n <- results.file$N
  id.index <- formatC(seq_len(n), flag="0", width=3)
  motdex <- seq(1, 3*n, by=3)
  fatdex <- seq(2, 3*n, by=3)
  offdex <- seq(3, 3*n, by=3)
  results.tbl1 <- tibble(m=results.file$TrioCNcall[motdex],
                        f=results.file$TrioCNcall[fatdex],
                        o=results.file$TrioCNcall[offdex],
                        id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="trio_calls", -id)
  results.tbl2 <- tibble(m=results.file$MBCNcall[motdex],
                        f=results.file$MBCNcall[fatdex],
                        o=results.file$MBCNcall[offdex],
                        id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="ind_calls", -id)
  results.tbl3 <- tibble(m=results.file$TruthCNcall[motdex],
                         f=results.file$TruthCNcall[fatdex],
                         o=results.file$TruthCNcall[offdex],
                         id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="truth_calls", -id)
   results.tbl <- left_join(results.tbl1, results.tbl2, by=c("id", "family_member")) %>%
    mutate(family_member=factor(family_member, levels=c("m", "f", "o"))) %>%
    arrange(id, family_member)
   results.tbl <- left_join(results.tbl, results.tbl3, by=c("id", "family_member")) %>%
     mutate(family_member=factor(family_member, levels=c("m", "f", "o"))) %>%
     arrange(id, family_member)
  results.tbl
}

summarise_results <- function(results){
  trios <- results %>% select(c(id, family_member, trio_calls)) %>%
    spread(family_member, trio_calls)

  ind <- results %>% select(c(id, family_member, ind_calls)) %>%
    spread(family_member, ind_calls) %>%
    set_colnames(c("id", "i.m", "i.f", "i.o"))

  truec <- results %>% select(c(id, family_member, truth_calls)) %>%
    spread(family_member, truth_calls) %>%
    set_colnames(c("id", "t.m", "t.f", "t.o"))

  tab <- left_join(trios, ind, by=c("id")) %>%
          left_join(truec, by=c("id"))

  unite2 <- select(tab, c(t.m, t.f)) %>%
    mutate(parents = apply(., 1, function(x) paste(sort(x), collapse="_"))) %>%
    mutate(id = truec$id) %>%
    select(-c(t.m,t.f))

  unite2$parents <- paste0("geno",unite2$parents)
  stables <- left_join(tab, unite2, by=c("id"))
  parent.levels <- expand.grid(0:2,0:2) %>% filter(Var1<=Var2) %>%
    mutate(parent.levels = paste0("geno", Var1, "_", Var2),
           parent.levels = factor(parent.levels)) %$%
    parent.levels
  stables$parents <- factor(stables$parents, levels = parent.levels)
  stables
}

overall_metric <- function(stables){
  correct <- stables %>% filter(o==t.o) %>% nrow
  total <- stables %>% nrow
  accuracy <- correct / total
  overall.metric <- list(c(correct = correct,
                    total = total,
                    accuracy = accuracy))
  overall <- list(overall = overall.metric)
  overall
}

overall_indmetric <- function(stables){
  correct <- stables %>% filter(i.o==t.o) %>% nrow
  total <- stables %>% nrow
  accuracy <- correct / total
  overall.metric <- list(c(correct = correct,
                           total = total,
                           accuracy = accuracy))
  overall <- list(ind.overall = overall.metric)
  overall
}

summarise_metrichild <- function(stables){
  cn2.correct <- stables %>% filter(t.o==2 & o==t.o) %>% nrow
  cn2.total <- stables %>% filter(t.o==2) %>% nrow
  cn2.accuracy <- cn2.correct / cn2.total
  cn1.correct <- stables %>% filter(t.o==1 & o==t.o) %>% nrow
  cn1.total <- stables %>% filter(t.o==1) %>% nrow
  cn1.accuracy <- cn1.correct / cn1.total
  cn0.correct <- stables %>% filter(t.o==0 & o==t.o) %>% nrow
  cn0.total <- stables %>% filter(t.o==0) %>% nrow
  cn0.accuracy <- cn0.correct / cn0.total

  cn2.metric <- list(c(correct = cn2.correct,
                     total = cn2.total,
                     accuracy = cn2.accuracy))
  cn1.metric <- list(c(correct = cn1.correct,
                     total = cn1.total,
                     accuracy = cn1.accuracy))
  cn0.metric <- list(c(correct = cn0.correct,
                     total = cn0.total,
                     accuracy = cn0.accuracy))

  summary.metric <- list(cn2 = cn2.metric,
                         cn1 = cn1.metric,
                         cn0 = cn0.metric)
  summary.metric
}

summarise_metrichildind <- function(stables){
  cn2.correct <- stables %>% filter(t.o==2 & i.o==t.o) %>% nrow
  cn2.total <- stables %>% filter(t.o==2) %>% nrow
  cn2.accuracy <- cn2.correct / cn2.total
  cn1.correct <- stables %>% filter(t.o==1 & i.o==t.o) %>% nrow
  cn1.total <- stables %>% filter(t.o==1) %>% nrow
  cn1.accuracy <- cn1.correct / cn1.total
  cn0.correct <- stables %>% filter(t.o==0 & i.o==t.o) %>% nrow
  cn0.total <- stables %>% filter(t.o==0) %>% nrow
  cn0.accuracy <- cn0.correct / cn0.total

  cn2.metric <- list(c(correct = cn2.correct,
                       total = cn2.total,
                       accuracy = cn2.accuracy))
  cn1.metric <- list(c(correct = cn1.correct,
                       total = cn1.total,
                       accuracy = cn1.accuracy))
  cn0.metric <- list(c(correct = cn0.correct,
                       total = cn0.total,
                       accuracy = cn0.accuracy))

  summary.metric <- list(ind.cn2 = cn2.metric,
                         ind.cn1 = cn1.metric,
                         ind.cn0 = cn0.metric)
  summary.metric
}

# by parent correct/ incorrect
filter_parent <- function(stables){
  parents.correct <- stables %>% filter(m==t.m & f==t.f)
  parents.incorrect <- stables %>% filter(m!=t.m | f!=t.f)

  filtered.results <- list(parents.correct = parents.correct,
                           parents.incorrect = parents.incorrect)
  filtered.results
}

filter_indparent <- function(stables){
  parents.correct <- stables %>% filter(i.m==t.m & i.f==t.f)
  parents.incorrect <- stables %>% filter(i.m!=t.m | i.f!=t.f)

  filtered.results <- list(ind.parents.correct = parents.correct,
                           ind.parents.incorrect = parents.incorrect)
  filtered.results
}

flattenlist <- function(x){
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){
    Recall(out)
  }else{
    return(out)
  }
}

# by combined parental genotype
filter_parentgeno <- function(stables){
  stab.list <- split(stables, stables$parents)
  stab.list
}

# by combined parental genotype and parent correct/incorrect
filter_doubleparent <- function(stables){
  stab.list <- split(stables, stables$parents)
  stab.list2 <- lapply(stab.list, filter_parent)
  stab.list2 <- flattenlist(stab.list2)
  stab.list2
}

filter_indoubleparent <- function(stables){
  stab.list <- split(stables, stables$parents)
  stab.list2 <- lapply(stab.list, filter_indparent)
  stab.list2 <- flattenlist(stab.list2)
  stab.list2
}

#overall_results <- function(stables){
#  stab.list1 <- filter_parent(stables)
#  stab.list2 <- filter_indparent(stables)
#  stab.list3 <- filter_parentgeno(stables)
#  stab.list4 <- filter_doubleparent(stables)
#  stab.list5 <- filter_indoubleparent(stables)
#  stab.result1 <- overall_metric(stables)
#  stab.result2 <- overall_indmetric(stables)
#  stab.result3 <- summarise_metrichild(stables)
#  stab.result4 <- summarise_metrichildind(stables)
#  stab.result5 <- lapply(stab.list1, summarise_metrichild)
#  stab.result6 <- lapply(stab.list2, summarise_metrichildind)
#  stab.result7 <- lapply(stab.list3, summarise_metrichild)
#  stab.result8 <- lapply(stab.list3, summarise_metrichildind)
#  stab.result9 <- lapply(stab.list4, summarise_metrichild)
#  stab.result10 <- lapply(stab.list5, summarise_metrichildind)
#  stab.results <- list(stab.result1, stab.result2, stab.result3,
#                       stab.result4, stab.result5,
#                       stab.result6, stab.result7, stab.result8,
#                       stab.result9, stab.result10)
#  stab.results <- flattenlist(stab.results)
#  stab.results <- data.frame(t(data.frame(stab.results)))
# #stab.results <- rownames_to_column(stab.results)
#  stab.results
#}

results_reduce <- function(results.list){
  results.table <- Reduce('+', results.list)
  results.table <- rownames_to_column(results.table)
  results.table <- results.table %>% select(-c(accuracy)) %>%
    mutate(accuracy = correct/total)
  results.table
}

other_metrics <- function(results.list){
  results.list.other <- aaply(laply(results.list, as.matrix), c(2, 3), quantile, na.rm=T)
  results.25 <- data.frame(results.list.other[,,2][,3])
  colnames(results.25) <- "quartile1"
  results.75 <- data.frame(results.list.other[,,4][,3])
  colnames(results.75) <- "quartile3"
  results.list.sd <- aaply(laply(results.list, as.matrix), c(2, 3), sd, na.rm=T)
  results.sd <- data.frame(results.list.sd[,3])
  colnames(results.sd) <- "sd"
  results.other <- cbind(results.25, results.75, results.sd)
  results.other
}

output_results <- function(results.list){
  results.table <- results_reduce(results.list)
  results.other <- other_metrics(results.list)
  rtable <- cbind(results.table, results.other)
  rtable
}


overall_results2 <- function(stables){
  stab.list1 <- filter_parent(stables)
  stab.list2 <- filter_indparent(stables)
  stab.list3 <- filter_parentgeno(stables)
  stab.list4 <- filter_doubleparent(stables)
  stab.list5 <- filter_indoubleparent(stables)
  stab.result1 <- overall_metric(stables)
  stab.result2 <- overall_indmetric(stables)
  stab.result3 <- summarise_metrichild(stables)
  stab.result4 <- summarise_metrichildind(stables)
  stab.result5 <- lapply(stab.list1, summarise_metrichild)
  stab.result6 <- lapply(stab.list2, summarise_metrichildind)
  stab.result7 <- lapply(stab.list3, summarise_metrichild)
  stab.result8 <- lapply(stab.list3, summarise_metrichildind)
  stab.result9 <- lapply(stab.list4, summarise_metrichild)
  stab.result10 <- lapply(stab.list5, summarise_metrichildind)
  stab.results1 <- list(stab.result1, stab.result3,
                        stab.result5, stab.result7,
                        stab.result9)
  stab.results2 <- list(stab.result2, stab.result4,
                        stab.result6, stab.result8,
                        stab.result10)
  stab.results1 <- flattenlist(stab.results1)
  stab.results2 <- flattenlist(stab.results2)
  stab.results1 <- data.frame(t(data.frame(stab.results1)))
  stab.results2 <- data.frame(t(data.frame(stab.results2)))
  stab.results <- cbind(stab.results1, stab.results2)
  colnames(stab.results) <- c("correct.trio", "total.trio", "accuracy.trio", "correct.ind", "total.ind", "accuracy.ind")
  stab.results
}

results_reduce2 <- function(results.list){
  results.table <- Reduce('+', results.list)
  results.table <- rownames_to_column(results.table)
  results.table <- results.table %>% select(-c(accuracy.trio)) %>%
    mutate(accuracy.trio = correct.trio/total.trio)
  results.table <- results.table %>% select(-c(accuracy.ind)) %>%
    mutate(accuracy.ind = correct.ind/total.ind)
  results.table
}

other_metrics2 <- function(results.list){
  results.list.other <- aaply(laply(results.list, as.matrix), c(2, 3), quantile, na.rm=T)
  results.25 <- data.frame(results.list.other[,,2][,3])
  colnames(results.25) <- "quartile1"
  results.75 <- data.frame(results.list.other[,,4][,3])
  colnames(results.75) <- "quartile3"
  results.list.sd <- aaply(laply(results.list, as.matrix), c(2, 3), sd, na.rm=T)
  results.sd <- data.frame(results.list.sd[,3])
  colnames(results.sd) <- "sd"
  results.25ind <- data.frame(results.list.other[,,2][,6])
  colnames(results.25ind) <- "quartile1.ind"
  results.75ind <- data.frame(results.list.other[,,4][,6])
  colnames(results.75ind) <- "quartile3.ind"
  results.list.sd <- aaply(laply(results.list, as.matrix), c(2, 3), sd, na.rm=T)
  results.sdind <- data.frame(results.list.sd[,6])
  colnames(results.sdind) <- "sd.ind"
  results.other <- cbind(results.25, results.75, results.sd, results.25ind,
                         results.75ind, results.sdind)
  results.other
}

output_results2 <- function(results.list){
  results.table <- results_reduce2(results.list)
  results.other <- other_metrics2(results.list)
  rtable <- cbind(results.table, results.other)
  rtable
}

rtable_overall <- function(rtable){
  rtable.overall <- rtable[c(1:4),c(2,3,6,10,4,5,7,13)]
  rtable.overall <- round(rtable.overall,3)
  #row.names(rtable.overall) <- c("Overall","CN2", "CN1", "CN0")
  colnames(rtable.overall) <- c("correct.trio", "total.trio", "accuracy.trio", "sd.trio", "correct.ind", "total.ind", "accuracy.ind", "sd.ind")
  rtable.overall
}

rtable_overalls <- function(rtable){
  rtable.overall <- rtable[c(1:4),c(6,10,7,13)]
  rtable.overall <- round(rtable.overall,3)
  #row.names(rtable.overall) <- c("Overall","CN2", "CN1", "CN0")
  colnames(rtable.overall) <- c("accuracy.trio", "sd.trio", "accuracy.ind", "sd.ind")
  rtable.overall
}

rtable_af <- function(rtable1, rtable2, rtable3){
  #rtable.af <- data.frame(matrix(nrow=9, ncol=8))
  rtable.fill1 <- rtable1[c(2:4),c(2,3,6,10,4,5,7,13)]
  rtable.fill2 <- rtable2[c(2:4),c(2,3,6,10,4,5,7,13)]
  rtable.fill3 <- rtable3[c(2:4),c(2,3,6,10,4,5,7,13)]
  rtable.fill1 <- round(rtable.fill1, 3)
  rtable.fill2 <- round(rtable.fill2, 3)
  rtable.fill3 <- round(rtable.fill3, 3)
  rtable.af <- rbind(rtable.fill1, rtable.fill2, rtable.fill3)

  #row.names(rtable.af) <- c("p0.5_CN2", "p0.3_CN2", "p0.1_CN2", "p0.5_CN1", "p0.3_CN1", "p0.1_CN1", "p0.5_CN0", "p0.3_CN0", "p0.1_CN0")
  row.names(rtable.af) <- c("p0.5_CN2", "p0.5_CN1", "p0.5_CN0", "p0.3_CN2",
                            "p0.3_CN1", "p0.3_CN0", "p0.1_CN2", "p0.1_CN1",
                            "p0.1_CN0")
  colnames(rtable.af) <- c("correct.trio", "total.trio", "accuracy.trio", "sd.trio", "correct.ind", "total.ind", "accuracy.ind", "sd.ind")
  rtable.af
}

rtable_afs <- function(rtable1, rtable2, rtable3){
  #rtable.af <- data.frame(matrix(nrow=9, ncol=8))
  rtable.fill1 <- rtable1[c(2:4),c(6,10,7,13)]
  rtable.fill2 <- rtable2[c(2:4),c(6,10,7,13)]
  rtable.fill3 <- rtable3[c(2:4),c(6,10,7,13)]
  rtable.fill1 <- round(rtable.fill1, 3)
  rtable.fill2 <- round(rtable.fill2, 3)
  rtable.fill3 <- round(rtable.fill3, 3)
  rtable.af <- rbind(rtable.fill1, rtable.fill2, rtable.fill3)

  #row.names(rtable.af) <- c("p0.5_CN2", "p0.3_CN2", "p0.1_CN2", "p0.5_CN1", "p0.3_CN1", "p0.1_CN1", "p0.5_CN0", "p0.3_CN0", "p0.1_CN0")
  row.names(rtable.af) <- c("p0.5_CN2", "p0.5_CN1", "p0.5_CN0", "p0.3_CN2",
                            "p0.3_CN1", "p0.3_CN0", "p0.1_CN2", "p0.1_CN1",
                            "p0.1_CN0")
  colnames(rtable.af) <- c("accuracy.trio", "sd.trio", "accuracy.ind", "sd.ind")
  rtable.af
}

rtable_pc <- function(rtable1, rtable2, rtable3){
  #rtable.af <- data.frame(matrix(nrow=9, ncol=8))
  rtable.fill1 <- rtable1[c(5:7),c(2,3,6,10,4,5,7,13)]
  rtable.fill2 <- rtable2[c(5:7),c(2,3,6,10,4,5,7,13)]
  rtable.fill3 <- rtable3[c(5:7),c(2,3,6,10,4,5,7,13)]
  rtable.fill1 <- round(rtable.fill1, 3)
  rtable.fill2 <- round(rtable.fill2, 3)
  rtable.fill3 <- round(rtable.fill3, 3)
  rtable.af <- rbind(rtable.fill1, rtable.fill2, rtable.fill3)

  #row.names(rtable.af) <- c("p0.5_CN2", "p0.3_CN2", "p0.1_CN2", "p0.5_CN1", "p0.3_CN1", "p0.1_CN1", "p0.5_CN0", "p0.3_CN0", "p0.1_CN0")
  row.names(rtable.af) <- c("p0.5_CN2", "p0.5_CN1", "p0.5_CN0", "p0.3_CN2",
                            "p0.3_CN1", "p0.3_CN0", "p0.1_CN2", "p0.1_CN1",
                            "p0.1_CN0")
  colnames(rtable.af) <- c("correct.trio", "total.trio", "accuracy.trio", "sd.trio", "correct.ind", "total.ind", "accuracy.ind", "sd.ind")
  rtable.af
}

rtable_pcs <- function(rtable1, rtable2, rtable3){
  #rtable.af <- data.frame(matrix(nrow=9, ncol=8))
  rtable.fill1 <- rtable1[c(5:7),c(6,10,7,13)]
  rtable.fill2 <- rtable2[c(5:7),c(6,10,7,13)]
  rtable.fill3 <- rtable3[c(5:7),c(6,10,7,13)]
  rtable.fill1 <- round(rtable.fill1, 3)
  rtable.fill2 <- round(rtable.fill2, 3)
  rtable.fill3 <- round(rtable.fill3, 3)
  rtable.af <- rbind(rtable.fill1, rtable.fill2, rtable.fill3)

  #row.names(rtable.af) <- c("p0.5_CN2", "p0.3_CN2", "p0.1_CN2", "p0.5_CN1", "p0.3_CN1", "p0.1_CN1", "p0.5_CN0", "p0.3_CN0", "p0.1_CN0")
  row.names(rtable.af) <- c("p0.5_CN2", "p0.5_CN1", "p0.5_CN0", "p0.3_CN2",
                            "p0.3_CN1", "p0.3_CN0", "p0.1_CN2", "p0.1_CN1",
                            "p0.1_CN0")
  colnames(rtable.af) <- c("accuracy.trio", "sd.trio", "accuracy.ind", "sd.ind")
  rtable.af
}

rtable_pic <- function(rtable1, rtable2, rtable3){
  #rtable.af <- data.frame(matrix(nrow=9, ncol=8))
  rtable.fill1 <- rtable1[c(8:10),c(2,3,6,10,4,5,7,13)]
  rtable.fill2 <- rtable2[c(8:10),c(2,3,6,10,4,5,7,13)]
  rtable.fill3 <- rtable3[c(8:10),c(2,3,6,10,4,5,7,13)]
  rtable.fill1 <- round(rtable.fill1, 3)
  rtable.fill2 <- round(rtable.fill2, 3)
  rtable.fill3 <- round(rtable.fill3, 3)
  rtable.af <- rbind(rtable.fill1, rtable.fill2, rtable.fill3)

  #row.names(rtable.af) <- c("p0.5_CN2", "p0.3_CN2", "p0.1_CN2", "p0.5_CN1", "p0.3_CN1", "p0.1_CN1", "p0.5_CN0", "p0.3_CN0", "p0.1_CN0")
  row.names(rtable.af) <- c("p0.5_CN2", "p0.5_CN1", "p0.5_CN0", "p0.3_CN2",
                            "p0.3_CN1", "p0.3_CN0", "p0.1_CN2", "p0.1_CN1",
                            "p0.1_CN0")
  colnames(rtable.af) <- c("correct.trio", "total.trio", "accuracy.trio", "sd.trio", "correct.ind", "total.ind", "accuracy.ind", "sd.ind")
  rtable.af
}

rtable_pics <- function(rtable1, rtable2, rtable3){
  #rtable.af <- data.frame(matrix(nrow=9, ncol=8))
  rtable.fill1 <- rtable1[c(8:10),c(6,10,7,13)]
  rtable.fill2 <- rtable2[c(8:10),c(6,10,7,13)]
  rtable.fill3 <- rtable3[c(8:10),c(6,10,7,13)]
  rtable.fill1 <- round(rtable.fill1, 3)
  rtable.fill2 <- round(rtable.fill2, 3)
  rtable.fill3 <- round(rtable.fill3, 3)
  rtable.af <- rbind(rtable.fill1, rtable.fill2, rtable.fill3)

  #row.names(rtable.af) <- c("p0.5_CN2", "p0.3_CN2", "p0.1_CN2", "p0.5_CN1", "p0.3_CN1", "p0.1_CN1", "p0.5_CN0", "p0.3_CN0", "p0.1_CN0")
  row.names(rtable.af) <- c("p0.5_CN2", "p0.5_CN1", "p0.5_CN0", "p0.3_CN2",
                            "p0.3_CN1", "p0.3_CN0", "p0.1_CN2", "p0.1_CN1",
                            "p0.1_CN0")
  colnames(rtable.af) <- c("accuracy.trio", "sd.trio", "accuracy.ind", "sd.ind")
  rtable.af
}

rtable_geno <- function(rtable){
  #rtable.geno <- data.frame(matrix(nrow=18, ncol=8))
  rtable.fill1 <- rtable[c(11:28),c(2,3,6,10,4,5,7,13)]
  rtable.geno <- round(rtable.fill1, 3)

  row.names(rtable.geno) <- c("geno00_CN2", "geno00_CN1", "geno00_CN0", "geno01_CN2", "geno01_CN1", "geno01_CN0",
                              "geno11_CN2", "geno11_CN1", "geno11_CN0", "geno02_CN2", "geno02_CN1", "geno02_CN0",
                              "geno12_CN2", "geno12_CN1", "geno12_CN0", "geno22_CN2", "geno22_CN1", "geno22_CN0")
  colnames(rtable.geno) <- c("correct.trio", "total.trio", "accuracy.trio", "sd.trio", "correct.ind", "total.ind", "accuracy.ind", "sd.ind")
  rtable.geno
}

rtable_genos <- function(rtable){
  #rtable.geno <- data.frame(matrix(nrow=18, ncol=8))
  rtable.fill1 <- rtable[c(11:28),c(6,10,7,13)]
  rtable.geno <- round(rtable.fill1, 3)

  row.names(rtable.geno) <- c("geno00_CN2", "geno00_CN1", "geno00_CN0", "geno01_CN2", "geno01_CN1", "geno01_CN0",
                              "geno11_CN2", "geno11_CN1", "geno11_CN0", "geno02_CN2", "geno02_CN1", "geno02_CN0",
                              "geno12_CN2", "geno12_CN1", "geno12_CN0", "geno22_CN2", "geno22_CN1", "geno22_CN0")
  colnames(rtable.geno) <- c("accuracy.trio", "sd.trio", "accuracy.ind", "sd.ind")
  rtable.geno
}

#  reps is the number of replicates
rtable_plotprep <- function(results.list, reps){
  results.listr <- lapply(results.list, rownames_to_column)
  rtable1mass <- bind_rows(results.listr, .id="grp")

  af.label <- c("p=0.5", "p=0.5", "p=0.5", "p=0.5",
                "p=0.3", "p=0.3", "p=0.3", "p=0.3",
                "p=0.1", "p=0.1", "p=0.1", "p=0.1")
  af.labels <- rep(af.label, reps)
  af.labels <- data.frame(af.labels)
  len <- length(results.list)
  af.label2 <- as.character(seq(1,len))
  size.label <- c("N50", "N100", "N250", "N1000")
  size.labels <- rep(size.label, reps)
  rlabels <- cbind(af.labels, size.labels, af.label2)

  rtable1m2 <- left_join(rtable1mass, rlabels, by=c("grp"="af.label2"))
  rtable1m2$geno <- str_sub(rtable1m2$rowname,-3,-1)
  colnames(rtable1m2)[3:8] <- c("correct.trio", "total.trio", "accuracy.trio", "correct.ind", "total.ind", "accuracy.ind")
  rtable1m2$accur.dif <- rtable1m2$accuracy.trio - rtable1m2$accuracy.ind
  rtable1m2
}

rtable_sg <- function(rtable.prep){
  len <- nrow(rtable.prep)
  extract.no1 <- seq(11,len, by=64)
  extract.no2 <- seq(12,len, by=64)
  extract.no3 <- seq(13,len, by=64)
  extract.no4 <- seq(14,len, by=64)
  extract.no5 <- seq(15,len, by=64)
  extract.no6 <- seq(16,len, by=64)
  extract.no7 <- seq(17,len, by=64)
  extract.no8 <- seq(18,len, by=64)
  extract.no9 <- seq(19,len, by=64)
  extract.no10 <- seq(20,len, by=64)
  extract.no11 <- seq(21,len, by=64)
  extract.no12 <- seq(22,len, by=64)
  extract.no13 <- seq(23,len, by=64)
  extract.no14 <- seq(24,len, by=64)
  extract.no15 <- seq(25,len, by=64)
  extract.no16 <- seq(26,len, by=64)
  extract.no17 <- seq(27,len, by=64)
  extract.no18 <- seq(28,len, by=64)
  extract.no <- sort(c(extract.no1, extract.no2, extract.no3, extract.no4, extract.no5, extract.no6, extract.no7, extract.no8, extract.no9, extract.no10, extract.no11, extract.no12, extract.no13, extract.no14, extract.no15, extract.no16, extract.no17, extract.no18))

  rtable.sg <- rtable.prep[extract.no,]
  rtable.sg
}

#rtable_genopic <- function(rtable){
#  rtable.geno <- data.frame(matrix(nrow=18, ncol=8))
#  rtable.fill1 <- rtable[c(57:92),c(2,3,4,7)]
#  rtable.fill2 <- rtable[c(93:128),c(2,3,4,7)]
#  rtable.fill1 <- round(rtable.fill1, 4)
#  rtable.fill2 <- round(rtable.fill2, 4)
#  rtable.geno[1,c(1:4)] <- rtable.fill1[4,]
#  rtable.geno[1,c(5:8)] <- rtable.fill2[4,]
#  rtable.geno[2,c(1:4)] <- rtable.fill1[5,]
#  rtable.geno[2,c(5:8)] <- rtable.fill2[5,]
#  rtable.geno[3,c(1:4)] <- rtable.fill1[6,]
#  rtable.geno[3,c(5:8)] <- rtable.fill2[6,]
#  rtable.geno[4,c(1:4)] <- rtable.fill1[10,]
#  rtable.geno[4,c(5:8)] <- rtable.fill2[10,]
#  rtable.geno[5,c(1:4)] <- rtable.fill1[11,]
#  rtable.geno[5,c(5:8)] <- rtable.fill2[11,]
#  rtable.geno[6,c(1:4)] <- rtable.fill1[12,]
#  rtable.geno[6,c(5:8)] <- rtable.fill2[12,]
#  rtable.geno[7,c(1:4)] <- rtable.fill1[16,]
#  rtable.geno[7,c(5:8)] <- rtable.fill2[16,]
# rtable.geno[8,c(1:4)] <- rtable.fill1[17,]
#  rtable.geno[8,c(5:8)] <- rtable.fill2[17,]
# rtable.geno[9,c(1:4)] <- rtable.fill1[18,]
 # rtable.geno[9,c(5:8)] <- rtable.fill2[18,]
#  rtable.geno[10,c(1:4)] <- rtable.fill1[22,]
#  rtable.geno[10,c(5:8)] <- rtable.fill2[22,]
#  rtable.geno[11,c(1:4)] <- rtable.fill1[23,]
#  rtable.geno[11,c(5:8)] <- rtable.fill2[23,]
#  rtable.geno[12,c(1:4)] <- rtable.fill1[24,]
#  rtable.geno[12,c(5:8)] <- rtable.fill2[24,]
#  rtable.geno[13,c(1:4)] <- rtable.fill1[28,]
#  rtable.geno[13,c(5:8)] <- rtable.fill2[28,]
#  rtable.geno[14,c(1:4)] <- rtable.fill1[29,]
#  rtable.geno[14,c(5:8)] <- rtable.fill2[29,]
#  rtable.geno[15,c(1:4)] <- rtable.fill1[30,]
#  rtable.geno[15,c(5:8)] <- rtable.fill2[30,]
# rtable.geno[16,c(1:4)] <- rtable.fill1[34,]
#  rtable.geno[16,c(5:8)] <- rtable.fill2[34,]
#  rtable.geno[17,c(1:4)] <- rtable.fill1[35,]
#  rtable.geno[17,c(5:8)] <- rtable.fill2[35,]
#  rtable.geno[18,c(1:4)] <- rtable.fill1[36,]
# rtable.geno[18,c(5:8)] <- rtable.fill2[36,]

 # row.names(rtable.geno) <- c("geno00_CN2", "geno00_CN1", "geno00_CN0", "geno01_CN2", "geno01_CN1", "geno01_CN0",
  #                            "geno11_CN2", "geno11_CN1", "geno11_CN0", "geno02_CN2", "geno02_CN1", "geno02_CN0",
   #                           "geno12_CN2", "geno12_CN1", "geno12_CN0", "geno22_CN2", "geno22_CN1", "geno22_CN0")
  #colnames(rtable.geno) <- c("correct.trio", "total.trio", "accuracy.trio", "sd.trio", "correct.ind", "total.ind", "accuracy.ind", "sd.ind")
#  rtable.geno
#}

params.set<- function(sizes, af, theta, sigma){
  reps <- length(sizes) * length(af)
  thetar <- matrix(rep(theta, reps),ncol=length(theta),
                   nrow=reps, byrow=T)
  sigmar <- matrix(rep(sigma, reps),ncol=length(sigma),
                   nrow=reps, byrow=T)
  sizer <- rep(sizes, length(af))
  af <- data.frame(af)
  afr <- af[rep(seq_len(nrow(af)), each=length(sizes)),]
  params <- cbind(afr, sizer, thetar, sigmar)
  params
}

read.data <- function(file.path){

  for (i in 1:L){
    ia <- file.no[i]
    file.names <- paste0("params_",ia,".rds")
    params.in <- readRDS(file.names)
    p <- round(sqrt(params.in$params$p[1]),2)
    N <- params.in$N
    delta <- params.in$Delta
    accur <- params.in$Accuracy
    accur.par <- params.in$AccuracyParents
    accur.child <- params.in$AccuracyChild
    accur.mb.indep <- params.in$AccuracyMB
    accur.mb.par <- params.in$AccuracyMBParents
    accur.mb.child <- params.in$AccuracyMBChild

    truth.call.par <- params.in$TruthCNpar
    model.call.par <- params.in$TrioCNpar
    parentsep <- seq(1,length(truth.call.par), by=2)
    parentsep1 <- seq(2,length(truth.call.par), by=2)
    truth.call.par1 <- truth.call.par[parentsep]
    truth.call.par2 <- truth.call.par[parentsep1]
    model.call.par1 <- model.call.par[parentsep]
    model.call.par2 <- model.call.par[parentsep1]
    parcheck1 <- truth.call.par1 == model.call.par1
    parcheck2 <- truth.call.par2 == model.call.par2
    parallcorrect <- parcheck1 == parcheck2

    # back to more generalised counts
    #accur.xtab <- table(truth.call.par, model.call.par)
    #accur.default <- matrix(rep(0,9), nrow=3, ncol=3)
    #row.names(accur.default) <- c("0","1","2")
    #accur.default <- data.frame(accur.default)
    #accur.default <- rownames_to_column(accur.default)
    #accur.default$rowname <- as.numeric(accur.default$rowname)
    #indexcol <- as.numeric(row.names(accur.xtab))
    #accur.xtab2 <- cbind(indexcol, accur.xtab)
    #accur.xtab2 <- as.data.frame(accur.xtab2)
    #accur.xtab3 <- left_join(accur.default, accur.xtab2, by=c("rowname"="indexcol"))
    #cl <- ncol(accur.xtab3)
    #rn <- colnames(accur.xtab3)[5:cl]
    #accur.xtab31 <- data.frame(accur.xtab3[,5:cl])
    #colnames(accur.xtab31) <- rn

    #accur.default2 <- accur.default[,2:4]
    #colnames(accur.default2)[1:3] <- c(0,1,2)
    #accur.xtab4 <- smartbind(accur.default2, accur.xtab31)
    #accur.xtab4 <- accur.xtab4[-c(1:3),]

    # note this section is specific to deletion
    #accur.cn2 <- round(accur.xtab4[3,3] / sum(accur.xtab4[3,], na.rm=T), 4)
    #accur.cn1 <- round(accur.xtab4[2,2] / sum(accur.xtab4[2,], na.rm=T), 4)
    #accur.cn0 <- round(accur.xtab4[1,1] / sum(accur.xtab4[1,], na.rm=T), 4)

    accurpartab <- data.frame(cbind(model.call.par, truth.call.par))
    accurpartab2 <- accurpartab[accurpartab$truth.call.par==2,]
    cn2p <- sum(accurpartab2$model.call.par==accurpartab2$truth.call.par)
    cn2pt <- nrow(accurpartab2)
    accur.cn2 <- cn2p / cn2pt

    accurpartab2 <- accurpartab[accurpartab$truth.call.par==1,]
    cn1p <- sum(accurpartab2$model.call.par==accurpartab2$truth.call.par)
    cn1pt <- nrow(accurpartab2)
    accur.cn1 <- cn1p / cn1pt
    accurpartab2 <- accurpartab[accurpartab$truth.call.par==0,]
    cn0p <- sum(accurpartab2$model.call.par==accurpartab2$truth.call.par)
    cn0pt <- nrow(accurpartab2)
    accur.cn0 <- cn0p / cn0pt

    # create accuracy counts by parental genotype
    truth.par <- paste0(truth.call.par1, truth.call.par2)
    truth.par <- gsub("10", "01", truth.par)
    truth.par <- gsub("20", "02", truth.par)
    truth.par <- gsub("21", "12", truth.par)

    if(any(truth.par=='00')){
      index <- truth.par=='00'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(model.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(model.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2p00 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pt00 <- nrow(accurpartab2)
      accur.cn2p00 <- cn2p00 / cn2pt00

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1p00 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pt00 <- nrow(accurpartab2)
      accur.cn1p00 <- cn1p00 / cn1pt00
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0p00 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pt00 <- nrow(accurpartab2)
      accur.cn0p00 <- cn0p00 / cn0pt00
    } else {
      cn2p00 <- NA
      cn2pt00 <- NA
      accur.cn2p00 <- cn2p00 / cn2pt00
      cn1p00 <- NA
      cn1pt00 <- NA
      accur.cn1p00 <- cn1p00 / cn1pt00
      cn0p00 <- NA
      cn0pt00 <- NA
      accur.cn0p00 <- cn0p00 / cn0pt00
    }

    if(any(truth.par=='01')){
      index <- truth.par=='01'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(model.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(model.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2p01 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pt01 <- nrow(accurpartab2)
      accur.cn2p01 <- cn2p01 / cn2pt01

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1p01 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pt01 <- nrow(accurpartab2)
      accur.cn1p01 <- cn1p01 / cn1pt01
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0p01 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pt01 <- nrow(accurpartab2)
      accur.cn0p01 <- cn0p01 / cn0pt01
    } else {
      cn2p01 <- NA
      cn2pt01 <- NA
      accur.cn2p01 <- cn2p01 / cn2pt01
      cn1p01 <- NA
      cn1pt01 <- NA
      accur.cn1p01 <- cn1p01 / cn1pt01
      cn0p01 <- NA
      cn0pt01 <- NA
      accur.cn0p01 <- cn0p01 / cn0pt01
    }


    if(any(truth.par=='02')){
      index <- truth.par=='02'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(model.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(model.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2p02 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pt02 <- nrow(accurpartab2)
      accur.cn2p02 <- cn2p02 / cn2pt02

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1p02 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pt02 <- nrow(accurpartab2)
      accur.cn1p02 <- cn1p02 / cn1pt02
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0p02 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pt02 <- nrow(accurpartab2)
      accur.cn0p02 <- cn0p02 / cn0pt02
    } else {
      cn2p02 <- NA
      cn2pt02 <- NA
      accur.cn2p02 <- cn2p02 / cn2pt02
      cn1p02 <- NA
      cn1pt02 <- NA
      accur.cn1p02 <- cn1p02 / cn1pt02
      cn0p02 <- NA
      cn0pt02 <- NA
      accur.cn0p02 <- cn0p02 / cn0pt02
    }


    if(any(truth.par=='11')){
      index <- truth.par=='11'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(model.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(model.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2p11 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pt11 <- nrow(accurpartab2)
      accur.cn2p11 <- cn2p11 / cn2pt11

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1p11 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pt11 <- nrow(accurpartab2)
      accur.cn1p11 <- cn1p11 / cn1pt11
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0p11 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pt11 <- nrow(accurpartab2)
      accur.cn0p11 <- cn0p11 / cn0pt11
    } else {
      cn2p11 <- NA
      cn2pt11 <- NA
      accur.cn2p11 <- cn2p11 / cn2pt11
      cn1p11 <- NA
      cn1pt11 <- NA
      accur.cn1p11 <- cn1p11 / cn1pt11
      cn0p11 <- NA
      cn0pt11 <- NA
      accur.cn0p11 <- cn0p11 / cn0pt11
    }


    if(any(truth.par=='12')){
      index <- truth.par=='12'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(model.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(model.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2p12 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pt12 <- nrow(accurpartab2)
      accur.cn2p12 <- cn2p12 / cn2pt12

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1p12 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pt12 <- nrow(accurpartab2)
      accur.cn1p12 <- cn1p12 / cn1pt12
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0p12 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pt12 <- nrow(accurpartab2)
      accur.cn0p12 <- cn0p12 / cn0pt12
    } else {
      cn2p12 <- NA
      cn2pt12 <- NA
      accur.cn2p12 <- cn2p12 / cn2pt12
      cn1p12 <- NA
      cn1pt12 <- NA
      accur.cn1p12 <- cn1p12 / cn1pt12
      cn0p12 <- NA
      cn0pt12 <- NA
      accur.cn0p12 <- cn0p12 / cn0pt12
    }


    if(any(truth.par=='22')){
      index <- truth.par=='22'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(model.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(model.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2p22 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pt22 <- nrow(accurpartab2)
      accur.cn2p22 <- cn2p22 / cn2pt22

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1p22 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pt22 <- nrow(accurpartab2)
      accur.cn1p22 <- cn1p22 / cn1pt22
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0p22 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pt22 <- nrow(accurpartab2)
      accur.cn0p22 <- cn0p22 / cn0pt22
    } else {
      cn2p22 <- NA
      cn2pt22 <- NA
      accur.cn2p22 <- cn2p22 / cn2pt22
      cn1p22 <- NA
      cn1pt22 <- NA
      accur.cn1p22 <- cn1p22 / cn1pt22
      cn0p22 <- NA
      cn0pt22 <- NA
      accur.cn0p22 <- cn0p22 / cn0pt22
    }

    # MB parents
    truth.call.par <- params.in$TruthCNpar
    mb.call.par <- params.in$MBCNpar
    parentsep <- seq(1,length(truth.call.par), by=2)
    parentsep1 <- seq(2,length(truth.call.par), by=2)
    truth.call.par1 <- truth.call.par[parentsep]
    truth.call.par2 <- truth.call.par[parentsep1]
    mb.call.par1 <- mb.call.par[parentsep]
    mb.call.par2 <- mb.call.par[parentsep1]
    parcheck1 <- truth.call.par1 == mb.call.par1
    parcheck2 <- truth.call.par2 == mb.call.par2
    parallcorrect <- parcheck1 == parcheck2

    #accur.xtab <- table(truth.call.par, mb.call.par)
    #accur.default <- matrix(rep(0,9), nrow=3, ncol=3)
    #row.names(accur.default) <- c("0","1","2")
    #accur.default <- data.frame(accur.default)
    #accur.default <- rownames_to_column(accur.default)
    #accur.default$rowname <- as.numeric(accur.default$rowname)
    #indexcol <- as.numeric(row.names(accur.xtab))
    #accur.xtab2 <- cbind(indexcol, accur.xtab)
    #accur.xtab2 <- as.data.frame(accur.xtab2)
    #accur.xtab3 <- left_join(accur.default, accur.xtab2, by=c("rowname"="indexcol"))
    #cl <- ncol(accur.xtab3)
    #rn <- colnames(accur.xtab3)[5:cl]
    #accur.xtab31 <- data.frame(accur.xtab3[,5:cl])
    #colnames(accur.xtab31) <- rn

    #accur.default2 <- accur.default[,2:4]
    #colnames(accur.default2)[1:3] <- c(0,1,2)
    #accur.xtab4 <- smartbind(accur.default2, accur.xtab31)
    #accur.xtab4 <- accur.xtab4[-c(1:3),]
    # note this section is specific to deletion
    #accur.mb.cn2 <- round(accur.xtab4[3,3] / sum(accur.xtab4[3,], na.rm=T), 4)
    # accur.mb.cn1 <- round(accur.xtab4[2,2] / sum(accur.xtab4[2,], na.rm=T), 4)
    # accur.mb.cn0 <- round(accur.xtab4[1,1] / sum(accur.xtab4[1,], na.rm=T), 4)

    accurpartab <- data.frame(cbind(mb.call.par, truth.call.par))
    accurpartab2 <- accurpartab[accurpartab$truth.call.par==2,]
    cn2mbp <- sum(accurpartab2$mb.call.par==accurpartab2$truth.call.par)
    cn2mbpt <- nrow(accurpartab2)
    accur.mb.cn2 <- cn2mbp / cn2mbpt
    accurpartab2 <- accurpartab[accurpartab$truth.call.par==1,]
    cn1mbp <- sum(accurpartab2$mb.call.par==accurpartab2$truth.call.par)
    cn1mbpt <- nrow(accurpartab2)
    accur.mb.cn1 <- cn1mbp / cn1mbpt
    accurpartab2 <- accurpartab[accurpartab$truth.call.par==0,]
    cn0mbp <- sum(accurpartab2$mb.call.par==accurpartab2$truth.call.par)
    cn0mbpt <- nrow(accurpartab2)
    accur.mb.cn0 <- cn0mbp / cn0mbpt

    if(any(truth.par=='00')){
      index <- truth.par=='00'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(mb.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(mb.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbp00 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpt00 <- nrow(accurpartab2)
      accur.cn2mbp00 <- cn2mbp00 / cn2mbpt00

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbp00 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpt00 <- nrow(accurpartab2)
      accur.cn1mbp00 <- cn1mbp00 / cn1mbpt00
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbp00 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpt00 <- nrow(accurpartab2)
      accur.cn0mbp00 <- cn0mbp00 / cn0mbpt00
    } else {
      cn2mbp00 <- NA
      cn2mbpt00 <- NA
      accur.cn2mbp00 <- cn2mbp00 / cn2mbpt00
      cn1mbp00 <- NA
      cn1mbpt00 <- NA
      accur.cn1mbp00 <- cn1mbp00 / cn1mbpt00
      cn0mbp00 <- NA
      cn0mbpt00 <- NA
      accur.cn0mbp00 <- cn0mbp00 / cn0mbpt00
    }

    if(any(truth.par=='01')){
      index <- truth.par=='01'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(mb.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(mb.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbp01 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpt01 <- nrow(accurpartab2)
      accur.cn2mbp01 <- cn2mbp01 / cn2mbpt01

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbp01 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpt01 <- nrow(accurpartab2)
      accur.cn1mbp01 <- cn1mbp01 / cn1mbpt01
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbp01 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpt01 <- nrow(accurpartab2)
      accur.cn0mbp01 <- cn0mbp01 / cn0mbpt01
    } else {
      cn2mbp01 <- NA
      cn2mbpt01 <- NA
      accur.cn2mbp01 <- cn2mbp01 / cn2mbpt01
      cn1mbp01 <- NA
      cn1mbpt01 <- NA
      accur.cn1mbp01 <- cn1mbp01 / cn1mbpt01
      cn0mbp01 <- NA
      cn0mbpt01 <- NA
      accur.cn0mbp01 <- cn0mbp01 / cn0mbpt01
    }


    if(any(truth.par=='02')){
      index <- truth.par=='02'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(mb.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(mb.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbp02 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpt02 <- nrow(accurpartab2)
      accur.cn2mbp02 <- cn2mbp02 / cn2mbpt02

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbp02 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpt02 <- nrow(accurpartab2)
      accur.cn1mbp02 <- cn1mbp02 / cn1mbpt02
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbp02 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpt02 <- nrow(accurpartab2)
      accur.cn0mbp02 <- cn0mbp02 / cn0mbpt02
    } else {
      cn2mbp02 <- NA
      cn2mbpt02 <- NA
      accur.cn2mbp02 <- cn2mbp02 / cn2mbpt02
      cn1mbp02 <- NA
      cn1mbpt02 <- NA
      accur.cn1mbp02 <- cn1mbp02 / cn1mbpt02
      cn0mbp02 <- NA
      cn0mbpt02 <- NA
      accur.cn0mbp02 <- cn0mbp02 / cn0mbpt02
    }


    if(any(truth.par=='11')){
      index <- truth.par=='11'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(mb.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(mb.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbp11 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpt11 <- nrow(accurpartab2)
      accur.cn2mbp11 <- cn2mbp11 / cn2mbpt11

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbp11 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpt11 <- nrow(accurpartab2)
      accur.cn1mbp11 <- cn1mbp11 / cn1mbpt11
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbp11 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpt11 <- nrow(accurpartab2)
      accur.cn0mbp11 <- cn0mbp11 / cn0mbpt11
    } else {
      cn2mbp11 <- NA
      cn2mbpt11 <- NA
      accur.cn2mbp11 <- cn2mbp11 / cn2mbpt11
      cn1mbp11 <- NA
      cn1mbpt11 <- NA
      accur.cn1mbp11 <- cn1mbp11 / cn1mbpt11
      cn0mbp11 <- NA
      cn0mbpt11 <- NA
      accur.cn0mbp11 <- cn0mbp11 / cn0mbpt11
    }


    if(any(truth.par=='12')){
      index <- truth.par=='12'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(mb.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(mb.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbp12 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpt12 <- nrow(accurpartab2)
      accur.cn2mbp12 <- cn2mbp12 / cn2mbpt12

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbp12 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpt12 <- nrow(accurpartab2)
      accur.cn1mbp12 <- cn1mbp12 / cn1mbpt12
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbp12 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpt12 <- nrow(accurpartab2)
      accur.cn0mbp12 <- cn0mbp12 / cn0mbpt12
    } else {
      cn2mbp12 <- NA
      cn2mbpt12 <- NA
      accur.cn2mbp12 <- cn2mbp12 / cn2mbpt12
      cn1mbp12 <- NA
      cn1mbpt12 <- NA
      accur.cn1mbp12 <- cn1mbp12 / cn1mbpt12
      cn0mbp12 <- NA
      cn0mbpt12 <- NA
      accur.cn0mbp12 <- cn0mbp12 / cn0mbpt12
    }


    if(any(truth.par=='22')){
      index <- truth.par=='22'
      truth.rel1 <- data.frame(truth.call.par1[index])
      names(truth.rel1) <- "truth.rel"
      truth.rel2 <- data.frame(truth.call.par2[index])
      names(truth.rel2) <- "truth.rel"
      model.rel1 <- data.frame(mb.call.par1[index])
      names(model.rel1) <- "model.rel"
      model.rel2 <- data.frame(mb.call.par2[index])
      names(model.rel2) <- "model.rel"
      truth.rel <- rbind(truth.rel1, truth.rel2)
      model.rel <- rbind(model.rel1, model.rel2)

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbp22 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpt22 <- nrow(accurpartab2)
      accur.cn2mbp22 <- cn2mbp22 / cn2mbpt22

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbp22 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpt22 <- nrow(accurpartab2)
      accur.cn1mbp22 <- cn1mbp22 / cn1mbpt22
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbp22 <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpt22 <- nrow(accurpartab2)
      accur.cn0mbp22 <- cn0mbp22 / cn0mbpt22
    } else {
      cn2mbp22 <- NA
      cn2mbpt22 <- NA
      accur.cn2mbp22 <- cn2mbp22 / cn2mbpt22
      cn1mbp22 <- NA
      cn1mbpt22 <- NA
      accur.cn1mbp22 <- cn1mbp22 / cn1mbpt22
      cn0mbp22 <- NA
      cn0mbpt22 <- NA
      accur.cn0mbp22 <- cn0mbp22 / cn0mbpt22
    }


    #Trio kids
    truth.call.off <- params.in$TruthCNoff
    model.call.off <- params.in$TrioCNoff
    accur.xtab <- table(truth.call.off, model.call.off)
    truth.call.off.par.correct <- truth.call.off[parallcorrect]
    model.call.off.par.correct <- model.call.off[parallcorrect]
    truth.call.off.par.wrong <- truth.call.off[!parallcorrect]
    model.call.off.par.wrong <- model.call.off[!parallcorrect]
    accur.parcor <- (sum(model.call.off.par.correct==truth.call.off.par.correct)) / length(truth.call.off.par.correct)
    accur.parincor <- (sum(model.call.off.par.wrong==truth.call.off.par.wrong)) / length(truth.call.off.par.wrong)

    accur.default <- matrix(rep(0,9), nrow=3, ncol=3)
    row.names(accur.default) <- c("0","1","2")
    accur.default <- data.frame(accur.default)
    accur.default <- rownames_to_column(accur.default)
    accur.default$rowname <- as.numeric(accur.default$rowname)
    indexcol <- as.numeric(row.names(accur.xtab))
    accur.xtab2 <- cbind(indexcol, accur.xtab)
    accur.xtab2 <- as.data.frame(accur.xtab2)
    accur.xtab3 <- left_join(accur.default, accur.xtab2, by=c("rowname"="indexcol"))
    cl <- ncol(accur.xtab3)
    rn <- colnames(accur.xtab3)[5:cl]
    accur.xtab31 <- data.frame(accur.xtab3[,5:cl])
    colnames(accur.xtab31) <- rn
    #accur.xtab3[is.na(accur.xtab3)] <- 0

    accur.default2 <- accur.default[,2:4]
    colnames(accur.default2)[1:3] <- c(0,1,2)
    accur.xtab4 <- smartbind(accur.default2, accur.xtab31)
    accur.xtab4 <- accur.xtab4[-c(1:3),]
    # note this section is specific to deletion
    accur.cn2.off <- round(accur.xtab4[3,3] / sum(accur.xtab4[3,], na.rm=T), 4)
    accur.cn1.off <- round(accur.xtab4[2,2] / sum(accur.xtab4[2,], na.rm=T), 4)
    accur.cn0.off <- round(accur.xtab4[1,1] / sum(accur.xtab4[1,], na.rm=T), 4)

    accurpartab <- data.frame(cbind(model.call.off, truth.call.off))
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==2,]
    cn2po <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn2pot <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==1,]
    cn1po <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn1pot <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==0,]
    cn0po <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn0pot <- nrow(accurpartab2)

    # create accuracy counts by parental genotype
    truth.par <- paste0(truth.call.par1, truth.call.par2)
    truth.par <- gsub("10", "01", truth.par)
    truth.par <- gsub("20", "02", truth.par)
    truth.par <- gsub("21", "12", truth.par)

    if(any(truth.par=='00')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='00'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po00pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot00pc <- nrow(accurpartab2)
      accur.cn2po00pc <- cn2po00pc / cn2pot00pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po00pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot00pc <- nrow(accurpartab2)
      accur.cn1po00pc <- cn1po00pc / cn1pot00pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po00pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot00pc <- nrow(accurpartab2)
      accur.cn0po00pc <- cn0po00pc / cn0pot00pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='00'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po00pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot00pno <- nrow(accurpartab2)
      accur.cn2po00pno <- cn2po00pno / cn2pot00pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po00pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot00pno <- nrow(accurpartab2)
      accur.cn1po00pno <- cn1po00pno / cn1pot00pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po00pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot00pno <- nrow(accurpartab2)
      accur.cn0po00pno <- cn0po00pno / cn0pot00pno
    } else {
      cn2po00pc <- NA
      cn2pot00pc <- NA
      accur.cn2po00pc <- cn2po00pc / cn2pot00pc
      cn1po00pc <- NA
      cn1pot00pc <- NA
      accur.cn1po00pc <- cn1po00pc / cn1pot00pc
      cn0po00pc <- NA
      cn0pot00pc <- NA
      accur.cn0po00pc <- cn0po00pc / cn0pot00pc
      cn2po00pno <- NA
      cn2pot00pno <- NA
      accur.cn2po00pno <- cn2po00pno / cn2pot00pno
      cn1po00pno <- NA
      cn1pot00pno <- NA
      accur.cn1po00pno <- cn1po00pno / cn1pot00pno
      cn0po00pno <- NA
      cn0pot00pno <- NA
      accur.cn0po00pno <- cn0po00pno / cn0pot00pno
    }

    if(any(truth.par=='01')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='01'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po01pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot01pc <- nrow(accurpartab2)
      accur.cn2po01pc <- cn2po01pc / cn2pot01pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po01pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot01pc <- nrow(accurpartab2)
      accur.cn1po01pc <- cn1po01pc / cn1pot01pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po01pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot01pc <- nrow(accurpartab2)
      accur.cn0po01pc <- cn0po01pc / cn0pot01pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='01'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po01pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot01pno <- nrow(accurpartab2)
      accur.cn2po01pno <- cn2po01pno / cn2pot01pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po01pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot01pno <- nrow(accurpartab2)
      accur.cn1po01pno <- cn1po01pno / cn1pot01pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po01pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot01pno <- nrow(accurpartab2)
      accur.cn0po01pno <- cn0po01pno / cn0pot01pno
    } else {
      cn2po01pc <- NA
      cn2pot01pc <- NA
      accur.cn2po01pc <- cn2po01pc / cn2pot01pc
      cn1po01pc <- NA
      cn1pot01pc <- NA
      accur.cn1po01pc <- cn1po01pc / cn1pot01pc
      cn0po01pc <- NA
      cn0pot01pc <- NA
      accur.cn0po01pc <- cn0po01pc / cn0pot01pc
      cn2po01pno <- NA
      cn2pot01pno <- NA
      accur.cn2po01pno <- cn2po01pno / cn2pot01pno
      cn1po01pno <- NA
      cn1pot01pno <- NA
      accur.cn1po01pno <- cn1po01pno / cn1pot01pno
      cn0po01pno <- NA
      cn0pot01pno <- NA
      accur.cn0po01pno <- cn0po01pno / cn0pot01pno
    }

    if(any(truth.par=='02')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='02'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po02pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot02pc <- nrow(accurpartab2)
      accur.cn2po02pc <- cn2po02pc / cn2pot02pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po02pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot02pc <- nrow(accurpartab2)
      accur.cn1po02pc <- cn1po02pc / cn1pot02pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po02pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot02pc <- nrow(accurpartab2)
      accur.cn0po02pc <- cn0po02pc / cn0pot02pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='02'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po02pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot02pno <- nrow(accurpartab2)
      accur.cn2po02pno <- cn2po02pno / cn2pot02pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po02pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot02pno <- nrow(accurpartab2)
      accur.cn1po02pno <- cn1po02pno / cn1pot02pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po02pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot02pno <- nrow(accurpartab2)
      accur.cn0po02pno <- cn0po02pno / cn0pot02pno
    } else {
      cn2po02pc <- NA
      cn2pot02pc <- NA
      accur.cn2po02pc <- cn2po02pc / cn2pot02pc
      cn1po02pc <- NA
      cn1pot02pc <- NA
      accur.cn1po02pc <- cn1po02pc / cn1pot02pc
      cn0po02pc <- NA
      cn0pot02pc <- NA
      accur.cn0po02pc <- cn0po02pc / cn0pot02pc
      cn2po02pno <- NA
      cn2pot02pno <- NA
      accur.cn2po02pno <- cn2po02pno / cn2pot02pno
      cn1po02pno <- NA
      cn1pot02pno <- NA
      accur.cn1po02pno <- cn1po02pno / cn1pot02pno
      cn0po02pno <- NA
      cn0pot02pno <- NA
      accur.cn0po02pno <- cn0po02pno / cn0pot02pno
    }

    if(any(truth.par=='11')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='11'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po11pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot11pc <- nrow(accurpartab2)
      accur.cn2po11pc <- cn2po11pc / cn2pot11pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po11pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot11pc <- nrow(accurpartab2)
      accur.cn1po11pc <- cn1po11pc / cn1pot11pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po11pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot11pc <- nrow(accurpartab2)
      accur.cn0po11pc <- cn0po11pc / cn0pot11pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='11'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po11pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot11pno <- nrow(accurpartab2)
      accur.cn2po11pno <- cn2po11pno / cn2pot11pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po11pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot11pno <- nrow(accurpartab2)
      accur.cn1po11pno <- cn1po11pno / cn1pot11pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po11pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot11pno <- nrow(accurpartab2)
      accur.cn0po11pno <- cn0po11pno / cn0pot11pno
    } else {
      cn2po11pc <- NA
      cn2pot11pc <- NA
      accur.cn2po11pc <- cn2po11pc / cn2pot11pc
      cn1po11pc <- NA
      cn1pot11pc <- NA
      accur.cn1po11pc <- cn1po11pc / cn1pot11pc
      cn0po11pc <- NA
      cn0pot11pc <- NA
      accur.cn0po11pc <- cn0po11pc / cn0pot11pc
      cn2po11pno <- NA
      cn2pot11pno <- NA
      accur.cn2po11pno <- cn2po11pno / cn2pot11pno
      cn1po11pno <- NA
      cn1pot11pno <- NA
      accur.cn1po11pno <- cn1po11pno / cn1pot11pno
      cn0po11pno <- NA
      cn0pot11pno <- NA
      accur.cn0po11pno <- cn0po11pno / cn0pot11pno
    }

    if(any(truth.par=='12')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='12'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po12pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot12pc <- nrow(accurpartab2)
      accur.cn2po12pc <- cn2po12pc / cn2pot12pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po12pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot12pc <- nrow(accurpartab2)
      accur.cn1po12pc <- cn1po12pc / cn1pot12pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po12pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot12pc <- nrow(accurpartab2)
      accur.cn0po12pc <- cn0po12pc / cn0pot12pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='12'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po12pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot12pno <- nrow(accurpartab2)
      accur.cn2po12pno <- cn2po12pno / cn2pot12pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po12pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot12pno <- nrow(accurpartab2)
      accur.cn1po12pno <- cn1po12pno / cn1pot12pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po12pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot12pno <- nrow(accurpartab2)
      accur.cn0po12pno <- cn0po12pno / cn0pot12pno
    } else {
      cn2po12pc <- NA
      cn2pot12pc <- NA
      accur.cn2po12pc <- cn2po12pc / cn2pot12pc
      cn1po12pc <- NA
      cn1pot12pc <- NA
      accur.cn1po12pc <- cn1po12pc / cn1pot12pc
      cn0po12pc <- NA
      cn0pot12pc <- NA
      accur.cn0po12pc <- cn0po12pc / cn0pot12pc
      cn2po12pno <- NA
      cn2pot12pno <- NA
      accur.cn2po12pno <- cn2po12pno / cn2pot12pno
      cn1po12pno <- NA
      cn1pot12pno <- NA
      accur.cn1po12pno <- cn1po12pno / cn1pot12pno
      cn0po12pno <- NA
      cn0pot12pno <- NA
      accur.cn0po12pno <- cn0po12pno / cn0pot12pno
    }

    if(any(truth.par=='22')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='22'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po22pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot22pc <- nrow(accurpartab2)
      accur.cn2po22pc <- cn2po22pc / cn2pot22pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po22pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot22pc <- nrow(accurpartab2)
      accur.cn1po22pc <- cn1po22pc / cn1pot22pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po22pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot22pc <- nrow(accurpartab2)
      accur.cn0po22pc <- cn0po22pc / cn0pot22pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='22'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(model.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2po22pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2pot22pno <- nrow(accurpartab2)
      accur.cn2po22pno <- cn2po22pno / cn2pot22pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1po22pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1pot22pno <- nrow(accurpartab2)
      accur.cn1po22pno <- cn1po22pno / cn1pot22pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0po22pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0pot22pno <- nrow(accurpartab2)
      accur.cn0po22pno <- cn0po22pno / cn0pot22pno
    } else {
      cn2po22pc <- NA
      cn2pot22pc <- NA
      accur.cn2po22pc <- cn2po22pc / cn2pot22pc
      cn1po22pc <- NA
      cn1pot22pc <- NA
      accur.cn1po22pc <- cn1po22pc / cn1pot22pc
      cn0po22pc <- NA
      cn0pot22pc <- NA
      accur.cn0po22pc <- cn0po22pc / cn0pot22pc
      cn2po22pno <- NA
      cn2pot22pno <- NA
      accur.cn2po22pno <- cn2po22pno / cn2pot22pno
      cn1po22pno <- NA
      cn1pot22pno <- NA
      accur.cn1po22pno <- cn1po22pno / cn1pot22pno
      cn0po22pno <- NA
      cn0pot22pno <- NA
      accur.cn0po22pno <- cn0po22pno / cn0pot22pno
    }

    # MB child
    truth.call.off <- params.in$TruthCNoff
    mb.call.off <- params.in$MBCNoff
    accur.xtab <- table(truth.call.off, mb.call.off)

    accur.default <- matrix(rep(0,9), nrow=3, ncol=3)
    row.names(accur.default) <- c("0","1","2")
    accur.default <- data.frame(accur.default)
    accur.default <- rownames_to_column(accur.default)
    accur.default$rowname <- as.numeric(accur.default$rowname)
    indexcol <- as.numeric(row.names(accur.xtab))
    accur.xtab2 <- cbind(indexcol, accur.xtab)
    accur.xtab2 <- as.data.frame(accur.xtab2)
    accur.xtab3 <- left_join(accur.default, accur.xtab2, by=c("rowname"="indexcol"))
    cl <- ncol(accur.xtab3)
    rn <- colnames(accur.xtab3)[5:cl]
    accur.xtab31 <- data.frame(accur.xtab3[,5:cl])
    colnames(accur.xtab31) <- rn
    #accur.xtab3[is.na(accur.xtab3)] <- 0

    accur.default2 <- accur.default[,2:4]
    colnames(accur.default2)[1:3] <- c(0,1,2)
    accur.xtab4 <- smartbind(accur.default2, accur.xtab31)
    accur.xtab4 <- accur.xtab4[-c(1:3),]
    # note this section is specific to deletion
    accur.mb.cn2.off <- round(accur.xtab4[3,3] / sum(accur.xtab4[3,], na.rm=T), 4)
    #if((2 %in% params.in$MBCNcall) == FALSE){
    #accur.mb.cn2.off[is.na(accur.mb.cn2.off)] <- 0
    #}
    accur.mb.cn1.off <- round(accur.xtab4[2,2] / sum(accur.xtab4[2,], na.rm=T), 4)
    #if((1 %in% params.in$MBCNcall) == FALSE){
    #accur.mb.cn1.off[is.na(accur.mb.cn1.off)] <- 0
    #}
    accur.mb.cn0.off <- round(accur.xtab4[1,1] / sum(accur.xtab4[1,], na.rm=T), 4)
    #if((0 %in% params.in$MBCNcall) == FALSE){
    #accur.mb.cn0.off[is.na(accur.mb.cn0.off)] <- 0
    #}

    accurpartab <- data.frame(cbind(mb.call.off, truth.call.off))
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==2,]
    cn2mbpo <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn2mbpot <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==1,]
    cn1mbpo <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn1mbpot <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==0,]
    cn0mbpo <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn0mbpot <- nrow(accurpartab2)

    truth.par <- paste0(truth.call.par1, truth.call.par2)
    truth.par <- gsub("10", "01", truth.par)
    truth.par <- gsub("20", "02", truth.par)
    truth.par <- gsub("21", "12", truth.par)

    if(any(truth.par=='00')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='00'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo00pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot00pc <- nrow(accurpartab2)
      accur.cn2mbpo00pc <- cn2mbpo00pc / cn2mbpot00pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo00pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot00pc <- nrow(accurpartab2)
      accur.cn1mbpo00pc <- cn1mbpo00pc / cn1mbpot00pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo00pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot00pc <- nrow(accurpartab2)
      accur.cn0mbpo00pc <- cn0mbpo00pc / cn0mbpot00pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='00'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo00pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot00pno <- nrow(accurpartab2)
      accur.cn2mbpo00pno <- cn2mbpo00pno / cn2mbpot00pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo00pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot00pno <- nrow(accurpartab2)
      accur.cn1mbpo00pno <- cn1mbpo00pno / cn1mbpot00pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo00pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot00pno <- nrow(accurpartab2)
      accur.cn0mbpo00pno <- cn0mbpo00pno / cn0mbpot00pno
    } else {
      cn2mbpo00pc <- NA
      cn2mbpot00pc <- NA
      accur.cn2mbpo00pc <- cn2mbpo00pc / cn2mbpot00pc
      cn1mbpo00pc <- NA
      cn1mbpot00pc <- NA
      accur.cn1mbpo00pc <- cn1mbpo00pc / cn1mbpot00pc
      cn0mbpo00pc <- NA
      cn0mbpot00pc <- NA
      accur.cn0mbpo00pc <- cn0mbpo00pc / cn0mbpot00pc
      cn2mbpo00pno <- NA
      cn2mbpot00pno <- NA
      accur.cn2mbpo00pno <- cn2mbpo00pno / cn2mbpot00pno
      cn1mbpo00pno <- NA
      cn1mbpot00pno <- NA
      accur.cn1mbpo00pno <- cn1mbpo00pno / cn1mbpot00pno
      cn0mbpo00pno <- NA
      cn0mbpot00pno <- NA
      accur.cn0mbpo00pno <- cn0mbpo00pno / cn0mbpot00pno
    }

    if(any(truth.par=='01')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='01'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo01pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot01pc <- nrow(accurpartab2)
      accur.cn2mbpo01pc <- cn2mbpo01pc / cn2mbpot01pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo01pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot01pc <- nrow(accurpartab2)
      accur.cn1mbpo01pc <- cn1mbpo01pc / cn1mbpot01pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo01pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot01pc <- nrow(accurpartab2)
      accur.cn0mbpo01pc <- cn0mbpo01pc / cn0mbpot01pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='01'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo01pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot01pno <- nrow(accurpartab2)
      accur.cn2mbpo01pno <- cn2mbpo01pno / cn2mbpot01pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo01pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot01pno <- nrow(accurpartab2)
      accur.cn1mbpo01pno <- cn1mbpo01pno / cn1mbpot01pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo01pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot01pno <- nrow(accurpartab2)
      accur.cn0mbpo01pno <- cn0mbpo01pno / cn0mbpot01pno
    } else {
      cn2mbpo01pc <- NA
      cn2mbpot01pc <- NA
      accur.cn2mbpo01pc <- cn2mbpo01pc / cn2mbpot01pc
      cn1mbpo01pc <- NA
      cn1mbpot01pc <- NA
      accur.cn1mbpo01pc <- cn1mbpo01pc / cn1mbpot01pc
      cn0mbpo01pc <- NA
      cn0mbpot01pc <- NA
      accur.cn0mbpo01pc <- cn0mbpo01pc / cn0mbpot01pc
      cn2mbpo01pno <- NA
      cn2mbpot01pno <- NA
      accur.cn2mbpo01pno <- cn2mbpo01pno / cn2mbpot01pno
      cn1mbpo01pno <- NA
      cn1mbpot01pno <- NA
      accur.cn1mbpo01pno <- cn1mbpo01pno / cn1mbpot01pno
      cn0mbpo01pno <- NA
      cn0mbpot01pno <- NA
      accur.cn0mbpo01pno <- cn0mbpo01pno / cn0mbpot01pno
    }

    if(any(truth.par=='02')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='02'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo02pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot02pc <- nrow(accurpartab2)
      accur.cn2mbpo02pc <- cn2mbpo02pc / cn2mbpot02pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo02pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot02pc <- nrow(accurpartab2)
      accur.cn1mbpo02pc <- cn1mbpo02pc / cn1mbpot02pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo02pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot02pc <- nrow(accurpartab2)
      accur.cn0mbpo02pc <- cn0mbpo02pc / cn0mbpot02pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='02'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo02pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot02pno <- nrow(accurpartab2)
      accur.cn2mbpo02pno <- cn2mbpo02pno / cn2mbpot02pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo02pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot02pno <- nrow(accurpartab2)
      accur.cn1mbpo02pno <- cn1mbpo02pno / cn1mbpot02pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo02pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot02pno <- nrow(accurpartab2)
      accur.cn0mbpo02pno <- cn0mbpo02pno / cn0mbpot02pno
    } else {
      cn2mbpo02pc <- NA
      cn2mbpot02pc <- NA
      accur.cn2mbpo02pc <- cn2mbpo02pc / cn2mbpot02pc
      cn1mbpo02pc <- NA
      cn1mbpot02pc <- NA
      accur.cn1mbpo02pc <- cn1mbpo02pc / cn1mbpot02pc
      cn0mbpo02pc <- NA
      cn0mbpot02pc <- NA
      accur.cn0mbpo02pc <- cn0mbpo02pc / cn0mbpot02pc
      cn2mbpo02pno <- NA
      cn2mbpot02pno <- NA
      accur.cn2mbpo02pno <- cn2mbpo02pno / cn2mbpot02pno
      cn1mbpo02pno <- NA
      cn1mbpot02pno <- NA
      accur.cn1mbpo02pno <- cn1mbpo02pno / cn1mbpot02pno
      cn0mbpo02pno <- NA
      cn0mbpot02pno <- NA
      accur.cn0mbpo02pno <- cn0mbpo02pno / cn0mbpot02pno
    }

    if(any(truth.par=='11')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='11'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo11pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot11pc <- nrow(accurpartab2)
      accur.cn2mbpo11pc <- cn2mbpo11pc / cn2mbpot11pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo11pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot11pc <- nrow(accurpartab2)
      accur.cn1mbpo11pc <- cn1mbpo11pc / cn1mbpot11pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo11pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot11pc <- nrow(accurpartab2)
      accur.cn0mbpo11pc <- cn0mbpo11pc / cn0mbpot11pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='11'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo11pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot11pno <- nrow(accurpartab2)
      accur.cn2mbpo11pno <- cn2mbpo11pno / cn2mbpot11pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo11pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot11pno <- nrow(accurpartab2)
      accur.cn1mbpo11pno <- cn1mbpo11pno / cn1mbpot11pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo11pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot11pno <- nrow(accurpartab2)
      accur.cn0mbpo11pno <- cn0mbpo11pno / cn0mbpot11pno
    } else {
      cn2mbpo11pc <- NA
      cn2mbpot11pc <- NA
      accur.cn2mbpo11pc <- cn2mbpo11pc / cn2mbpot11pc
      cn1mbpo11pc <- NA
      cn1mbpot11pc <- NA
      accur.cn1mbpo11pc <- cn1mbpo11pc / cn1mbpot11pc
      cn0mbpo11pc <- NA
      cn0mbpot11pc <- NA
      accur.cn0mbpo11pc <- cn0mbpo11pc / cn0mbpot11pc
      cn2mbpo11pno <- NA
      cn2mbpot11pno <- NA
      accur.cn2mbpo11pno <- cn2mbpo11pno / cn2mbpot11pno
      cn1mbpo11pno <- NA
      cn1mbpot11pno <- NA
      accur.cn1mbpo11pno <- cn1mbpo11pno / cn1mbpot11pno
      cn0mbpo11pno <- NA
      cn0mbpot11pno <- NA
      accur.cn0mbpo11pno <- cn0mbpo11pno / cn0mbpot11pno
    }

    if(any(truth.par=='12')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='12'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo12pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot12pc <- nrow(accurpartab2)
      accur.cn2mbpo12pc <- cn2mbpo12pc / cn2mbpot12pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo12pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot12pc <- nrow(accurpartab2)
      accur.cn1mbpo12pc <- cn1mbpo12pc / cn1mbpot12pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo12pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot12pc <- nrow(accurpartab2)
      accur.cn0mbpo12pc <- cn0mbpo12pc / cn0mbpot12pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='12'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo12pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot12pno <- nrow(accurpartab2)
      accur.cn2mbpo12pno <- cn2mbpo12pno / cn2mbpot12pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo12pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot12pno <- nrow(accurpartab2)
      accur.cn1mbpo12pno <- cn1mbpo12pno / cn1mbpot12pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo12pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot12pno <- nrow(accurpartab2)
      accur.cn0mbpo12pno <- cn0mbpo12pno / cn0mbpot12pno
    } else {
      cn2mbpo12pc <- NA
      cn2mbpot12pc <- NA
      accur.cn2mbpo12pc <- cn2mbpo12pc / cn2mbpot12pc
      cn1mbpo12pc <- NA
      cn1mbpot12pc <- NA
      accur.cn1mbpo12pc <- cn1mbpo12pc / cn1mbpot12pc
      cn0mbpo12pc <- NA
      cn0mbpot12pc <- NA
      accur.cn0mbpo12pc <- cn0mbpo12pc / cn0mbpot12pc
      cn2mbpo12pno <- NA
      cn2mbpot12pno <- NA
      accur.cn2mbpo12pno <- cn2mbpo12pno / cn2mbpot12pno
      cn1mbpo12pno <- NA
      cn1mbpot12pno <- NA
      accur.cn1mbpo12pno <- cn1mbpo12pno / cn1mbpot12pno
      cn0mbpo12pno <- NA
      cn0mbpot12pno <- NA
      accur.cn0mbpo12pno <- cn0mbpo12pno / cn0mbpot12pno
    }

    if(any(truth.par=='22')){
      truth.par.pc <- truth.par[parallcorrect]
      index <- truth.par.pc=='22'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo22pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot22pc <- nrow(accurpartab2)
      accur.cn2mbpo22pc <- cn2mbpo22pc / cn2mbpot22pc

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo22pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot22pc <- nrow(accurpartab2)
      accur.cn1mbpo22pc <- cn1mbpo22pc / cn1mbpot22pc
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo22pc <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot22pc <- nrow(accurpartab2)
      accur.cn0mbpo22pc <- cn0mbpo22pc / cn0mbpot22pc

      truth.par.pno <- truth.par[!parallcorrect]
      index <- truth.par.pno=='22'
      truth.rel <- data.frame(truth.call.off[index])
      names(truth.rel) <- "truth.rel"
      model.rel <- data.frame(mb.call.off[index])
      names(model.rel) <- "model.rel"

      accurpartab <- data.frame(cbind(model.rel, truth.rel))
      accurpartab2 <- accurpartab[accurpartab$truth.rel==2,]
      cn2mbpo22pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn2mbpot22pno <- nrow(accurpartab2)
      accur.cn2mbpo22pno <- cn2mbpo22pno / cn2mbpot22pno

      accurpartab2 <- accurpartab[accurpartab$truth.rel==1,]
      cn1mbpo22pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn1mbpot22pno <- nrow(accurpartab2)
      accur.cn1mbpo22pno <- cn1mbpo22pno / cn1mbpot22pno
      accurpartab2 <- accurpartab[accurpartab$truth.rel==0,]
      cn0mbpo22pno <- sum(accurpartab2$model.rel==accurpartab2$truth.rel)
      cn0mbpot22pno <- nrow(accurpartab2)
      accur.cn0mbpo22pno <- cn0mbpo22pno / cn0mbpot22pno
    } else {
      cn2mbpo22pc <- NA
      cn2mbpot22pc <- NA
      accur.cn2mbpo22pc <- cn2mbpo22pc / cn2mbpot22pc
      cn1mbpo22pc <- NA
      cn1mbpot22pc <- NA
      accur.cn1mbpo22pc <- cn1mbpo22pc / cn1mbpot22pc
      cn0mbpo22pc <- NA
      cn0mbpot22pc <- NA
      accur.cn0mbpo22pc <- cn0mbpo22pc / cn0mbpot22pc
      cn2mbpo22pno <- NA
      cn2mbpot22pno <- NA
      accur.cn2mbpo22pno <- cn2mbpo22pno / cn2mbpot22pno
      cn1mbpo22pno <- NA
      cn1mbpot22pno <- NA
      accur.cn1mbpo22pno <- cn1mbpo22pno / cn1mbpot22pno
      cn0mbpo22pno <- NA
      cn0mbpot22pno <- NA
      accur.cn0mbpo22pno <- cn0mbpo22pno / cn0mbpot22pno
    }

    #Trio Child Parents Correct
    truth.call.off <- params.in$TruthCNoff[parallcorrect]
    model.call.off <- params.in$TrioCNoff[parallcorrect]
    accur.xtab <- table(truth.call.off, model.call.off)
    accurpartab <- data.frame(cbind(model.call.off, truth.call.off))
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==2,]
    cn2popc <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn2potpc <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==1,]
    cn1popc <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn1potpc <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==0,]
    cn0popc <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn0potpc <- nrow(accurpartab2)

    # MB child Parents Correct
    truth.call.off <- params.in$TruthCNoff[parallcorrect]
    mb.call.off <- params.in$MBCNoff[parallcorrect]
    accur.xtab <- table(truth.call.off, mb.call.off)

    accurpartab <- data.frame(cbind(mb.call.off, truth.call.off))
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==2,]
    cn2mbpopc <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn2mbpotpc <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==1,]
    cn1mbpopc <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn1mbpotpc <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==0,]
    cn0mbpopc <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn0mbpotpc <- nrow(accurpartab2)

    #Trio Child Parents Wrong
    truth.call.off <- params.in$TruthCNoff[!parallcorrect]
    model.call.off <- params.in$TrioCNoff[!parallcorrect]
    accur.xtab <- table(truth.call.off, model.call.off)
    accurpartab <- data.frame(cbind(model.call.off, truth.call.off))
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==2,]
    cn2popw <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn2potpw <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==1,]
    cn1popw <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn1potpw <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==0,]
    cn0popw <- sum(accurpartab2$model.call.off==accurpartab2$truth.call.off)
    cn0potpw <- nrow(accurpartab2)

    # MB child Parents Wrong
    truth.call.off <- params.in$TruthCNoff[!parallcorrect]
    mb.call.off <- params.in$MBCNoff[!parallcorrect]
    accur.xtab <- table(truth.call.off, mb.call.off)

    accurpartab <- data.frame(cbind(mb.call.off, truth.call.off))
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==2,]
    cn2mbpopw <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn2mbpotpw <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==1,]
    cn1mbpopw <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn1mbpotpw <- nrow(accurpartab2)
    accurpartab2 <- accurpartab[accurpartab$truth.call.off==0,]
    cn0mbpopw <- sum(accurpartab2$mb.call.off==accurpartab2$truth.call.off)
    cn0mbpotpw <- nrow(accurpartab2)

    accuracy.table[i,1] <- p
    accuracy.table[i,2] <- N
    accuracy.table[i,3] <- delta
    accuracy.table[i,4] <- accur
    accuracy.table[i,5] <- accur.par
    accuracy.table[i,6] <- accur.child
    accuracy.table[i,7] <- accur.mb.par
    accuracy.table[i,8] <- accur.mb.indep
    accuracy.table[i,9] <- accur.cn2
    accuracy.table[i,10] <- accur.cn1
    accuracy.table[i,11] <- accur.cn0
    accuracy.table[i,12] <- accur.mb.cn2
    accuracy.table[i,13] <- accur.mb.cn1
    accuracy.table[i,14] <- accur.mb.cn0
    accuracy.table[i,15] <- accur.cn2.off
    accuracy.table[i,16] <- accur.cn1.off
    accuracy.table[i,17] <- accur.cn0.off
    accuracy.table[i,18] <- accur.mb.cn2.off
    accuracy.table[i,19] <- accur.mb.cn1.off
    accuracy.table[i,20] <- accur.mb.cn0.off
    accuracy.table[i,21] <- accur.parcor
    accuracy.table[i,22] <- accur.parincor
    accuracy.table[i,23] <- cn2p
    accuracy.table[i,24] <- cn2pt
    accuracy.table[i,25] <- cn1p
    accuracy.table[i,26] <- cn1pt
    accuracy.table[i,27] <- cn0p
    accuracy.table[i,28] <- cn0pt
    accuracy.table[i,29] <- cn2mbp
    accuracy.table[i,30] <- cn2mbpt
    accuracy.table[i,31] <- cn1mbp
    accuracy.table[i,32] <- cn1mbpt
    accuracy.table[i,33] <- cn0mbp
    accuracy.table[i,34] <- cn0mbpt
    accuracy.table[i,35] <- cn2po
    accuracy.table[i,36] <- cn2pot
    accuracy.table[i,37] <- cn1po
    accuracy.table[i,38] <- cn1pot
    accuracy.table[i,39] <- cn0po
    accuracy.table[i,40] <- cn0pot
    accuracy.table[i,41] <- cn2mbpo
    accuracy.table[i,42] <- cn2mbpot
    accuracy.table[i,43] <- cn1mbpo
    accuracy.table[i,44] <- cn1mbpot
    accuracy.table[i,45] <- cn0mbpo
    accuracy.table[i,46] <- cn0mbpot
    accuracy.table[i,47] <- cn2popc
    accuracy.table[i,48] <- cn2potpc
    accuracy.table[i,49] <- cn1popc
    accuracy.table[i,50] <- cn1potpc
    accuracy.table[i,51] <- cn0popc
    accuracy.table[i,52] <- cn0potpc
    accuracy.table[i,53] <- cn2mbpopc
    accuracy.table[i,54] <- cn2mbpotpc
    accuracy.table[i,55] <- cn1mbpopc
    accuracy.table[i,56] <- cn1mbpotpc
    accuracy.table[i,57] <- cn0mbpopc
    accuracy.table[i,58] <- cn0mbpotpc
    accuracy.table[i,59] <- cn2popw
    accuracy.table[i,60] <- cn2potpw
    accuracy.table[i,61] <- cn1popw
    accuracy.table[i,62] <- cn1potpw
    accuracy.table[i,63] <- cn0popw
    accuracy.table[i,64] <- cn0potpw
    accuracy.table[i,65] <- cn2mbpopw
    accuracy.table[i,66] <- cn2mbpotpw
    accuracy.table[i,67] <- cn1mbpopw
    accuracy.table[i,68] <- cn1mbpotpw
    accuracy.table[i,69] <- cn0mbpopw
    accuracy.table[i,70] <- cn0mbpotpw
    accuracy.table[i,71] <- cn2po00pc
    accuracy.table[i,72] <- cn2pot00pc
    accuracy.table[i,73] <- accur.cn2po00pc
    accuracy.table[i,74] <- cn1po00pc
    accuracy.table[i,75] <- cn1pot00pc
    accuracy.table[i,76] <- accur.cn1po00pc
    accuracy.table[i,77] <- cn0po00pc
    accuracy.table[i,78] <- cn0pot00pc
    accuracy.table[i,79] <- accur.cn0po00pc
    accuracy.table[i,80] <- cn2po01pc
    accuracy.table[i,81] <- cn2pot01pc
    accuracy.table[i,82] <- accur.cn2po01pc
    accuracy.table[i,83] <- cn1po01pc
    accuracy.table[i,84] <- cn1pot01pc
    accuracy.table[i,85] <- accur.cn1po01pc
    accuracy.table[i,86] <- cn0po01pc
    accuracy.table[i,87] <- cn0pot01pc
    accuracy.table[i,88] <- accur.cn0po01pc
    accuracy.table[i,89] <- cn2po02pc
    accuracy.table[i,90] <- cn2pot02pc
    accuracy.table[i,91] <- accur.cn2po02pc
    accuracy.table[i,92] <- cn1po02pc
    accuracy.table[i,93] <- cn1pot02pc
    accuracy.table[i,94] <- accur.cn1po02pc
    accuracy.table[i,95] <- cn0po02pc
    accuracy.table[i,96] <- cn0pot02pc
    accuracy.table[i,97] <- accur.cn0po02pc
    accuracy.table[i,98] <- cn2po11pc
    accuracy.table[i,99] <- cn2pot11pc
    accuracy.table[i,100] <- accur.cn2po11pc
    accuracy.table[i,101] <- cn1po11pc
    accuracy.table[i,102] <- cn1pot11pc
    accuracy.table[i,103] <- accur.cn1po11pc
    accuracy.table[i,104] <- cn0po11pc
    accuracy.table[i,105] <- cn0pot11pc
    accuracy.table[i,106] <- accur.cn0po11pc
    accuracy.table[i,107] <- cn2po12pc
    accuracy.table[i,108] <- cn2pot12pc
    accuracy.table[i,109] <- accur.cn2po12pc
    accuracy.table[i,110] <- cn1po12pc
    accuracy.table[i,111] <- cn1pot12pc
    accuracy.table[i,112] <- accur.cn1po12pc
    accuracy.table[i,113] <- cn0po12pc
    accuracy.table[i,114] <- cn0pot12pc
    accuracy.table[i,115] <- accur.cn0po12pc
    accuracy.table[i,116] <- cn2po22pc
    accuracy.table[i,117] <- cn2pot22pc
    accuracy.table[i,118] <- accur.cn2po22pc
    accuracy.table[i,119] <- cn1po22pc
    accuracy.table[i,120] <- cn1pot22pc
    accuracy.table[i,121] <- accur.cn1po22pc
    accuracy.table[i,122] <- cn0po22pc
    accuracy.table[i,123] <- cn0pot22pc
    accuracy.table[i,124] <- accur.cn0po22pc
    accuracy.table[i,125] <- cn2po00pno
    accuracy.table[i,126] <- cn2pot00pno
    accuracy.table[i,127] <- accur.cn2po00pno
    accuracy.table[i,128] <- cn1po00pno
    accuracy.table[i,129] <- cn1pot00pno
    accuracy.table[i,130] <- accur.cn1po00pno
    accuracy.table[i,131] <- cn0po00pno
    accuracy.table[i,132] <- cn0pot00pno
    accuracy.table[i,133] <- accur.cn0po00pno
    accuracy.table[i,134] <- cn2po01pno
    accuracy.table[i,135] <- cn2pot01pno
    accuracy.table[i,136] <- accur.cn2po01pno
    accuracy.table[i,137] <- cn1po01pno
    accuracy.table[i,138] <- cn1pot01pno
    accuracy.table[i,139] <- accur.cn1po01pno
    accuracy.table[i,140] <- cn0po01pno
    accuracy.table[i,141] <- cn0pot01pno
    accuracy.table[i,142] <- accur.cn0po01pno
    accuracy.table[i,143] <- cn2po02pno
    accuracy.table[i,144] <- cn2pot02pno
    accuracy.table[i,145] <- accur.cn2po02pno
    accuracy.table[i,146] <- cn1po02pno
    accuracy.table[i,147] <- cn1pot02pno
    accuracy.table[i,148] <- accur.cn1po02pno
    accuracy.table[i,149] <- cn0po02pno
    accuracy.table[i,150] <- cn0pot02pno
    accuracy.table[i,151] <- accur.cn0po02pno
    accuracy.table[i,152] <- cn2po11pno
    accuracy.table[i,153] <- cn2pot11pno
    accuracy.table[i,154] <- accur.cn2po11pno
    accuracy.table[i,155] <- cn1po11pno
    accuracy.table[i,156] <- cn1pot11pno
    accuracy.table[i,157] <- accur.cn1po11pno
    accuracy.table[i,158] <- cn0po11pno
    accuracy.table[i,159] <- cn0pot11pno
    accuracy.table[i,160] <- accur.cn0po11pno
    accuracy.table[i,161] <- cn2po12pno
    accuracy.table[i,162] <- cn2pot12pno
    accuracy.table[i,163] <- accur.cn2po12pno
    accuracy.table[i,164] <- cn1po12pno
    accuracy.table[i,165] <- cn1pot12pno
    accuracy.table[i,166] <- accur.cn1po12pno
    accuracy.table[i,167] <- cn0po12pno
    accuracy.table[i,168] <- cn0pot12pno
    accuracy.table[i,169] <- accur.cn0po12pno
    accuracy.table[i,170] <- cn2po22pno
    accuracy.table[i,171] <- cn2pot22pno
    accuracy.table[i,172] <- accur.cn2po22pno
    accuracy.table[i,173] <- cn1po22pno
    accuracy.table[i,174] <- cn1pot22pno
    accuracy.table[i,175] <- accur.cn1po22pno
    accuracy.table[i,176] <- cn0po22pno
    accuracy.table[i,177] <- cn0pot22pno
    accuracy.table[i,178] <- accur.cn0po22pno
    accuracy.table[i,179] <- cn2mbpo00pc
    accuracy.table[i,180] <- cn2mbpot00pc
    accuracy.table[i,181] <- accur.cn2mbpo00pc
    accuracy.table[i,182] <- cn1mbpo00pc
    accuracy.table[i,183] <- cn1mbpot00pc
    accuracy.table[i,184] <- accur.cn1mbpo00pc
    accuracy.table[i,185] <- cn0mbpo00pc
    accuracy.table[i,186] <- cn0mbpot00pc
    accuracy.table[i,187] <- accur.cn0mbpo00pc
    accuracy.table[i,188] <- cn2mbpo01pc
    accuracy.table[i,189] <- cn2mbpot01pc
    accuracy.table[i,190] <- accur.cn2mbpo01pc
    accuracy.table[i,191] <- cn1mbpo01pc
    accuracy.table[i,192] <- cn1mbpot01pc
    accuracy.table[i,193] <- accur.cn1mbpo01pc
    accuracy.table[i,194] <- cn0mbpo01pc
    accuracy.table[i,195] <- cn0mbpot01pc
    accuracy.table[i,196] <- accur.cn0mbpo01pc
    accuracy.table[i,197] <- cn2mbpo02pc
    accuracy.table[i,198] <- cn2mbpot02pc
    accuracy.table[i,199] <- accur.cn2mbpo02pc
    accuracy.table[i,200] <- cn1mbpo02pc
    accuracy.table[i,201] <- cn1mbpot02pc
    accuracy.table[i,202] <- accur.cn1mbpo02pc
    accuracy.table[i,203] <- cn0mbpo02pc
    accuracy.table[i,204] <- cn0mbpot02pc
    accuracy.table[i,205] <- accur.cn0mbpo02pc
    accuracy.table[i,206] <- cn2mbpo11pc
    accuracy.table[i,207] <- cn2mbpot11pc
    accuracy.table[i,208] <- accur.cn2mbpo11pc
    accuracy.table[i,209] <- cn1mbpo11pc
    accuracy.table[i,210] <- cn1mbpot11pc
    accuracy.table[i,211] <- accur.cn1mbpo11pc
    accuracy.table[i,212] <- cn0mbpo11pc
    accuracy.table[i,213] <- cn0mbpot11pc
    accuracy.table[i,214] <- accur.cn0mbpo11pc
    accuracy.table[i,215] <- cn2mbpo12pc
    accuracy.table[i,216] <- cn2mbpot12pc
    accuracy.table[i,217] <- accur.cn2mbpo12pc
    accuracy.table[i,218] <- cn1mbpo12pc
    accuracy.table[i,219] <- cn1mbpot12pc
    accuracy.table[i,220] <- accur.cn1mbpo12pc
    accuracy.table[i,221] <- cn0mbpo12pc
    accuracy.table[i,222] <- cn0mbpot12pc
    accuracy.table[i,223] <- accur.cn0mbpo12pc
    accuracy.table[i,224] <- cn2mbpo22pc
    accuracy.table[i,225] <- cn2mbpot22pc
    accuracy.table[i,226] <- accur.cn2mbpo22pc
    accuracy.table[i,227] <- cn1mbpo22pc
    accuracy.table[i,228] <- cn1mbpot22pc
    accuracy.table[i,229] <- accur.cn1mbpo22pc
    accuracy.table[i,230] <- cn0mbpo22pc
    accuracy.table[i,231] <- cn0mbpot22pc
    accuracy.table[i,232] <- accur.cn0mbpo22pc
    accuracy.table[i,233] <- cn2mbpo00pno
    accuracy.table[i,234] <- cn2mbpot00pno
    accuracy.table[i,235] <- accur.cn2mbpo00pno
    accuracy.table[i,236] <- cn1mbpo00pno
    accuracy.table[i,237] <- cn1mbpot00pno
    accuracy.table[i,238] <- accur.cn1mbpo00pno
    accuracy.table[i,239] <- cn0mbpo00pno
    accuracy.table[i,240] <- cn0mbpot00pno
    accuracy.table[i,241] <- accur.cn0mbpo00pno
    accuracy.table[i,242] <- cn2mbpo01pno
    accuracy.table[i,243] <- cn2mbpot01pno
    accuracy.table[i,244] <- accur.cn2mbpo01pno
    accuracy.table[i,245] <- cn1mbpo01pno
    accuracy.table[i,246] <- cn1mbpot01pno
    accuracy.table[i,247] <- accur.cn1mbpo01pno
    accuracy.table[i,248] <- cn0mbpo01pno
    accuracy.table[i,249] <- cn0mbpot01pno
    accuracy.table[i,250] <- accur.cn0mbpo01pno
    accuracy.table[i,251] <- cn2mbpo02pno
    accuracy.table[i,252] <- cn2mbpot02pno
    accuracy.table[i,253] <- accur.cn2mbpo02pno
    accuracy.table[i,254] <- cn1mbpo02pno
    accuracy.table[i,255] <- cn1mbpot02pno
    accuracy.table[i,256] <- accur.cn1mbpo02pno
    accuracy.table[i,257] <- cn0mbpo02pno
    accuracy.table[i,258] <- cn0mbpot02pno
    accuracy.table[i,259] <- accur.cn0mbpo02pno
    accuracy.table[i,260] <- cn2mbpo11pno
    accuracy.table[i,261] <- cn2mbpot11pno
    accuracy.table[i,262] <- accur.cn2mbpo11pno
    accuracy.table[i,263] <- cn1mbpo11pno
    accuracy.table[i,264] <- cn1mbpot11pno
    accuracy.table[i,265] <- accur.cn1mbpo11pno
    accuracy.table[i,266] <- cn0mbpo11pno
    accuracy.table[i,267] <- cn0mbpot11pno
    accuracy.table[i,268] <- accur.cn0mbpo11pno
    accuracy.table[i,269] <- cn2mbpo12pno
    accuracy.table[i,270] <- cn2mbpot12pno
    accuracy.table[i,271] <- accur.cn2mbpo12pno
    accuracy.table[i,272] <- cn1mbpo12pno
    accuracy.table[i,273] <- cn1mbpot12pno
    accuracy.table[i,274] <- accur.cn1mbpo12pno
    accuracy.table[i,275] <- cn0mbpo12pno
    accuracy.table[i,276] <- cn0mbpot12pno
    accuracy.table[i,277] <- accur.cn0mbpo12pno
    accuracy.table[i,278] <- cn2mbpo22pno
    accuracy.table[i,279] <- cn2mbpot22pno
    accuracy.table[i,280] <- accur.cn2mbpo22pno
    accuracy.table[i,281] <- cn1mbpo22pno
    accuracy.table[i,282] <- cn1mbpot22pno
    accuracy.table[i,283] <- accur.cn1mbpo22pno
    accuracy.table[i,284] <- cn0mbpo22pno
    accuracy.table[i,285] <- cn0mbpot22pno
    accuracy.table[i,286] <- accur.cn0mbpo22pno

  }

 # return accuracy.table
}
