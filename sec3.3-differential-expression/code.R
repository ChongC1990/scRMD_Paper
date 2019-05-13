setwd("C:/Users/CC/OneDrive/Algorithm")
library("kernlab")
library("Rmagic")
library("reticulate")
library("scImpute")
library("doParallel")
library("limma")
library("ggplot2")
source("Scimpute/simulation.R")
source("Scimpute/calculate_weight.R")
source("Scimpute/dmix.R")
source("Scimpute/get_mix_parameters.R")
source("Scimpute/imputation_model.R")
source("Scimpute/rmix.R")
source("function.R")

library(corpcor)
library(RSpectra)
library(expm)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(limma)
#library(Cairo)


HSC.data <- read.csv("RealData/GSE59114_DBA_GEO_all.csv", row.names = 1)
# from colnames to label
HSC.label <- colnames(HSC.data)
pos.uline <- gregexpr("_", HSC.label)
HSC.nobs <- length(HSC.label)
HSC.label.df <- data.frame(age = character(HSC.nobs), celltype = character(HSC.nobs), 
                           stringsAsFactors=FALSE)
for (i in 1:length(HSC.label)) {
  HSC.label.df$age[i] <- substr(HSC.label[[i]], 1, pos.uline[[i]][1] -1)
  HSC.label.df$celltype[i] <- substr(HSC.label[[i]], pos.uline[[i]][2] + 1, pos.uline[[i]][3] -1)
}

# focus on LTHSC
HSC.LT <- HSC.data[, HSC.label.df$celltype == "LTHSC"]
HSC.LT.label <- HSC.label.df$age[HSC.label.df$celltype == "LTHSC"]
HSC.LT.label <- HSC.LT.label[colSums(HSC.LT) > 1e4]
HSC.LT <- HSC.LT[, colSums(HSC.LT) > 1e4]
HSC.LT <- gene.filter(HSC.LT)
LTHSC.lcpm <- log10(apply(HSC.LT, 2, function(x) x / sum(x) * 1000000) + 1)
#LTHSC.lcpm <- gene.filter(LTHSC.lcpm)
LTHSC.coldata <- data.frame(condition = HSC.LT.label,row.names = colnames(HSC.LT))


# imputation methods
# setwd("C:/Users/pkuwc/Desktop/Study/Research/Single-cell/Code&Data/R")
# library("scImpute")
# source("Scimpute/get_mix_parameters.R")
# source("Scimpute/calculate_weight.R")
# source("Scimpute/dmix.R")
# source("Scimpute/rmix.R")
# source("Scimpute/imputation_model.R")
# plist = get_mix_parameters(LTHSC.lcpm + log10(1.01))
# res.sci = imputation_model1(count = LTHSC.lcpm + log10(1.01), point = log10(1.01), plist, drop_thre = 0.5)
# res.sci = res.sci - log10(1.01)
# 
# write.csv(res.sci, "LTHSC.scimpute.csv")

# no need to run scimpute, read results directly
# LTHSC.res.sci = imputation_model8(count = LTHSC.lcpm+log10(1.01), labeled = False, 
#   			point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = 10, out_dir = "./")$count_imp
# LTHSC.res.sci = LTHSC.res.sci - log10(1.01)
LTHSC.tpm <- apply(HSC.LT, 2, function(x) x / sum(x) * 1000000)
#write.csv(LTHSC.tpm, file = "RealData/LTHSC_tpm.csv")
#scimpute(count_path = "RealData/LTHSC_tpm.csv", out_dir = "sci_res/", ncores = 1, Kcluster = 2)
#LTHSC.res.sci_full <- log10(read.csv("sci_res/scimpute_count.csv", row.names = 1) + 1 )
LTHSC.res.sci_full <- imputation_model8(count = LTHSC.lcpm + log10(1.01), labeled = FALSE, 
         point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = 2, out_dir = "./")$count_imp
LTHSC.res.sci_full = LTHSC.res.sci_full - log10(1.01) 
LTHSC.res.rmd_full <- rmd(LTHSC.lcpm, candidate = quantile(LTHSC.lcpm[LTHSC.lcpm > 0], 0.05))
LTHSC.res.magic_full <- t(magic(t(LTHSC.lcpm), genes="all_genes")$result)

limma_test <- function(df, foldchange = 1.5, label = LTHSC.coldata){
  if (class(label) == "data.frame"){
    label = label[,1]
  }
  old_mean = rowMeans(df[, label == "old"])
  young_mean = rowMeans(df[, label == "young"])
  folds = rbind(old_mean/young_mean, young_mean/old_mean)
  fold_changes = apply(folds, 2, max)
  limma.res <- lmFit(df, design = model.matrix(~label))
  limma.res <- eBayes(limma.res)
  limma.padj <- p.adjust(limma.res$p.value[,2], method = "fdr")
  limma.padj[fold_changes < foldchange] = 1
  return (limma.padj)
}

# limma检验！！！
#foldchange <- 0
{
  
  # LTHSC.sc.absfc <- abs(rowMeans(LTHSC.lcpm[,LTHSC.coldata[,1] == "old"]) - 
  #                         rowMeans(LTHSC.lcpm[,LTHSC.coldata[,1] == "young"]))
  # LTHSC.raw.limma.res <- lmFit(LTHSC.lcpm, design = model.matrix(~LTHSC.coldata[,1]))
  # LTHSC.raw.limma.res <- eBayes(LTHSC.raw.limma.res)
  # LTHSC.raw.limma.padj <- p.adjust(LTHSC.raw.limma.res$p.value[,2], method = "fdr")
  # LTHSC.raw.limma.padj[LTHSC.sc.absfc < log10(foldchange)] <- 1
  
  # LTHSC.sci.limma.res <- lmFit(LTHSC.res.sci, design = model.matrix(~LTHSC.coldata[,1]))
  # LTHSC.sci.limma.res <- eBayes(LTHSC.sci.limma.res)
  # LTHSC.sci.limma.padj <- p.adjust(LTHSC.sci.limma.res$p.value[,2], method = "fdr")
  # LTHSC.sci.limma.padj[LTHSC.sc.absfc < log10(foldchange)] <- 1
  
  # LTHSC.rmd.limma.res <- lmFit(LTHSC.res.rmd$exprs, design = model.matrix(~LTHSC.coldata[,1]))
  # LTHSC.rmd.limma.res <- eBayes(LTHSC.rmd.limma.res)
  # LTHSC.rmd.limma.padj <- p.adjust(LTHSC.rmd.limma.res$p.value[,2], method = "fdr")
  # LTHSC.rmd.limma.padj[LTHSC.sc.absfc < log10(foldchange)] <- 1
  
  # LTHSC.magic.limma.res <- lmFit(LTHSC.res.magic, design = model.matrix(~LTHSC.coldata[,1]))
  # LTHSC.magic.limma.res <- eBayes(LTHSC.magic.limma.res)
  # LTHSC.magic.limma.padj <- p.adjust(LTHSC.magic.limma.res$p.value[,2], method = "fdr")
  # LTHSC.magic.limma.padj[LTHSC.sc.absfc < log10(foldchange)] <- 1
  LTHSC.raw.p.adj_fc15 <- limma_test(df = LTHSC.lcpm, foldchange = 1.5)
  LTHSC.sci.p.adj_fc15 <- limma_test(df = LTHSC.res.sci_full, foldchange = 1.5)
  LTHSC.rmd.p.adj_fc15 <- limma_test(df = LTHSC.res.rmd_full$exprs, foldchange = 1.5)
  LTHSC.magic.p.adj_fc15 <- limma_test(df = LTHSC.res.magic_full, foldchange = 1.5)

  LTHSC.raw.p.adj_fc12 <- limma_test(df = LTHSC.lcpm, foldchange = 1.2)
  LTHSC.sci.p.adj_fc12 <- limma_test(df = LTHSC.res.sci_full, foldchange = 1.2)
  LTHSC.rmd.p.adj_fc12 <- limma_test(df = LTHSC.res.rmd_full$exprs, foldchange = 1.2)
  LTHSC.magic.p.adj_fc12 <- limma_test(df = LTHSC.res.magic_full, foldchange = 1.2)

  LTHSC.raw.p.adj_fc10 <- limma_test(df = LTHSC.lcpm, foldchange = 1.0)
  LTHSC.sci.p.adj_fc10 <- limma_test(df = LTHSC.res.sci_full, foldchange = 1.0)
  LTHSC.rmd.p.adj_fc10 <- limma_test(df = LTHSC.res.rmd_full$exprs, foldchange = 1.0)
  LTHSC.magic.p.adj_fc10 <- limma_test(df = LTHSC.res.magic_full, foldchange = 1.0)

}

test.eval <- function(p1, p2, top = FALSE, thresh = 0.05){
  if (top == FALSE){
    r_pos <- names(p1[p1 <= thresh])
    r_neg <- names(p1[p1 > thresh])
  }
  else{
    r_pos <- names(sort(p1)[1:top])
    r_neg <- names(sort(p1)[(top+1):length(p1)])
  }
  # p_pos <- names(sort(p2)[1:top])
  # p_neg <- names(sort(p2)[(top+1):length(p1)])
  p_pos <- names(p2[p2 <= thresh])
  p_neg <- names(p2[p2 > thresh])
  tp = length(intersect(r_pos, p_pos))
  tn = length(intersect(r_neg, p_neg))
  fp = length(setdiff(p_pos, r_pos))
  fn = length(setdiff(p_neg, r_neg))
  precision <- tp/(tp + fp)
  recall <- tp/(tp + fn)
  Fscore <- 2 * tp / (2 * tp + fp + fn)
  return(list(TP = tp, Fscore = Fscore, TPR = recall, precision = precision))
}
{
# limma reproducibility！！！执行需要花一点点时间，画原图的话可以跳过这一步，直接读表，已经存在文件里了
  #set.seed(2019)
  thresh <- 0.05
  npercentage <- 8; ntrial <- 50; nobs <- npercentage * ntrial * 4; ngold <- 1000;
  # res.reproduce <- data.frame(method = character(nobs), percentage = numeric(nobs), trial = numeric(nobs),
  #                             detection = numeric(nobs), precision = numeric(nobs), recall = numeric(nobs), 
  #                             Fscore = numeric(nobs), MCC = numeric(nobs), stringsAsFactors = F)
  res.reproduce_fc15 <- data.frame(method = character(nobs), percentage = numeric(nobs), trial = numeric(nobs),
                              detection = numeric(nobs), precision = numeric(nobs), recall = numeric(nobs), Fscore = numeric(nobs), stringsAsFactors = F)
  res.reproduce_fc12 <- data.frame(method = character(nobs), percentage = numeric(nobs), trial = numeric(nobs),
                              detection = numeric(nobs), precision = numeric(nobs), recall = numeric(nobs), Fscore = numeric(nobs), stringsAsFactors = F)
  res.reproduce_fc10 <- data.frame(method = character(nobs), percentage = numeric(nobs), trial = numeric(nobs),
                              detection = numeric(nobs), precision = numeric(nobs), recall = numeric(nobs), Fscore = numeric(nobs), stringsAsFactors = F)
  young_index = which(HSC.LT.label == "young")
  old_index = which(HSC.LT.label == "old")
  for (i in 1:npercentage){
    print(i)
    percentage = 10 - i
    sub.young.n <- ceiling(percentage/10 * sum(HSC.LT.label == "young"))
    sub.old.n <- ceiling(percentage/10 * sum(HSC.LT.label == "old"))
    for (k in 1:ntrial) {
      set.seed(k)
      sub.young.idx <- sample(young_index, sub.young.n)
      sub.old.idx <- sample(old_index, sub.old.n)
      df_raw = LTHSC.lcpm[,c(sub.young.idx,sub.old.idx)]
      df_label = LTHSC.coldata[c(sub.young.idx,sub.old.idx), ]
      df_magic = t(magic(t(df_raw), genes="all_genes", npca=min(100, dim(df_raw)[2]))$result)
      df_rmd = rmd(df_raw, candidate = quantile(df_raw[df_raw > 0], 0.05))
      df_sci = imputation_model8(count = df_raw + log10(1.01), labeled = FALSE, 
         point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = 2, out_dir = "./")$count_imp
      df_sci = df_sci - log10(1.01)   
      # raw.sub.test <- lmFit(LTHSC.lcpm[,c(sub.young.idx,sub.old.idx)], design = model.matrix(~c(rep("young",sub.young.n),rep("old",sub.old.n))))
      # raw.sub.test <- eBayes(raw.sub.test)
      # raw.sub.p <- p.adjust(raw.sub.test$p.value[,2], method = "fdr")
      # rmd.sub.test <- lmFit(LTHSC.res.rmd$exprs[,c(sub.young.idx,sub.old.idx)], design = model.matrix(~c(rep("young",sub.young.n),rep("old",sub.old.n))))
      # rmd.sub.test <- eBayes(rmd.sub.test)
      # rmd.sub.p <- p.adjust(rmd.sub.test$p.value[,2], method = "fdr")
      # sci.sub.test <- lmFit(LTHSC.res.sci[,c(sub.young.idx,sub.old.idx)], design = model.matrix(~c(rep("young",sub.young.n),rep("old",sub.old.n))))
      # sci.sub.test <- eBayes(sci.sub.test)
      # sci.sub.p <- p.adjust(sci.sub.test$p.value[,2], method = "fdr")
      # magic.sub.test <- lmFit(LTHSC.res.magic[,c(sub.young.idx,sub.old.idx)], design = model.matrix(~c(rep("young",sub.young.n),rep("old",sub.old.n))))
      # magic.sub.test <- eBayes(magic.sub.test)
      # magic.sub.p <- p.adjust(magic.sub.test$p.value[,2], method = "fdr")
      
      # raw.test.eval <- test.eval(LTHSC.raw.p.adj, raw.sub.p, sort(LTHSC.raw.p.adj)[ngold], sort(raw.sub.p)[ngold])
      # rmd.test.eval <- test.eval(LTHSC.rmd.p.adj, rmd.sub.p, sort(LTHSC.rmd.p.adj)[ngold], sort(rmd.sub.p)[ngold])
      # sci.test.eval <- test.eval(LTHSC.sci.p.adj, sci.sub.p, sort(LTHSC.sci.p.adj)[ngold], sort(sci.sub.p)[ngold])
      # magic.test.eval <- test.eval(LTHSC.magic.p.adj, magic.sub.p, sort(LTHSC.magic.p.adj)[ngold], sort(magic.sub.p)[ngold])
      raw.sub.p_fc15 <- limma_test(df = df_raw, label = df_label, foldchange = 1.5)
      rmd.sub.p_fc15 <- limma_test(df = df_rmd$exprs, label = df_label, foldchange = 1.5)
      sci.sub.p_fc15 <- limma_test(df = df_sci, label = df_label, foldchange = 1.5)
      magic.sub.p_fc15 <- limma_test(df = df_magic, label = df_label, foldchange = 1.5)

      raw.sub.p_fc12 <- limma_test(df = df_raw, label = df_label, foldchange = 1.2)
      rmd.sub.p_fc12 <- limma_test(df = df_rmd$exprs, label = df_label, foldchange = 1.2)
      sci.sub.p_fc12 <- limma_test(df = df_sci, label = df_label, foldchange = 1.2)
      magic.sub.p_fc12 <- limma_test(df = df_magic, label = df_label, foldchange = 1.2)

      raw.sub.p_fc10 <- limma_test(df = df_raw, label = df_label, foldchange = 1.0)
      rmd.sub.p_fc10 <- limma_test(df = df_rmd$exprs, label = df_label, foldchange = 1.0)
      sci.sub.p_fc10 <- limma_test(df = df_sci, label = df_label, foldchange = 1.0)
      magic.sub.p_fc10 <- limma_test(df = df_magic, label = df_label, foldchange = 1.0)

      raw.test.eval_fc15 <- test.eval(LTHSC.raw.p.adj_fc15, raw.sub.p_fc15)
      rmd.test.eval_fc15 <- test.eval(LTHSC.rmd.p.adj_fc15, rmd.sub.p_fc15)
      sci.test.eval_fc15 <- test.eval(LTHSC.sci.p.adj_fc15, sci.sub.p_fc15)
      magic.test.eval_fc15 <- test.eval(LTHSC.magic.p.adj_fc15, magic.sub.p_fc15)

      raw.test.eval_fc12 <- test.eval(LTHSC.raw.p.adj_fc12, raw.sub.p_fc12)
      rmd.test.eval_fc12 <- test.eval(LTHSC.rmd.p.adj_fc12, rmd.sub.p_fc12)
      sci.test.eval_fc12 <- test.eval(LTHSC.sci.p.adj_fc12, sci.sub.p_fc12)
      magic.test.eval_fc12 <- test.eval(LTHSC.magic.p.adj_fc12, magic.sub.p_fc12)

      raw.test.eval_fc10 <- test.eval(LTHSC.raw.p.adj_fc10, raw.sub.p_fc10)
      rmd.test.eval_fc10 <- test.eval(LTHSC.rmd.p.adj_fc10, rmd.sub.p_fc10)
      sci.test.eval_fc10 <- test.eval(LTHSC.sci.p.adj_fc10, sci.sub.p_fc10)
      magic.test.eval_fc10 <- test.eval(LTHSC.magic.p.adj_fc10, magic.sub.p_fc10)

      res.reproduce_fc15$method[ntrial * (i-1) * 4 + (k-1) * 4 + 1] <- "Raw"
      res.reproduce_fc15$method[ntrial * (i-1) * 4 + (k-1) * 4 + 2] <- "RMD"
      res.reproduce_fc15$method[ntrial * (i-1) * 4 + (k-1) * 4 + 3] <- "SCI"
      res.reproduce_fc15$method[ntrial * (i-1) * 4 + (k-1) * 4 + 4] <- "MAGIC"

      res.reproduce_fc12$method[ntrial * (i-1) * 4 + (k-1) * 4 + 1] <- "Raw"
      res.reproduce_fc12$method[ntrial * (i-1) * 4 + (k-1) * 4 + 2] <- "RMD"
      res.reproduce_fc12$method[ntrial * (i-1) * 4 + (k-1) * 4 + 3] <- "SCI"
      res.reproduce_fc12$method[ntrial * (i-1) * 4 + (k-1) * 4 + 4] <- "MAGIC"

      res.reproduce_fc10$method[ntrial * (i-1) * 4 + (k-1) * 4 + 1] <- "Raw"
      res.reproduce_fc10$method[ntrial * (i-1) * 4 + (k-1) * 4 + 2] <- "RMD"
      res.reproduce_fc10$method[ntrial * (i-1) * 4 + (k-1) * 4 + 3] <- "SCI"
      res.reproduce_fc10$method[ntrial * (i-1) * 4 + (k-1) * 4 + 4] <- "MAGIC"
      
      res.reproduce_fc15[ntrial*(i-1)*4 + (k-1)*4 + 1, -1] <- c(percentage, k, sum(raw.sub.p_fc15<thresh), raw.test.eval_fc15$precision,
                                                           raw.test.eval_fc15$TPR, raw.test.eval_fc15$Fscore)
      res.reproduce_fc15[ntrial*(i-1)*4 + (k-1)*4 + 2, -1] <- c(percentage, k, sum(rmd.sub.p_fc15<thresh), rmd.test.eval_fc15$precision,
                                                           rmd.test.eval_fc15$TPR, rmd.test.eval_fc15$Fscore)
      res.reproduce_fc15[ntrial*(i-1)*4 + (k-1)*4 + 3, -1] <- c(percentage, k, sum(sci.sub.p_fc15<thresh), sci.test.eval_fc15$precision,
                                                           sci.test.eval_fc15$TPR, sci.test.eval_fc15$Fscore)
      res.reproduce_fc15[ntrial*(i-1)*4 + (k-1)*4 + 4, -1] <- c(percentage, k, sum(magic.sub.p_fc15<thresh), magic.test.eval_fc15$precision,
                                                           magic.test.eval_fc15$TPR, magic.test.eval_fc15$Fscore)

      res.reproduce_fc12[ntrial*(i-1)*4 + (k-1)*4 + 1, -1] <- c(percentage, k, sum(raw.sub.p_fc12<thresh), raw.test.eval_fc12$precision,
                                                           raw.test.eval_fc12$TPR, raw.test.eval_fc12$Fscore)
      res.reproduce_fc12[ntrial*(i-1)*4 + (k-1)*4 + 2, -1] <- c(percentage, k, sum(rmd.sub.p_fc12<thresh), rmd.test.eval_fc12$precision,
                                                           rmd.test.eval_fc12$TPR, rmd.test.eval_fc12$Fscore)
      res.reproduce_fc12[ntrial*(i-1)*4 + (k-1)*4 + 3, -1] <- c(percentage, k, sum(sci.sub.p_fc12<thresh), sci.test.eval_fc12$precision,
                                                           sci.test.eval_fc12$TPR, sci.test.eval_fc12$Fscore)
      res.reproduce_fc12[ntrial*(i-1)*4 + (k-1)*4 + 4, -1] <- c(percentage, k, sum(magic.sub.p_fc12<thresh), magic.test.eval_fc12$precision,
                                                           magic.test.eval_fc12$TPR, magic.test.eval_fc12$Fscore)

      res.reproduce_fc10[ntrial*(i-1)*4 + (k-1)*4 + 1, -1] <- c(percentage, k, sum(raw.sub.p_fc10<thresh), raw.test.eval_fc10$precision,
                                                           raw.test.eval_fc10$TPR, raw.test.eval_fc10$Fscore)
      res.reproduce_fc10[ntrial*(i-1)*4 + (k-1)*4 + 2, -1] <- c(percentage, k, sum(rmd.sub.p_fc10<thresh), rmd.test.eval_fc10$precision,
                                                           rmd.test.eval_fc10$TPR, rmd.test.eval_fc10$Fscore)
      res.reproduce_fc10[ntrial*(i-1)*4 + (k-1)*4 + 3, -1] <- c(percentage, k, sum(sci.sub.p_fc10<thresh), sci.test.eval_fc10$precision,
                                                           sci.test.eval_fc10$TPR, sci.test.eval_fc10$Fscore)
      res.reproduce_fc10[ntrial*(i-1)*4 + (k-1)*4 + 4, -1] <- c(percentage, k, sum(magic.sub.p_fc10<thresh), magic.test.eval_fc10$precision,
                                                           magic.test.eval_fc10$TPR, magic.test.eval_fc10$Fscore)
    }
  }
  write.csv(res.reproduce_fc15, "Result/LTHSC-limma-reproducibility_fc1.5.csv")
  write.csv(res.reproduce_fc12, "Result/LTHSC-limma-reproducibility_fc1.2.csv")
  write.csv(res.reproduce_fc10, "Result/LTHSC-limma-reproducibility_fc1.0.csv")
}
# D图的绘制！！！
# plot reproducibility
scRMD.col = "#e41a1c";raw.col = "#377eb8";sci.col = "#4daf4a";magic.col = "#984ea3";te.col = "#ff7f00"
{
  npercentage <- 8;
  LTHSC.res.reproduce <- read.csv("Result/LTHSC-limma-reproducibility_fc1.5.csv",
                                  row.names = 1)
  LTHSC.raw.reproduce <- sapply(1:npercentage, function(x) {
    colMeans(LTHSC.res.reproduce[LTHSC.res.reproduce$method == "Raw" & (LTHSC.res.reproduce$percentage == 10-x), -c(1,3) ], na.rm = T)
  } )
  LTHSC.rmd.reproduce <- sapply(1:npercentage, function(x) {
    colMeans(LTHSC.res.reproduce[LTHSC.res.reproduce$method == "RMD" & LTHSC.res.reproduce$percentage == 10-x, -c(1,3) ], na.rm = T)
  } )
  LTHSC.sci.reproduce <- sapply(1:npercentage, function(x) {
    colMeans(LTHSC.res.reproduce[LTHSC.res.reproduce$method == "SCI" & LTHSC.res.reproduce$percentage == 10-x, -c(1,3) ], na.rm = T)
  } )
  LTHSC.magic.reproduce <- sapply(1:npercentage, function(x) {
    colMeans(LTHSC.res.reproduce[LTHSC.res.reproduce$method == "MAGIC" & LTHSC.res.reproduce$percentage == 10-x, -c(1,3) ], na.rm = T)
  } )
  
  LTHSC.res.reproduce <- data.frame(Method = rep(c("Raw","scRMD","scImpute","MAGIC"), each = npercentage),
                                    t(cbind(LTHSC.raw.reproduce, LTHSC.rmd.reproduce, LTHSC.sci.reproduce, LTHSC.magic.reproduce)))
  LTHSC.res.reproduce$Method <- factor(LTHSC.res.reproduce$Method, levels = c("scRMD", "scImpute", "MAGIC", "Raw"))
  LTHSC.detection.plot <- ggplot(data = LTHSC.res.reproduce, mapping = aes(x = percentage / 10, y = detection)) + 
    geom_line(aes(colour = Method), size = 1) +  scale_y_log10(breaks = c(10,100,1000,10000)) + 
    scale_x_reverse( limits =c(0.9,0.2), breaks = 9:2/10) + 
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15), 
          plot.title = element_text(hjust = 0.5, size = 25),
          legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + 
    ggtitle("Number of detections") +scale_color_manual(values=c(scRMD.col,sci.col,magic.col,raw.col))+ 
    theme(panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
          axis.line.y = element_line(colour = "black", size = 1),
          legend.key = element_blank(), strip.background = element_blank(),
          legend.position = c(0.1,0.22),legend.justification = c("left"),
          legend.background = element_rect(fill="transparent")) +
    xlab("") + ylab("Value")
  LTHSC.detection.plot
  
  LTHSC.res.reproduce.nomagic <- LTHSC.res.reproduce[LTHSC.res.reproduce$Method != "MAGIC", ]
  LTHSC.precision.plot <- ggplot(data = LTHSC.res.reproduce.nomagic, mapping = aes(x = percentage / 10, y = precision)) + 
    geom_line(aes(colour = Method), size = 1)  + 
    scale_x_reverse( limits =c(0.9,0.2), breaks = 9:2/10) + 
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15), 
          plot.title = element_text(hjust = 0.5, size = 25)) + 
    ggtitle("Precision") +scale_color_manual(values=c(scRMD.col,sci.col,raw.col))+ 
    theme(panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
          axis.line.y = element_line(colour = "black", size = 1),
          legend.key = element_blank(), strip.background = element_blank(),
          legend.position = "none") +
    xlab("Size of subset") + ylab("")
  LTHSC.precision.plot
  
  LTHSC.recall.plot <- ggplot(data = LTHSC.res.reproduce.nomagic, mapping = aes(x = percentage / 10, y = recall)) + 
    geom_line(aes(colour = Method), size = 1)  + 
    scale_x_reverse( limits =c(0.9,0.2), breaks = 9:2/10) + 
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15), 
          plot.title = element_text(hjust = 0.5, size = 25)) + 
    ggtitle("Recall") +scale_color_manual(values=c(scRMD.col,sci.col,raw.col))+ 
    theme(panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
          axis.line.y = element_line(colour = "black", size = 1),
          legend.key = element_blank(), strip.background = element_blank(),
          legend.position = "none") +
    xlab("") + ylab("")
  LTHSC.recall.plot
  
  LTHSC.Fscore.plot <- ggplot(data = LTHSC.res.reproduce.nomagic, mapping = aes(x = percentage / 10, y = Fscore)) + 
    geom_line(aes(colour = Method), size = 1)  + 
    scale_x_reverse( limits =c(0.9,0.2), breaks = 9:2/10) + 
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15), 
          plot.title = element_text(hjust = 0.5, size = 25)) + 
    ggtitle("F-score") +scale_color_manual(values=c(scRMD.col,sci.col,raw.col))+ 
    theme(panel.background = element_blank(), axis.line.x = element_line(colour = "black", size = 1),
          axis.line.y = element_line(colour = "black", size = 1),
          legend.key = element_blank(), strip.background = element_blank(),
          legend.position = "none") +
    xlab("") + ylab("")
  LTHSC.Fscore.plot
  
  
  cairo_pdf("LTHSC-reproducibility_fc1.5.pdf", width = 16, height = 4)
  multiplot(LTHSC.detection.plot, LTHSC.precision.plot, LTHSC.recall.plot, LTHSC.Fscore.plot, cols = 4)
  dev.off()
}