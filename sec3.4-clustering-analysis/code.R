setwd("/data1/XiLab/chenchong/scRMD")
library("kernlab")
library("Rmagic")
library("reticulate")
library("scImpute")
library("data.table")
library("doParallel")
library("cidr")
library("SIMLR")
source("Scimpute/simulation.R")
source("Scimpute/calculate_weight.R")
source("Scimpute/dmix.R")
source("Scimpute/get_mix_parameters.R")
source("Scimpute/imputation_model.R")
source("Scimpute/rmix.R")
source("function.R")



################## Basic Functions ####
subsample <- function(A,j,prob){# A is p by n
  print(j)
  p <- ncol(A)
  sample.mtx <- rep(0,p)
  nonzero.idx <- A[j,] != 0
  sample.mtx[nonzero.idx] <- rbinom(sum(nonzero.idx), as.matrix(round(A[j,nonzero.idx])), prob)
  return(sample.mtx)
}

GeneNorm <- function(Y){
  Y = Y/rowSums(Y)*1000000
  return (log10(Y+1))
}


############### Compare Clustering Result ########################################
candidate_thres <- c(0, 0.02, 0.05, 0.1, 0.2)
datasets <-c("Usoskin","Ting","Pollen","Deng")
methods <- c("raw","rmd","sci","mag")
repeats <- 10
set.seed(2017)
for (cid in 2:length(datasets)){
	path <- paste("RealData/",datasets[cid],"_RAW.txt",sep="")
	RAW <- read.table(path,header=F)
	true_labs <- RAW[1,]
	RAW <- t(RAW[-1,])
	Counts <- round(10^RAW-1)
	N=dim(RAW)[1];P=dim(RAW)[2]
	K = max(as.integer(true_labs))
	#### Full Data ###
	res.raw <- RAW
    fwrite(data.table(t(res.raw)), paste("real_data/raw_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)
  	res.sci = imputation_model8(count = t(res.raw+log10(1.01)), labeled = False, 
  			point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  	res.sci = t(res.sci - log10(1.01))
    fwrite(data.table(t(res.sci)), paste("real_data/sci_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)
    for (qq in candidate_thres){
    	res.rmd <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw>0], qq))$exprs
    	fwrite(data.table(t(res.rmd)), paste("real_data/rmd_", datasets[cid], "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    } 		
	res.mag <-magic(res.raw, genes="all_genes")$result
    fwrite(data.table(t(res.mag)), paste("real_data/mag_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)	
	#### Drop Model ###
	Pro <- exp(-0.2*RAW*RAW)
	for(i in 1:repeats){
		res.raw <- RAW * matrix(rbinom(N*P, 1, as.matrix(1 - Pro)), N, P)
		res.raw <- res.raw[,colSums(res.raw)!=0]
    	fwrite(data.table(t(res.raw)), paste("real_data/drop_raw_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)
  		res.sci = imputation_model8(count = t(res.raw+log10(1.01)), labeled = False, 
  				point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  		res.sci = t(res.sci - log10(1.01))
    	fwrite(data.table(t(res.sci)), paste("real_data/drop_sci_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    	for (qq in candidate_thres){
    		res.rmd <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw>0], qq))$exprs
    		fwrite(data.table(t(res.rmd)), paste("real_data/drop_rmd_", datasets[cid], "_rep", i, "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    	} 		
		res.mag <-magic(res.raw, genes="all_genes")$result
    	fwrite(data.table(t(res.mag)), paste("real_data/drop_mag_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)	
	}

	#### Down Sampling ####
	for (i in 1:repeats){
		num.cores = 10
		registerDoParallel(cores=num.cores)
		subcounts <- Counts
		n = dim(Counts)[1]
		result <- foreach(j=1:n)%dopar%subsample(Counts,j,0.1)
		for (j in 1:n){
			subcounts[j,] = result[[j]]
		}
		res.raw <- GeneNorm(subcounts)
		res.raw <- res.raw[,colSums(res.raw)!=0]
    	fwrite(data.table(t(res.raw)), paste("real_data/down_raw_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)
  		res.sci = imputation_model8(count = t(res.raw+log10(1.01)), labeled = False, 
  				point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  		res.sci = t(res.sci - log10(1.01))
    	fwrite(data.table(t(res.sci)), paste("real_data/down_sci_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    	for (qq in candidate_thres){
    		res.rmd <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw>0], qq))$exprs
    		fwrite(data.table(t(res.rmd)), paste("real_data/down_rmd_", datasets[cid], "_rep", i, "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    	} 		
		res.mag <-magic(res.raw, genes="all_genes")$result
    	fwrite(data.table(t(res.mag)), paste("real_data/down_mag_", datasets[cid], "_rep", i, ".txt", sep=""),row.names = F, col.names = F, quote = F)	
	}
	########### Low Expression ###
	MeanExpression <- colSums(RAW)/colSums(RAW!=0)
	res.raw <- RAW[,MeanExpression<quantile(MeanExpression)[2]]
    fwrite(data.table(t(res.raw)), paste("real_data/low_raw_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)
  	res.sci = imputation_model8(count = t(res.raw+log10(1.01)), labeled = False, 
  			point = log10(1.01), drop_thre = 0.5, ncores = 1, Kcluster = K, out_dir = "./")$count_imp
  	res.sci = t(res.sci - log10(1.01))
    fwrite(data.table(t(res.sci)), paste("real_data/low_sci_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)
    for (qq in candidate_thres){
    	res.rmd <- rmd(as.matrix(res.raw), candidate = quantile(res.raw[res.raw>0], qq))$exprs
    	fwrite(data.table(t(res.rmd)), paste("real_data/low_rmd_", datasets[cid], "_candidate", qq, ".txt", sep=""),row.names = F, col.names = F, quote = F)
    } 		
	res.mag <-magic(res.raw, genes="all_genes")$result
    fwrite(data.table(t(res.mag)), paste("real_data/low_mag_", datasets[cid], ".txt", sep=""),row.names = F, col.names = F, quote = F)	
}


##################### plot the result #######################
library("ggplot2")

colors = c("#ff7f00","#377eb8","#e41a1c","#4daf4a","#984ea3")
methods <- c("raw","rmd","sci","mag")
METHODs <- c("RAW", "scRMD", "scImpute", "MAGIC")
quantiles <- c(0, 0.02, 0.05, 0.1, 0.2)
datasets <-c("Usoskin","Ting","Pollen","Deng")
repeats = 10

################## Estimated K ##################################
datasets <-c("Deng")
for (data in datasets){
	path <- paste("RealData/",data,"_RAW.txt",sep="")
	RAW <- read.table(path,header=F)
	true_labs <- RAW[1,]
	#K = max(as.integer(true_labs))
	fulldata.df <- data.frame()
	for (i in c(1:4)){
		res_tmp = data.frame(model = rep("FullData", 4), method = rep(METHODs[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = numeric(4), nclustsd = rep(0, 4), ari = numeric(4), arisd = rep(0, 4))
		if (i ==2){
			df = t(read.table(paste("real_data/", methods[i], "_", data, "_candidate0.05.txt", sep=""), header = F, sep = ","))
		}
		else{
			df = t(read.table(paste("real_data/", methods[i], "_", data, ".txt", sep=""), header = F, sep = ","))
		}
		if (data == "Usoskin"){
			tags = t(round(2^df-1))
		}
		else{
			tags = t(round(10^df-1))
		}
		sData <- CIDR(tags)
		K = sData@nCluster
		res_tmp$nclust = K
		res_tmp$ari[1] = adjustedRandIndex(true_labs, sData@clusters)
		res_tmp$ari[2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
		res_tmp$ari[3] = adjustedRandIndex(true_labs, tSNE(df, K))
		res_tmp$ari[4] = adjustedRandIndex(true_labs, SparsePca(df, K))
		fulldata.df <- rbind(fulldata.df, res_tmp)
	}

	dropmodel.df <- data.frame()
	for (i in c(1:3)){
		ari <- matrix(0, repeats, 4)
		nClust <- matrix(0, repeats, 4)
		for (j in c(1:repeats)){
			if (i ==2){
				df = t(read.table(paste("real_data/drop_", methods[i], "_", data, "_rep", j, "_candidate0.05.txt", sep=""), header = F, sep = ","))
			}
			else{
				df = t(read.table(paste("real_data/drop_", methods[i], "_", data, "_rep", j, ".txt", sep=""), header = F, sep = ","))
			}
			if (data == "Usoskin"){
				tags = t(round(2^df-1))
			}
			else{
				tags = t(round(10^df-1))
			}
			sData <- CIDR(tags)
			K = sData@nCluster
			nClust[j,] = K
			ari[j,1] = adjustedRandIndex(true_labs, sData@clusters)
			ari[j,2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
			ari[j,3] = adjustedRandIndex(true_labs, tSNE(df, K))
			ari[j,4] = adjustedRandIndex(true_labs, SparsePca(df, K))	
		}
		res_tmp = data.frame(model = rep("DropModel", 4), method = rep(METHODs[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = apply(nClust, 2, mean), nclustsd = apply(nClust, 2, sd), ari = apply(ari, 2, mean), arisd = apply(ari, 2, sd))
		dropmodel.df <- rbind(dropmodel.df, res_tmp)
	}

	downsampling.df <- data.frame()
	for (i in c(1:3)){
		ari <- matrix(0, repeats, 4)
		nClust <- matrix(0, repeats, 4)
		for (j in c(1:repeats)){
			if (i ==2){
				df = t(read.table(paste("real_data/down_", methods[i], "_", data, "_rep", j, "_candidate0.05.txt", sep=""), header = F, sep = ","))
			}
			else{
				df = t(read.table(paste("real_data/down_", methods[i], "_", data, "_rep", j, ".txt", sep=""), header = F, sep = ","))
			}
			tags = t(round(10^df-1))
			sData <- CIDR(tags)
			K = sData@nCluster
			nClust[j,] = K
			ari[j,1] = adjustedRandIndex(true_labs, sData@clusters)
			ari[j,2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
			ari[j,3] = adjustedRandIndex(true_labs, tSNE(df, K))
			ari[j,4] = adjustedRandIndex(true_labs, SparsePca(df, K))	
		}
		res_tmp = data.frame(model = rep("DownSampling", 4), method = rep(METHODs[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = apply(nClust, 2, mean), nclustsd = apply(nClust, 2, sd), ari = apply(ari, 2, mean), arisd = apply(ari, 2, sd))
		downsampling.df <- rbind(downsampling.df, res_tmp)
	}

	lowexpression.df <- data.frame()
	for (i in c(1:3)){
		res_tmp = data.frame(model = rep("LowExpression", 4), method = rep(METHODs[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = numeric(4), nclustsd = rep(0, 4), ari = numeric(4), arisd = rep(0, 4))
		if (i ==2){
			df = t(read.table(paste("real_data/low_", methods[i], "_", data, "_candidate0.05.txt", sep=""), header = F, sep = ","))
		}
		else{
			df = t(read.table(paste("real_data/low_", methods[i], "_", data, ".txt", sep=""), header = F, sep = ","))
		}
		if (data == "Usoskin"){
			tags = t(round(2^df))
		}
		else{
			tags = t(round(10^df))
		}
		index = which(apply(df, 1, sum)==0)
		df[index, ] = rnorm(length(index) * dim(df)[2], 0.1, 0.1)
		sData <- CIDR(tags)
		K = sData@nCluster
		res_tmp$nclust = K
		res_tmp$ari[1] = adjustedRandIndex(true_labs, sData@clusters)
		res_tmp$ari[2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
		res_tmp$ari[3] = adjustedRandIndex(true_labs, tSNE(df, K))
		res_tmp$ari[4] = adjustedRandIndex(true_labs, SparsePca(df, K))
		lowexpression.df <- rbind(lowexpression.df, res_tmp)
	}
	all_df = rbind(fulldata.df, dropmodel.df, downsampling.df, lowexpression.df)
	write.table(all_df, file = paste("real_data/", data, "_res.txt", sep=""), quote = F, col.names = T, row.names = F)	
}

########################################## plot result ##############################


############################ sensitivity analysis ####################################
QUANRILES = c("cutoff_0", "cutoff_0.02", "cutoff_0.05", "cutoff_0.1", "cutoff_0.2")
datasets <-c("Ting","Pollen","Deng")
for (data in datasets){
	path <- paste("RealData/",data,"_RAW.txt",sep="")
	RAW <- read.table(path,header=F)
	true_labs <- RAW[1,]
	K = max(as.integer(true_labs))
	fulldata.df <- data.frame()
	for (i in c(1:5)){
		res_tmp = data.frame(model = rep("FullData", 4), method = rep(QUANRILES[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = numeric(4), nclustsd = rep(0, 4), ari = numeric(4), arisd = rep(0, 4))
		df = t(read.table(paste("real_data/rmd_", data, "_candidate", quantiles[i], ".txt", sep=""), header = F, sep = ","))
		if (data == "Usoskin"){
			tags = t(round(2^df-1))
		}
		else{
			tags = t(round(10^df-1))
		}
		sData <- CIDR(tags)
		#K = sData@nCluster
		res_tmp$nclust = K
		res_tmp$ari[1] = adjustedRandIndex(true_labs, sData@clusters)
		res_tmp$ari[2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
		res_tmp$ari[3] = adjustedRandIndex(true_labs, tSNE(df, K))
		res_tmp$ari[4] = adjustedRandIndex(true_labs, SparsePca(df, K))
		fulldata.df <- rbind(fulldata.df, res_tmp)
	}

	dropmodel.df <- data.frame()
	for (i in c(1:5)){
		ari <- matrix(0, repeats, 4)
		nClust <- matrix(0, repeats, 4)
		for (j in c(1:repeats)){
			df = t(read.table(paste("real_data/drop_rmd_", data, "_rep", j, "_candidate", quantiles[i], ".txt", sep=""), header = F, sep = ","))
			if (data == "Usoskin"){
				tags = t(round(2^df-1))
			}
			else{
				tags = t(round(10^df-1))
			}
			sData <- CIDR(tags)
			#K = sData@nCluster
			nClust[j,] = K
			ari[j,1] = adjustedRandIndex(true_labs, sData@clusters)
			ari[j,2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
			ari[j,3] = adjustedRandIndex(true_labs, tSNE(df, K))
			ari[j,4] = adjustedRandIndex(true_labs, SparsePca(df, K))	
		}
		res_tmp = data.frame(model = rep("DropModel", 4), method = rep(QUANRILES[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = apply(nClust, 2, mean), nclustsd = apply(nClust, 2, sd), ari = apply(ari, 2, mean), arisd = apply(ari, 2, sd))
		dropmodel.df <- rbind(dropmodel.df, res_tmp)
	}

	downsampling.df <- data.frame()
	for (i in c(1:5)){
		ari <- matrix(0, repeats, 4)
		nClust <- matrix(0, repeats, 4)
		for (j in c(1:repeats)){
			df = t(read.table(paste("real_data/down_rmd_", data, "_rep", j, "_candidate", quantiles[i], ".txt", sep=""), header = F, sep = ","))
			tags = t(round(10^df-1))
			sData <- CIDR(tags)
			#K = sData@nCluster
			nClust[j,] = K
			ari[j,1] = adjustedRandIndex(true_labs, sData@clusters)
			ari[j,2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
			ari[j,3] = adjustedRandIndex(true_labs, tSNE(df, K))
			ari[j,4] = adjustedRandIndex(true_labs, SparsePca(df, K))	
		}
		res_tmp = data.frame(model = rep("DownSampling", 4), method = rep(QUANRILES[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = apply(nClust, 2, mean), nclustsd = apply(nClust, 2, sd), ari = apply(ari, 2, mean), arisd = apply(ari, 2, sd))
		downsampling.df <- rbind(downsampling.df, res_tmp)
	}

	lowexpression.df <- data.frame()
	for (i in c(1:5)){
		res_tmp = data.frame(model = rep("LowExpression", 4), method = rep(QUANRILES[i], 4), clusm = c("CIDR", "SIMLR", "tSNE", "PCA"),
			nclust = numeric(4), nclustsd = rep(0, 4), ari = numeric(4), arisd = rep(0, 4))
		df = t(read.table(paste("real_data/low_rmd_", data, "_candidate", quantiles[i], ".txt", sep=""), header = F, sep = ","))
		if (data == "Usoskin"){
			tags = t(round(2^df))
		}
		else{
			tags = t(round(10^df))
		}
		index = which(apply(df, 1, sum)==0)
		df[index, ] = rnorm(length(index) * dim(df)[2], 0.1, 0.1)
		sData <- CIDR(tags)
		#K = sData@nCluster
		res_tmp$nclust = K
		res_tmp$ari[1] = adjustedRandIndex(true_labs, sData@clusters)
		res_tmp$ari[2] = adjustedRandIndex(true_labs, SIMLR(t(df), K, cores.ratio = 0.2)$y$cluster)
		res_tmp$ari[3] = adjustedRandIndex(true_labs, tSNE(df, K))
		res_tmp$ari[4] = adjustedRandIndex(true_labs, SparsePca(df, K))
		lowexpression.df <- rbind(lowexpression.df, res_tmp)
	}
	all_df = rbind(fulldata.df, dropmodel.df, downsampling.df, lowexpression.df)
	write.table(all_df, file = paste("real_data/", data, "_sensitivity.txt", sep=""), quote = F, col.names = T, row.names = F)	
}

######
library(ggplot2)
colors = c("#ff7f00","#377eb8","#e41a1c","#4daf4a","#984ea3")
datasets <-c("Usoskin","Ting","Pollen","Deng")

for (data in datasets){
	df <- read.table(paste(data, "_sensitivity.txt", sep=""), header = T)
	df$xx = paste(df$model, df$clusm)
	pdf(paste(data,"_sensitivity.pdf",sep=""),height=5,width=16)
		p<-ggplot(data=df, aes(x=xx, y=ari, fill=method)) +
			geom_bar(stat="identity", position=position_dodge())+ 
			scale_fill_manual(values=colors)+
			ggtitle(data)+theme(plot.title = element_text(hjust = 0.5))+
  			theme_minimal()+geom_errorbar(aes(ymin=ari-arisd, ymax=ari+arisd), width=.2,position=position_dodge(.9)) 		
  		print(p)

	dev.off()
}


