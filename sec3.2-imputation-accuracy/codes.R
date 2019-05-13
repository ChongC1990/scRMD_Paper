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


########################################################### Zeigenhain #########################################################################

Methods <- c("CELseq","DropSeq","MARSseq","SCRBseq","SmartSeq2","SmartSeq")
threshold = 0.05
Zeig <- read.table("RealData/Zeig_complete.txt",header=T)
ERCC <- read.table("RealData/ERCC.txt",header=F)
for (i in 1:6){
	if (i==6){
		ind = setdiff(grep(Methods[i],colnames(Zeig)),grep(Methods[i-1],colnames(Zeig)))
	}
	else{	
		ind <- grep(Methods[i],colnames(Zeig))
	}
	
	correlation <- data.frame(RAW=numeric(length(ind)),scRMD=numeric(length(ind)),scImpute=numeric(length(ind)),MAGIC=numeric(length(ind)))
	data <- Zeig[,ind]
	data <- data[apply(data!=0,1,sum)/length(ind)>threshold,]
	data <- as.matrix(data)%*%diag(1/apply(data,2,sum))
	data <- data*1000000
	logdata <- log10(data+1)
	res.rmd <- rmd(as.matrix(logdata))
	plist = get_mix_parameters(logdata + log10(1.01),ncore=20); 
	res.sci = imputation_model1(count = logdata+log10(1.01), point = log10(1.01), plist, drop_thre = 0.5,ncores=20);
	res.sci = res.sci - log10(1.01)
	res.mag <- t(magic(t(as.matrix(logdata))))

	ercc_index = grep("ERCC",rownames(logdata))
	index <- GetIndex(rownames(logdata)[ercc_index],as.character(ERCC[,2]))
	for (j in 1:length(ind)){
		correlation[j,1]=cor(logdata[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,2]=cor(res.rmd$exprs[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,3]=cor(res.sci[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,4]=cor(res.mag[ercc_index,j],log10(ERCC[index,4]))
	}
	corre <- melt(correlation)
	colnames(corre)=c("method","correlation")
	pdf(paste("barplot",Methods[i],".pdf",sep=""),height=4,width=4)
	p <- ggplot(corre, aes(x=method, y=correlation,fill=method))+ 
  	geom_boxplot(outlier.colour="red", outlier.shape=8,
                outlier.size=4) +scale_fill_manual(values=c("#377eb8","#e41a1c","#4daf4a","#984ea3"))
	print(p)
	dev.off()

}


##### Wilcoxon signed-rank test ####

Wilcoxon <- data.frame(RAW.scRMD=numeric(5),RAW.scImpute=numeric(5),
	RAW.MAGIC=numeric(5),scRMD.scImpute=numeric(5),
	scRMD.MAGIC=numeric(5),scImpute.MAGIC=numeric(5))

Methods <- c("CELseq","MARSseq","SCRBseq","SmartSeq2","SmartSeq")
threshold = 0.05
Zeig <- read.table("RealData/Zeig_complete.txt",header=T)
ERCC <- read.table("RealData/ERCC.txt",header=F)
for (i in 1:5){
	if (i==5){
		ind = setdiff(grep(Methods[i],colnames(Zeig)),grep(Methods[i-1],colnames(Zeig)))
	}
	if (i<5){	
		ind <- grep(Methods[i],colnames(Zeig))
	}
	
	correlation <- data.frame(RAW=numeric(length(ind)),scRMD=numeric(length(ind)),scImpute=numeric(length(ind)),MAGIC=numeric(length(ind)))
	data <- Zeig[,ind]
	data <- data[apply(data!=0,1,sum)/length(ind)>threshold,]
	data <- as.matrix(data)%*%diag(1/apply(data,2,sum))
	data <- data*1000000
	logdata <- log10(data+1)
	res.rmd <- rmd(as.matrix(logdata))
	plist = get_mix_parameters(logdata + log10(1.01),ncore=20); 
	res.sci = imputation_model1(count = logdata+log10(1.01), point = log10(1.01), plist, drop_thre = 0.5,ncores=20);
	res.sci = res.sci - log10(1.01)
	res.mag <- t(magic(t(as.matrix(logdata))))

	ercc_index = grep("ERCC",rownames(logdata))
	index <- GetIndex(rownames(logdata)[ercc_index],as.character(ERCC[,2]))
	for (j in 1:length(ind)){
		correlation[j,1]=cor(logdata[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,2]=cor(res.rmd$exprs[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,3]=cor(res.sci[ercc_index,j],log10(ERCC[index,4]))
		correlation[j,4]=cor(res.mag[ercc_index,j],log10(ERCC[index,4]))
	}
	Wilcoxon[i,1] = wilcox.test(correlation$RAW,correlation$scRMD, paired=FALSE)$p.value
	Wilcoxon[i,2] = wilcox.test(correlation$RAW,correlation$scImpute, paired=FALSE)$p.value
	Wilcoxon[i,3] = wilcox.test(correlation$RAW,correlation$MAGIC, paired=FALSE)$p.value
	Wilcoxon[i,4] = wilcox.test(correlation$scRMD,correlation$scImpute, paired=FALSE)$p.value
	Wilcoxon[i,5] = wilcox.test(correlation$scRMD,correlation$MAGIC, paired=FALSE)$p.value
	Wilcoxon[i,6] = wilcox.test(correlation$scImpute,correlation$MAGIC, paired=FALSE)$p.value
}

write.table(t(Wilcoxon),file="Wilcoxon.txt",quote=F,col.names=F)
############ Mean Error ################

Methods <- c("CELseq","DropSeq","MARSseq","SCRBseq","SmartSeq2","SmartSeq")

threshold = 0.05
Zeig <- read.table("RealData/Zeig_complete.txt",header=T)
ERCC <- read.table("RealData/ERCC.txt",header=F)
MeanError <- data.frame(scRMD=numeric(6),scImpute=numeric(6),MAGIC=numeric(6))
for (i in 1:6){
	if (i==6){
		ind = setdiff(grep(Methods[i],colnames(Zeig)),grep(Methods[i-1],colnames(Zeig)))
	}
	if(i!=6){	
		ind <- grep(Methods[i],colnames(Zeig))
	}
	data <- Zeig[,ind]
	data <- data[apply(data!=0,1,sum)/length(ind)>threshold,]
	data <- as.matrix(data)%*%diag(1/apply(data,2,sum))
	data <- data*1000000
	logdata <- log10(data+1)
	Pro <- exp(-0.2*logdata*logdata)
	N = dim(logdata)[1]; P = dim(logdata)[2]
	set.seed(2017)
	raw <- logdata * matrix(rbinom(N*P, 1, 1 - Pro), N, P)
	drop <- which(raw==0)
	res.rmd <- rmd(as.matrix(raw))
	plist = get_mix_parameters(raw + log10(1.01),ncore=20); 
	res.sci = imputation_model1(count = raw+log10(1.01), point = log10(1.01), plist, drop_thre = 0.5,ncores=20);
	res.sci = res.sci - log10(1.01)
	res.mag <- t(magic(t(as.matrix(raw))))
	MeanError[i,1] = mean(abs((res.rmd$exprs-logdata)[drop]))
	MeanError[i,2] = mean(abs((res.sci-logdata)[drop]))
	MeanError[i,3] = mean(abs((res.mag-logdata)[drop]))
}


######### Cell Correlation ################
Methods <- c("CELseq","DropSeq","MARSseq","SCRBseq","SmartSeq2","SmartSeq")
alpha = 0.1
ind <- grep(Methods[1],colnames(Zeig))
data <- Zeig[,ind]
N = dim(data)[2]
RowSum <- apply(data>0,1,sum)
index = which(RowSum/N>=alpha)

for (i in 2:6){
	if (i==6){
		ind = setdiff(grep(Methods[i],colnames(Zeig)),grep(Methods[i-1],colnames(Zeig)))
	}
	if(i!=6){	
		ind <- grep(Methods[i],colnames(Zeig))
	}
	data <- Zeig[,ind]
	N = dim(data)[2]
	RowSum <- apply(data>0,1,sum)
	index <- intersect(index,which(RowSum/N>=alpha))
}

HighZeig <- Zeig[index,]
ImputedData <- list(raw=list(),rmd=list(),sci=list(),mag=list())

for (i in 1:6){
	if (i==6){
		ind = setdiff(grep(Methods[i],colnames(HighZeig)),grep(Methods[i-1],colnames(HighZeig)))
	}
	if(i!=6){	
		ind <- grep(Methods[i],colnames(HighZeig))
	}
	data <- HighZeig[,ind]
	#data <- data[apply(data!=0,1,sum)/length(ind)>threshold,]
	data <- as.matrix(data)%*%diag(1/apply(data,2,sum))
	data <- data*1000000
	raw <- log10(data+1)
	ImputedData$raw[[i]] <- raw
	ImputedData$rmd[[i]] <- rmd(as.matrix(raw))$exprs
	plist = get_mix_parameters(raw + log10(1.01),ncore=20); 
	res.sci = imputation_model1(count = raw+log10(1.01), point = log10(1.01), plist, drop_thre = 0.5,ncores=20);
	ImputedData$sci[[i]] = res.sci - log10(1.01)
	ImputedData$mag[[i]] <- t(magic(t(as.matrix(raw))))
}

correlation <- data.frame(RAW=numeric(length(15)),scRMD=numeric(length(15)),scImpute=numeric(length(15)),MAGIC=numeric(length(15))) 

k=1
for (i in 1:5){
	for (j in (i+1):6){
		correlation[k,1] = median(cor(ImputedData$raw[[i]],ImputedData$raw[[j]]))
		correlation[k,2] = median(cor(ImputedData$rmd[[i]],ImputedData$rmd[[j]]))
		correlation[k,3] = median(cor(ImputedData$sci[[i]],ImputedData$sci[[j]]))
		correlation[k,4] = median(cor(ImputedData$mag[[i]],ImputedData$mag[[j]]))
		k=k+1
	}
}


corre <- melt(correlation)
colnames(corre)=c("method","correlation")
pdf("barplot.pdf",height=4,width=3)
p <- ggplot(corre, aes(x=method, y=correlation)) + 
geom_boxplot(outlier.colour="red", outlier.shape=8,
    outlier.size=4)+geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)
p
dev.off()
