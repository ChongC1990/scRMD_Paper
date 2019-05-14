####################################################################### 
#######    Imputation for Ultra-high-throughput scRNA-seq data
####### Requirement: conda 4.5.4, R 3.4.3, Rtsne 0.13, Rcpp 0.12.18, 
#######              Ramgic 1.4.0, scImpute 0.0.9
#######################################################################
source('function.R')

### Karaiskos
ImputationAndDimReduction(Raw_fileName='Karaiskos.txt',Cluster_fileName='Karaiskos_clust.txt',
	file='Karaiskos',Kcluster=4,cores=20,seed.value=12345)
plotFig(file='Karaiskos',Cluster_fileName='Karaiskos_clust.txt')

### PBMC
ImputationAndDimReduction(Raw_fileName='PBMC.txt',Cluster_fileName='PBMC_clust.txt',
	file='PBMC',Kcluster=1,cores=20,seed.value=12345)
plotFig(file='PBMC',Cluster_fileName='PBMC_clust.txt')

### Hrvatin
ImputationAndDimReduction(Raw_fileName='Hrvatin.txt',Cluster_fileName='Hrvatin_clust.txt',
	file='Hrvatin',,Kcluster=1,cores=20,seed.value=16431)
plotFig(file='Hrvatin',Cluster_fileName='Hrvatin_clust.txt')

### Alles
ImputationAndDimReduction(Raw_fileName='Alles.txt',Cluster_fileName='Alles_clust.txt',
	file='Alles',Kcluster=1,cores=20,seed.value=3454)
plotFig2(file='Alles',Cluster_fileName='Alles_clust.txt')
