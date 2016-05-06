## install required packages (only at first time)
install.packages(c("tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit","vegan"))

## load class definition and functions
source("RaceID2_StemID_class.R")

## input data
x <- read.csv("transcript_counts_intestine_5days_YFP.xls",sep="\t",header=TRUE)
rownames(x) <- x$GENEID
# prdata: data.frame with transcript counts for all genes (rows) in all cells (columns); with rownames == gene ids; remove ERCC spike-ins 
prdata <- x[grep("ERCC",rownames(x),invert=TRUE),-1]

## RaceID2
# initialize SCseq object with transcript counts
sc <- SCseq(prdata)
# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=500, downsample=TRUE, dsn=1, rseed=17000)
# k-medoids clustering
sc <- clustexp(sc,clustnr=30,bootnr=50,metric="pearson",do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids")
# compute t-SNE map
sc <- comptsne(sc,rseed=15555)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.95)

## diagnostic plots
# gap statistics: only if do.gap == TRUE
##plotgap(sc)
# plot within-cluster dispersion as a function of the cluster number: only if sat == TRUE
plotsaturation(sc,disp=TRUE)
# plot change of the within-cluster dispersion as a function of the cluster number: only if sat == TRUE
plotsaturation(sc)
# silhouette of k-medoids clusters
plotsilhouette(sc)
# Jaccard's similarity of k-medoids clusters
plotjaccard(sc)
# barchart of outlier probabilities
plotoutlierprobs(sc)
# regression of background model
plotbackground(sc)
# dependence of outlier number on probability threshold (probthr)
plotsensitivity(sc)
# heatmap of k-medoids cluster
clustheatmap(sc,final=FALSE,hmethod="single")
# heatmap of final cluster
clustheatmap(sc,final=TRUE,hmethod="single")
# highlight k-medoids clusters in t-SNE map
plottsne(sc,final=FALSE)
# highlight final clusters in t-SNE map
plottsne(sc,final=TRUE)
# highlight cell labels in t-SNE map
plotlabelstsne(sc,labels=sub("(\\_\\d+)","",names(sc@ndata)))
# highlight groups of cells by symbols in t-SNE map
plotsymbolstsne(sc,types=sub("(\\_\\d+)$","", names(sc@ndata)))
# highlight transcirpt counts of a set of genes in t-SNE map, e. g. all Apoa genes
g <- c("Apoa1__chr9", "Apoa1bp__chr3", "Apoa2__chr1", "Apoa4__chr9", "Apoa5__chr9")
plotexptsne(sc,g,n="Apoa genes",logsc=TRUE)

## identification of marker genes
# differentially regulated genes in each cluster compared to the full ensemble
cdiff <- clustdiffgenes(sc,pvalue=.01)

## write results to text files
# final clusters 
x <- data.frame(CELLID=names(sc@cpart),cluster=sc@cpart)
write.table(x[order(x$cluster,decreasing=FALSE),],"cell_clust.xls",row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
# differentially expressed genes in cluster
for ( n in names(cdiff) ) write.table(data.frame(GENEID=rownames(cdiff[[n]]),cdiff[[n]]),paste(paste("cell_clust_diff_genes",sub("\\.","\\_",n),sep="_"),".xls",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

# differentially expressed genes between two sets of clusters, e. g. cluster 1 and clusters 2,3
d <- diffgenes(sc,cl1=1,cl2=c(2,3),mincount=5)
plotdiffgenes(d,gene=names(d$z)[1])


## StemID

# initialization
ltr <- Ltree(sc)
# computation of the entropy
ltr <- compentropy(ltr)
# computation of the projections for all cells
ltr <- projcells(ltr,cthr=2,nmode=FALSE)
# computation of the projections for all cells after randomization
ltr <- projback(ltr,pdishuf=2000,nmode=FALSE,rseed=17000)
# assembly of the lineage tree
ltr <- lineagetree(ltr,pthr=0.01,nmode=FALSE)
# determination of significant differentiation trajectories
ltr <- comppvalue(ltr,pethr=0.01,nmode=FALSE)

## diagnostic plots
# histogram of ratio between cell-to-cell distances in the embedded and the input space
plotdistanceratio(ltr)
# t-SNE map of the clusters with more than cthr cells including a minimum spanning tree for the cluster medoids
plotmap(ltr)
# visualization of the projections in t-SNE space overlayed with a minimum spanning tree connecting the cluster medoids
plotmapprojections(ltr)
# lineage tree showing the projections of all cells in t-SNE space
plottree(ltr,showCells=TRUE,nmode=FALSE,scthr=0)
# lineage tree without showing the projections of all cells
plottree(ltr,showCells=FALSE,nmode=FALSE,scthr=0)
# heatmap of the enrichment p-values for all inter-cluster links
plotlinkpv(ltr)
# heatmap of the link score for all inter-cluster links
plotlinkscore(ltr)
# heatmap showing the fold enrichment (or depletion) for significantly enriched or depleted links
projenrichment(ltr)

## extract projections onto all links for all cells in a given cluster i
x <- getproj(ltr,i=1)
# heatmap of all projections for cluster i
pheatmap(x$pr)
# heatmap of z-score for all projections for cluster i
pheatmap(x$prz)

## extracting all cells on two branches sharing the same cluster and computing differentially expressed genes between these two branches
x <- branchcells(ltr,list("1.3","1.2"))
# z-scores for differentially expressed genes
head(x$diffgenes$z)
# plotting the cells on the two branches as additional clusters in the t-SNE map
plottsne(x$scl)


## computing the StemID score
x <- compscore(ltr,nn=1)
#plotting the StemID score
plotscore(ltr,1)

 
