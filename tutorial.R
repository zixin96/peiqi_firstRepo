# scatter plot matrix
dat <- read.table("spellman.txt",header=T,row.names=1)
dat <- as.data.frame(dat);

# other data sets in R to use
library(Biobase);	library(annotate);	library(golubEsets);
data(golubTrain);	 data(golubTest); 	data(geneData);
dat <- geneData or dat <- exprs(golubTrain) or dat <- exprs(golubTest)

# Boxplot of MPG by Car Cylinders
boxplot(mpg~cyl,data=mtcars, main="Car Milage Data",
        xlab="Number of Cylinders", ylab="Miles Per Gallon")



# box plots
boxplot(log2(dat),cex=0.3,col="lightblue",las=2,cex.axis=0.6,ylab="log2(intensity)",main="Alon Colon Cancer dataset - all 2K probesets")

# random selection of 5 samples
rand.sams <- sample(names(dat),5,replace=F)
# plot trellis
pairs(log2(dat[,rand.sams]),pch=21,col=1,bg='salmon',main="Alon colon cancer dataset\n5 random patient samples",cex=0.4)

# Pearson's correlation matrix
library(gplots)
dat.cor <- cor(dat)

layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(dat.cor,na.rm=T),max(dat.cor,na.rm=T),length=10)
image(dat.cor,main="Correlation plot Normal/Tumor data",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(dat.cor)),label=dimnames(dat.cor)[[2]],cex.axis=0.9,las=2)

image(as.matrix(leg),col=cx,axes=F)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)

# random sample of 8 genes
rand.genes <- sample(dimnames(dat)[[1]],8,replace=F)

# center the gene values to be on similar scales
dat<-t(scale(t(dat)))




# profile plot
plot(c(1,ncol(dat)),range(dat[rand.genes,]),type='n',main="Profile plot of 8 random genes",xlab="Samples",ylab="Expression",axes=F)
axis(side=1,at=c(1:62),labels=dimnames(dat)[[2]],cex.axis=0.4,las=2)
axis(side=2)
for(i in 1:length(rand.genes)) {
  dat.y <- as.numeric(dat[rand.genes[i],])
  lines(c(1:ncol(dat)),dat.y,col=i,lwd=2)
}

# pca biplot
biplot(prcomp(t(log2(dat[1:10,]))),cex=0.6,main="PCA biplot - genes and patient samples",col=c("black","grey"),expand=0.8)

# k-means cluster profiles
d.k <- kmeans(log2(dat),6)
par(mfrow=c(2,3))
for(i in 1:6) {
  tmp <- dat[d.k$cluster==i,]
  matplot(c(1:ncol(dat)),log2(t(tmp)),type='l',col="lightgrey",xlab='',ylab='log2(intensity)',axes=F)
  me <- as.numeric(apply(log2(tmp),2,mean))
  lines(c(1:ncol(dat)),me,lwd=2,col='red')
  axis(1,at=c(1:ncol(dat)),dimnames(dat)[[2]],las=2,cex.axis=0.6)
  axis(2)
  title(main=paste('Cluster',i))
}

# cv vs. mean plot
dat.mean <- apply(log2(dat),2,mean)		# calculate mean for each sample
dat.sd <- sqrt(apply(log2(dat),2,var))		# calculate st.deviation for each sample
dat.cv <- dat.sd/dat.mean			#calculate cv

plot(dat.mean,dat.cv,main="Alon colon cancer dataset\nSample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.5)

# average correlation plot
dat.avg <- apply(dat.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of Tumor/Normal samples",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

# 2D sample pca plot
dat.pca <- prcomp(t(dat))
dat.loads <- dat.pca$x[,1:2]
plot(dat.loads[,1],dat.loads[,2],main="Sample PCA plot",xlab="p1",ylab="p2",col='red',cex=1,pch=15)
text(dat.loads,label=dimnames(dat)[[2]],pos=1,cex=0.5)

# k-means clustering for missing value imputation
dat <- dat[2:30,]				# only use 29 genes for example
cl <- kmeans(dat[,-1],centers=5, iter.max=20)		# cluster into 5 groups
# we pretend to be missing a value at sample#1 gene #2
groups <- cl$cluster			# get cluster membership for each gene
groups				# look at groups to see where gene 2 is
group.2 <- groups==2			# since gene 2 is in group 2, get all other members
genes.cluster <- dimnames(dat)[[1]][group.2]
genes.cluster			# look at all other genes in cluster #2

gene.dist <- dist(dat[genes.cluster,-1],method="euclidean")	# get distances from genes in cluster 2 to 					# gene #2
gene.dist <- as.matrix(gene.dist)
gene.dist <- gene.dist[2:5,1]
gene.weight <- as.numeric(gene.dist/sum(gene.dist))	# get weights for each gene

weight.mean <-  weighted.mean(dat[genes.cluster[-1],1], gene.weight)	# calculate weighted mean for 						# gene #2

# perspective plot
data(volcano)	# load volcano data set
persp(volcano, theta=45, phi=30, col="red")


# parameter histogram
dat.25 <- apply(dat,2,quantile,probs=0.25)
dat.75 <- apply(dat,2,quantile,probs=0.75)
dat.iqr <- dat.75-dat.25 
l1 <- median(dat.iqr,na.rm=TRUE) - 2*mad(dat.iqr,na.rm=TRUE)
l2 <- median(dat.iqr,na.rm=TRUE) + 2*mad(dat.iqr,na.rm=TRUE)
hist(dat.iqr,main="Distribution of IQRs for Tumor Samples",cex.main=1.5,col="salmon")
abline(v=c(l1,l2),lty=2,col="red")

# calculate mean for some genes, with respect to class
library(multtest)
data(golub)
dat <- as.data.frame(golub)
ann <- golub.cl
dat.aml <- apply(dat[,ann==1],1,mean)
dat.all <- apply(dat[,ann==0],1,mean)
tab <- data.frame(rbind(dat.aml[1:20],dat.all[1:20]))
dimnames(tab)[[1]] <- c("AML","ALL")
names(tab) <- dimnames(dat)[[1]][1:20]
mp <- barplot(tab)
tot <- colMeans(tab)
text(mp, tot + 3, format(tot), xpd = TRUE, col = "blue")
barplot(as.matrix(tab),beside=T,col=c("red","yellow"),legend=rownames(as.matrix(tab)),ylim=c(-5,5),ylab="Expression")
title(main = "Mean Expression Levels of first 20 genes")

# cluster tree
dat <- t(dat)			#transpose dat
dat.dist <- dist(dat,method="euclidean")	# calculate distance
dat.clust <- hclust(dat.dist,method="single")	# calculate clusters
plot(dat.clust,labels=names(dat),cex=0.75)	# plot cluster tree
