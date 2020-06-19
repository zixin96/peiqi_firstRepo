dataset_peiqi = read.table("spellman.txt", header = T, row.names = 1)
dat = dataset_peiqi[, 23:46]

dat = as.data.frame(dat)


gene = "YAL002W"
gene.numeric = as.numeric(dat[gene,])
gene.numeric[is.na(gene.numeric)] = mean(gene.numeric, na.rm = T)
gene.numeric

# profile plot
plot(c(1,ncol(dat)),range(gene.numeric),type='n',main="Profile plot of 8 random genes",xlab="Samples",ylab="Expression",axes=F)
axis(side=1,at=c(1:24),labels=dimnames(dat)[[2]],cex.axis=0.4,las=2)
axis(side=2)
lines(c(1:ncol(dat)),gene.numeric,lwd=2)

