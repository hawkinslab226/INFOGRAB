##
# setwd("C:/Users/aolim/Desktop/env_aol_dataAnalysis")

##Count Normlization Gene Expression
##RPKM norm
###

##libraries
library(limma)
library(edgeR)

##input files 
gene=read.table("gene2.tab", header=T)
head(gene)
len=read.table("gene_len.tab")
colnames(len)=c("gene", "start" , "end", "size")
head(len)
cnt=read.table("all_cnt_0625.txt", header=T)  
cnt[1:5,1:5]
cnt$gene=sub("^gene-", "", cnt$gene)
cnt$gene=sub("^rna-", "", cnt$gene)
#
cnt[1:5,1:5]
##

##subset the data
len.cnt=as.data.frame(len[which(len$gene %in%
                                  cnt$gene),])

##Gene expression matrix

##filter the genes with size
rm = (cnt$gene %in% len.cnt$gene)
table(rm)
data =cnt[rm,]
dim(data)
data[1:5,1:5]

#transform rownames
row.names(data) = data[,1] 
data = data[,-1]          
str(data)
data[1:5,1:5]

##order
data[1:10,1:2]
data=data[order(rownames(data)),]
data[1:10,1:2]
#
len.cnt[1:5,1:2]
len.cnt=len.cnt[order(len.cnt$gene),]



##remove the zero
## total counts per gene
Totalcounts = rowSums(data)
## genes with zero count?
table(Totalcounts==0)
## filter genes with 0 counts
rm = rowMeans(data)==0
data.expr = data [!rm,]
dim(data.expr)#22647    40

#filter by cnt.expr
len.cnt1=as.data.frame(len.cnt[which(len.cnt$gene %in%
                                  rownames(data.expr)),])
dim(len.cnt1)

#check the order
data.expr[1:5,1:5]
x=c(len.cnt1$gene)
y=c(rownames(data.expr))
order=as.data.frame(cbind(x,y))


##normalization data in RPKM
y = DGEList(counts=data.expr,genes=data.frame(Length=len.cnt1$size))
y = calcNormFactors(y)
RPKM = rpkm(y)
head(RPKM)
RPKM[1:5,1:5]

log.rpkm=log(RPKM+1, base=2)
log.rpkm[1:5,1:5]

##saving the data
write.table(len.cnt1, "len0625_filter.txt",sep = "\t", row.names = F, quote = F)
write.table(log.rpkm, "logRPKM0625_filter.txt",sep = "\t", row.names = T, quote = F)