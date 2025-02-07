setwd("/net/hawkins/vol1/home/aolima/data_analysisR/gene_express")
#Genome annoation Figures
##FAANG
##Andressa O. de Lima
##******************************************************************************
##colors  vector
mycolors1 <- c(
  `Kidney`="#43009A",
  `Trachea`="#990099",
  `BCell`="#0DFAFA",
  `Spleen_T_Cells`="#13B7B7",
  `Bursa`="#004C99",
  `Thymus`="#2685E4",
  `d0_macrophage`="#F02B6D",
  `d3_macrophage`="#FF0091",
  `d6_macrophage`="#F572BC",
  `Macrophage.lung`="#D06AAA",
  `Monocyte.blood,`="#91155B",
  `Ileum`="#98D55A",
  `Jejunum`="#4C9900",
  `Proximal.Cecum`="#CCFF99",
  `Dark.meat`="#A1122A",
  `White.meat`="#FFCCCC",
  `Isthmus.od2`="#FF8000",
  `Magnum.od1`="#DE5100",
  `Shell.Gland.od3`="#FFC78E",
  `Ovary`="#CCCC00",
  `NCBI` ="#5C4B99",
  `Ensembl` ="#FFDBC3")
#function colors
mycols1 <- function(...) {
  cols <- c(...)
  if (is.null(cols))
    return (mycolors1)
  mycolors1[cols]
}

###To request the color tissue use the names
color_names= names(mycolors1)
color_names
###*****************************************************************************
###*RHodor Colors
###*https://github.com/AndressaOL/RHodor_colors
library(devtools)
library(scales)
source("/net/hawkins/vol1/home/aolima/tools/RHodor_colors/R/RHodor_palettes.R")
colors1=hodor_pal("hawkins_s")(15)
colors1
show_col(colors1)

##end**************************************************************************
##library colors
library(RColorBrewer)
library(limma)
library(edgeR)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(dendextend)
library(Biobase)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(uwot)

##output
out=file.path("output",fsep="/")

##******************************************************************************
##*
#fpkm
fpkm=read.table("all_transcript_fpkm_raw0718.txt", header=T, row.names = 1)
fpkm[1:5,1:5]

#cnt
cnt=read.delim("all_cnt_0718.txt", header=T, row.names = 1)
cnt[1:5,1:5]
##
##pheno
pheno=read.csv("Pheno_data071523.csv", header=T)
pheno[1:4,1:4]
##
##domain
domain=read.table("top5.domain_code07025.txt", header=T)
predict=read.table("merge.predictedlncRNA.txt", header=T)

##filter iso with top5 domain
top=as.data.frame(unique(domain$name))
colnames(top)= "transcript"
write.table(top, file= "top_domain_iso0725.txt", row.names = F, sep = "\t",
            na= "NA", quote = F)

##gene expression cluster

##cnt
##filter zero
#remove genes non-expressed and low expressed
## total counts per gene
Totalcounts = rowSums(cnt)
## genes with zero count?
table(Totalcounts==0)
## filter genes with 0 counts
rm = rowMeans(cnt)==0
cnt.filter = cnt[!rm,]
dim(cnt.filter) #22648    40

###normalization using TMM
d = DGEList(cnt.filter)
TMM = calcNormFactors(d, method = "TMM")
head(TMM$samples)
#get the normalized counts:
cpm = cpm(TMM, log=FALSE)
cpm=log(cpm+1, base=2)  #log scale

##FPKM
## total counts per gene
Totalcounts = rowSums(fpkm)
## genes with zero count?
table(Totalcounts==0)
## filter genes with 0 counts
rm = rowMeans(fpkm)==0
##
fpkm.filter = fpkm[!rm,]
dim(fpkm.filter) #95127    40

#dendogram
datExpr = t(as.data.frame(cpm))
sampleTree = hclust(dist(datExpr), method = "average")
lab=labels(sampleTree)

##convert
dend <- as.dendrogram(sampleTree)

##plot theme
#select theme
theme01= theme(
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  axis.title.y = element_text(size = 14),
  axis.text.y = element_text(size = 12),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

png(file = "gene_expressed_cluster.png", width = 16, height = 12,
    units = "in", res = 600)
par(mar = c(10, 10, 4, 2)) 
##dendogram
dend %>%
  set("labels_col", value = c("#2D2727", "#FC6736","#0C2D57"), k=3) %>%
  set("branches_k_color", value = c("#2D2727", "#FC6736","#0C2D57"), k = 3) %>%
  #set("branches_lty", 4) %>%
  plot(horiz=FALSE, axes=FALSE)
##tittles
title(main = "Cluster Dendrogram")
xlabel <- ""
ylabel <- ""
mtext(xlabel, side = 1, line = 8, cex = 1.5)  # Adjust line and size for x-axis label
mtext(ylabel, side = 2, line = 8, cex = 1.5)  # Adjust line and size for y-axis label
dev.off()

#plot2
##tissue colors: mycolors1
datExpr = t(as.data.frame(cpm))
datExpr[1:5,1:5]

head(pheno)
head(datExpr)
#randomly select columns
#set.seed(123)
#random_cols= sample(ncol(datExpr), 2000)
#data_selected = datExpr[, random_cols]
#
##top variable genes
variances= apply(datExpr, 2, var)
top_genes = order(variances, decreasing = TRUE)[1:2000]
data_selected= datExpr[, top_genes]
##
data.t=t(data_selected)
#
sampleTree1 = hclust(dist(data.t), method = "average")
lab=labels(sampleTree1)
#
row_dend = as.dendrogram(sampleTree1)
##
# Subset pheno to match selected samples
pheno.filter=pheno[,c("SampleID.", "Tissue","System")]
pheno_selected= pheno.filter[match(rownames(data_selected), pheno.filter$SampleID.), ]
row_anno= rowAnnotation(System = pheno_selected$System)
#
column_anno= HeatmapAnnotation(
  df = data.frame(System = pheno_selected$System), 
  col = list(System = c("Immune" = "#31E1F7", "Respiratory" = "#400D51", 
                        "Excretory"="#FF7777", "Muscular" = "#D800A6",
                        "Intestine" ="#6499E9",  "Reproductive" = "#836FFF"))
)

# Generate the heatmap
heatmap= Heatmap(data.t, 
                   col = colorRampPalette(wes_palette("Zissou1"))(100),
                   name = "Expression", 
                   cluster_rows = row_dend, 
                   show_row_dend = TRUE, 
                   row_dend_reorder = FALSE, 
                   show_row_names = FALSE,
                   show_column_names = TRUE, 
                   heatmap_legend_param = list(title = "Expression levels"),
                  top_annotation = column_anno)

png(file = "gene_expressed_to2000.png", width = 16, height = 8,
   units = "in", res = 300)
par(mar = c(10, 10, 4, 2)) 
heatmap
dev.off()

##PCA 
cpm[1:5,1:5]
head(pheno.filter)

# Perform UMAP using uwot
datExpr = t(as.data.frame(cpm))
datExpr[1:5,1:5]
umap_result= uwot::umap(datExpr)
umap_data= as.data.frame(umap_result)
#umap_data$ = rownames(umap_data)
umap_data1= merge(umap_data, pheno.filter, by.x = "row.names", by.y = "SampleID.")

#plot
# Plot UMAP
head(pheno.filter)
mycolors1
col.tissue=c("BCell" ="#0DFAFA", "Bursa"="#004C99","Proximal.Cecum"="#CCFF99",
        "Macrophage(D0)" = "#F02B6D" , "Macrophage(D3)" = "#FF0091",
        "Macrophage(D6)" = "#F572BC", "Iliotibial.major"= "#A1122A",
        "Ileum" = "#98D55A", "Jejunum" ="#4C9900", "Kidney" ="#43009A",
        "Macrophage(lung)" = "#43009A", "Monocyte(blood)" = "#91155B", 
        "Magnum" ="#DE5100", "Isthmus"= "#FF8000", "Shell.Gland" = "#FFC78E",
        "Ovary"= "#CCCC00", "TCell(spleen)" = "#13B7B7", "Thymus" ="#2685E4",
        "Trachea"="#990099", "Pectoralis.major" = "#FFCCCC")


umap_plot= ggplot(umap_data1, aes(x = V1, y = V2, color = Tissue, 
                                  shape = System)) +
  geom_point(size = 2) +
  labs(title = " Gene Expression Profile", x = "UMAP1", y = "UMAP2") +
  theme01 +
  scale_color_manual(values = col.tissue) +
  ggthemes::theme_few() 

ggsave(here::here(out,"UMAP_gene.svg"), plot = umap_plot, width = 12, height = 10,
       units = "in", dpi = 300)

##

#explore isoforms
transcript=read.table("all40_genomeAnnoatation_ncbi_copy.tab", header=T)
head(transcript)

transcript_code=as.data.frame(transcript[,c("qry_id","class_code")])
head(transcript_code)
###

head(pheno.filter)
head(predict)
dim(predict)



fpkm.filter[1:5,1:5]

##filter data
fpkm.predict= which(row.names(fpkm) %in% predict$name)
filtered_fpkm=fpkm[fpkm.predict, ]
filtered_fpkm[1:5,1:5]
dim(filtered_fpkm) #6058   40

##expressed
##filter zero
## total counts per gene
Totalcounts = rowSums(filtered_fpkm)
## genes with zero count?
table(Totalcounts==0)
## filter genes with 0 counts
rm = rowMeans(filtered_fpkm)==0
fpkm.expressed = filtered_fpkm[!rm,]
dim(fpkm.expressed) #5680   40

#plot isoforms 
library(pheatmap)

dat.t=t(fpkm.expressed)
# Annotation data frame
annotation= data.frame(Tissues = pheno.filter$Tissue)
rownames(annotation)= colnames(fpkm.expressed)

# Annotation colors
annotation_colors= list(
   Tissues= c("BCell" ="#0DFAFA", "Bursa"="#004C99","Proximal.Cecum"="#CCFF99",
    "Macrophage(D0)" = "#F02B6D" , "Macrophage(D3)" = "#FF0091",
    "Macrophage(D6)" = "#F572BC", "Iliotibial.major"= "#A1122A",
    "Ileum" = "#98D55A", "Jejunum" ="#4C9900", "Kidney" ="#43009A",
    "Macrophage(lung)" = "#43009A", "Monocyte(blood)" = "#91155B", 
    "Magnum" ="#DE5100", "Isthmus"= "#FF8000", "Shell.Gland" = "#FFC78E",
    "Ovary"= "#CCCC00", "TCell(spleen)" = "#13B7B7", "Thymus" ="#2685E4",
    "Trachea"="#990099", "Pectoralis.major" = "#FFCCCC")
  )

# Generate the heatmap
heatmap=pheatmap(
  fpkm.expressed, 
  color = colorRampPalette(wes_palette("Zissou1"))(100),
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  scale="row",
  drop_levels = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  legend_breaks = c(-6, 0, 6),
  legend_labels = c("low", "medium" ,"high"),
  legend = TRUE,
  main = "Predicted Transcript Isoforms Expressed",
  fontsize_col = 14,
  fontsize_row = 14
)

# Save the heatmap to a file
png(file = "isoforms_expressed.png", width = 16, height = 8, units = "in", res = 300)
heatmap
dev.off()

#plots
head(fpkm.expressed)
head(predict)

#boxplot mRNA vcs LncRNA
##filter data
lncRNA.predict=subset(predict, type == "lncRNA")
mRNA.predict=subset(predict, type == "mRNA")

mRNA.expressed= which(row.names(fpkm.expressed) %in% mRNA.predict$name)
mRNA.expressed.filter=fpkm.expressed[mRNA.expressed, ]
mRNA.expressed.filter[1:5,1:5]
dim(mRNA.expressed.filter) #1358   40

lncRNA.expressed= which(row.names(fpkm.expressed) %in% lncRNA.predict$name)
lncRNA.expressed.filter=fpkm.expressed[lncRNA.expressed, ]
lncRNA.expressed.filter[1:5,1:5]
dim(lncRNA.expressed.filter) #4322   40


#vector tissues 
tissue_groups= list(
  Bcell = c("Bcell.01","BCell.02"),
  Bursa = c("Bursa.01","Bursa.02"),
  Proximal_Cecum = c("P.cecum.01","P.cecum.02"),
  D0 = c("Macrophage.D0.01","Macrophage.D0.02"),
  D3 = c("Macrophage.D3.01","Macrophage.D3.02"),
  D6 = c("Macrophage.D6.01","Macrophage.D6.02"),
  Iliotibial = c("T.muscle.01", "T.muscle.02"),
  Ileum = c("Ileum.01","Ileum.02"),
  Jejunum = c("Jejunum.01","Jejunum.02"),
  Kidney = c("Kidney.01","Kidney.02"),
  Macrophage_lung = c("Macrophage.lung.01","Macrophage.lung.02"),
  Monocyte = c("Monocyte.01","Monocyte.02"),
  Magnum = c("Magnum.od1.01","Magnum.od1.01"),
  Isthmus = c("Isthmus.od2.01","Isthmus.od2.02"),
  Shell_Gland = c("Shell.Gland.od3.01","Shell.Gland.od3.02"),
  Ovary = c("Ovary.01","Ovary.02"),
  TCell = c("TCell.01","TCell.02"),
  Thymus = c("Thymus.01","Thymus.02"),
  Trachea = c("Trachea.01","Trachea.02"),
  Pectoralis = c("B.muscle.01","B.muscle.02")
)

dat.expr=t(fpkm.expressed)
dat.expr[1:5,1:5]
fpkm.expressed[1:5,1:5]

##calculate the sum
# 
data_sum= list()
for (tissue in names(tissue_groups)) {
  columns=tissue_groups[[tissue]]
  col=as.data.frame(columns)
  expr= (colnames(fpkm.expressed) %in% col$columns)
  dt =as.matrix(fpkm.expressed[,expr])
  sum_tissue=rowSums(dt)
  data_sum[[tissue]]= sum_tissue
}

#create the data
colnames(dat.expr)
data_sum.all=as.data.frame(t(bind_cols(data_sum)))
colnames(data_sum.all)=colnames(dat.expr)
dt.sum.tissue=as.data.frame(t(data_sum.all))
##

head(dt.sum.tissue)
head(predict)

data.merge.sum=merge(dt.sum.tissue,predict, by.x="row.names", by.y="name")
head(data.merge.sum)

##mRNA
data.sum.mRNA=subset(data.merge.sum, type == "mRNA")
names(data.sum.mRNA)
data.sum.mRNA.filter=data.sum.mRNA[,-(22:25)]
row.names(data.sum.mRNA.filter) = data.sum.mRNA.filter[,1]
data.sum.mRNA.filter = data.sum.mRNA.filter[,-1]          
dt.mRNA2=t(data.sum.mRNA.filter)
dt.mRNA=as.data.frame(rowSums(dt.mRNA2))
colnames(dt.mRNA)= "mRNA"
##

#lnCRNA
data.sum.lRNA=subset(data.merge.sum, type == "lncRNA")
names(data.sum.lRNA)
data.sum.lRNA.filter=data.sum.lRNA[,-(22:25)]
row.names(data.sum.lRNA.filter) = data.sum.lRNA.filter[,1]
data.sum.lRNA.filter = data.sum.lRNA.filter[,-1]          
dt.lRNA2=t(data.sum.lRNA.filter)
dt.lRNA=as.data.frame(rowSums(dt.lRNA2))
colnames(dt.lRNA)= "lncRNA"
#

data.sum2=cbind(dt.mRNA, dt.lRNA)
head(data.sum2)
data.sum2=dplyr::as_tibble(data.sum2, rownames = "Tissue")
head(data.sum2)

##boxplot
col.d=c("BCell" ="#0DFAFA", "Bursa"="#004C99","Proximal.Cecum"="#CCFF99",
             "Macrophage(D0)" = "#F02B6D" , "Macrophage(D3)" = "#FF0091",
             "Macrophage(D6)" = "#F572BC", "Iliotibial.major"= "#A1122A",
             "Ileum" = "#98D55A", "Jejunum" ="#4C9900", "Kidney" ="#43009A",
             "Macrophage(lung)" = "#D06AAA", "Monocyte(blood)" = "#91155B", 
             "Magnum" ="#DE5100", "Isthmus"= "#FF8000", "Shell.Gland" = "#FFC78E",
             "Ovary"= "#CCCC00", "TCell(spleen)" = "#13B7B7", "Thymus" ="#2685E4",
             "Trachea"="#990099", "Pectoralis.major" = "#FFCCCC")


melted_df=reshape2::melt(data.sum2, id.vars = "Tissue", variable.name = "Isoform_type", 
                value.name = "Sum.Expression")

melted_df$Tissue[melted_df$Tissue == "Bcell"] = "BCell"
melted_df$Tissue[melted_df$Tissue == "Proximal_Cecum"] = "Proximal.Cecum"
melted_df$Tissue[melted_df$Tissue == "D0"] = "Macrophage(D0)"
melted_df$Tissue[melted_df$Tissue == "D3"] = "Macrophage(D3)"
melted_df$Tissue[melted_df$Tissue == "D6"] = "Macrophage(D6)"
melted_df$Tissue[melted_df$Tissue == "Pectoralis"] = "Pectoralis.major"
melted_df$Tissue[melted_df$Tissue == "Macrophage_lung"] = "Macrophage(lung)"
melted_df$Tissue[melted_df$Tissue == "Shell_Gland"] = "Shell.Gland"
melted_df$Tissue[melted_df$Tissue == "Monocyte"] = "Monocyte(blood)"
melted_df$Tissue[melted_df$Tissue == "TCell"] = "TCell(spleen)"
melted_df$Tissue[melted_df$Tissue == "Iliotibial"] = "Iliotibial.major"

head(melted_df)

#color iso
isoform_colors= c("#7D8ABC", "#304463")

p2= ggplot(melted_df, aes(x = Isoform_type, y = Sum.Expression, fill = Isoform_type)) +
  geom_boxplot() +
 ggtitle("Total of Transcript Isoform Expresison") +
  scale_fill_manual(values = isoform_colors) +
  ggthemes::theme_few() +
  theme01 +
  ylab("Transcript Isoform Expression") +
  xlab("") 

ggsave(here::here(out,"mRNA_lncRNA.svg"), plot = p2, width = 12, height = 10,
       units = "in", dpi = 300)
##

##heatmap mRNA and lnRNA
library(pheatmap)
dim(mRNA.expressed.filter) #1358   40
mRNA.expressed.filter[1:5,1:5]
dim(lncRNA.expressed.filter) #4322   40
lncRNA.expressed.filter [1:5,1:5]

# Annotation colors
annotation_colors= list(
  Tissues= c("BCell" ="#0DFAFA", "Bursa"="#004C99","Proximal.Cecum"="#CCFF99",
             "Macrophage(D0)" = "#F02B6D" , "Macrophage(D3)" = "#FF0091",
             "Macrophage(D6)" = "#F572BC", "Iliotibial.major"= "#A1122A",
             "Ileum" = "#98D55A", "Jejunum" ="#4C9900", "Kidney" ="#43009A",
             "Macrophage(lung)" = "#43009A", "Monocyte(blood)" = "#91155B", 
             "Magnum" ="#DE5100", "Isthmus"= "#FF8000", "Shell.Gland" = "#FFC78E",
             "Ovary"= "#CCCC00", "TCell(spleen)" = "#13B7B7", "Thymus" ="#2685E4",
             "Trachea"="#990099", "Pectoralis.major" = "#FFCCCC"),
  
  System = c("Immune" = "#31E1F7", "Respiratory" = "#400D51", "Excretory"="#FF7777",
           "Muscular" = "#D800A6","Intestine" ="#6499E9", "Reproductive" = "#836FFF")
  )

##
heat.cols=colorRampPalette(hodor_pal("black_m")(7))(250)

#heatmap 

#lncRNA
dat.lncRNA=t(lncRNA.expressed.filter)
#annotation= data.frame(Tissues = pheno.filter$Tissue)
annotation= data.frame(Tissues = pheno.filter$Tissue, System = pheno.filter$System)
rownames(annotation)= colnames(lncRNA.expressed.filter)

plncRNA=pheatmap(
   lncRNA.expressed.filter, 
    color =heat.cols ,
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    scale="row",
    drop_levels = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    #legend_breaks = c(-6, 0, 6),
    #legend_labels = c("-6", "0" ,"6"),
    legend = TRUE,
    main = "lncRNAs Expressed",
    fontsize_col = 14,
    fontsize_row = 14
  )

png(file = "isoforms_lncRNA01.png", width = 16, height = 8, units = "in", res = 300)
plncRNA
dev.off()

##
#mRNA heatmap
#lncRNA
dat.mRNA=t(mRNA.expressed.filter)
#annotation= data.frame(Tissues = pheno.filter$Tissue)
annotation= data.frame(Tissues = pheno.filter$Tissue, System = pheno.filter$System)
rownames(annotation)= colnames(mRNA.expressed.filter)

pmRNA=pheatmap(
  mRNA.expressed.filter, 
  color =heat.cols ,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  scale="row",
  drop_levels = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  #legend_breaks = c(-6, 0, 6),
  #legend_labels = c("-6", "0" ,"6"),
  legend = TRUE,
  main = "mRNAs Expressed",
  fontsize_col = 14,
  fontsize_row = 14
)

png(file = "isoforms_mRNA01.png", width = 16, height = 8, units = "in", res = 300)
pmRNA
dev.off()

##end

##heatmap and boxplot U
##
#boxplot mRNA vcs LncRNA
##filter data
head(predict)
lncRNA.predict.u=subset(predict, type == "lncRNA" &  code == "u")
mRNA.predict.u=subset(predict, type == "mRNA"&  code == "u")

mRNA.expressed.u= which(row.names(fpkm.expressed) %in% mRNA.predict.u$name)
mRNA.expressed.filter.u=fpkm.expressed[mRNA.expressed.u, ]
mRNA.expressed.filter.u[1:5,1:5]
dim(mRNA.expressed.filter.u) #872  40

lncRNA.expressed.u= which(row.names(fpkm.expressed) %in% lncRNA.predict.u$name)
lncRNA.expressed.filter.u=fpkm.expressed[lncRNA.expressed.u, ]
lncRNA.expressed.filter.u[1:5,1:5]
dim(lncRNA.expressed.filter.u) #1334   40

#vector tissues 
tissue_groups= list(
  Bcell = c("Bcell.01","BCell.02"),
  Bursa = c("Bursa.01","Bursa.02"),
  Proximal_Cecum = c("P.cecum.01","P.cecum.02"),
  D0 = c("Macrophage.D0.01","Macrophage.D0.02"),
  D3 = c("Macrophage.D3.01","Macrophage.D3.02"),
  D6 = c("Macrophage.D6.01","Macrophage.D6.02"),
  Iliotibial = c("T.muscle.01", "T.muscle.02"),
  Ileum = c("Ileum.01","Ileum.02"),
  Jejunum = c("Jejunum.01","Jejunum.02"),
  Kidney = c("Kidney.01","Kidney.02"),
  Macrophage_lung = c("Macrophage.lung.01","Macrophage.lung.02"),
  Monocyte = c("Monocyte.01","Monocyte.02"),
  Magnum = c("Magnum.od1.01","Magnum.od1.01"),
  Isthmus = c("Isthmus.od2.01","Isthmus.od2.02"),
  Shell_Gland = c("Shell.Gland.od3.01","Shell.Gland.od3.02"),
  Ovary = c("Ovary.01","Ovary.02"),
  TCell = c("TCell.01","TCell.02"),
  Thymus = c("Thymus.01","Thymus.02"),
  Trachea = c("Trachea.01","Trachea.02"),
  Pectoralis = c("B.muscle.01","B.muscle.02")
)


##heatmap mRNA and lnRNA
library(pheatmap)
dim(mRNA.expressed.filter.u) #872  40
mRNA.expressed.filter.u[1:5,1:5]
dim(lncRNA.expressed.filter.u) #1334   40
lncRNA.expressed.filter.u [1:5,1:5]

# Annotation colors
annotation_colors= list(
  Tissues= c("BCell" ="#0DFAFA", "Bursa"="#004C99","Proximal.Cecum"="#CCFF99",
             "Macrophage(D0)" = "#F02B6D" , "Macrophage(D3)" = "#FF0091",
             "Macrophage(D6)" = "#F572BC", "Iliotibial.major"= "#A1122A",
             "Ileum" = "#98D55A", "Jejunum" ="#4C9900", "Kidney" ="#43009A",
             "Macrophage(lung)" = "#43009A", "Monocyte(blood)" = "#91155B", 
             "Magnum" ="#DE5100", "Isthmus"= "#FF8000", "Shell.Gland" = "#FFC78E",
             "Ovary"= "#CCCC00", "TCell(spleen)" = "#13B7B7", "Thymus" ="#2685E4",
             "Trachea"="#990099", "Pectoralis.major" = "#FFCCCC"),
  System = c("Immune" = "#31E1F7", "Respiratory" = "#400D51", "Excretory"="#FF7777",
             "Muscular" = "#D800A6","Intestine" ="#6499E9", "Reproductive" = "#836FFF")
)

##
heat.cols=colorRampPalette(hodor_pal("black_m")(7))(250)

#heatmap 

#lncRNA
#annotation= data.frame(Tissues = pheno.filter$Tissue)
annotation= data.frame(Tissues = pheno.filter$Tissue, System = pheno.filter$System)
rownames(annotation)= colnames(lncRNA.expressed.filter.u)

plncRNA.u=pheatmap(
  lncRNA.expressed.filter.u, 
  color =heat.cols ,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  scale="row",
  drop_levels = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  #legend_breaks = c(-6, 0, 6),
  #legend_labels = c("-6", "0" ,"6"),
  legend = TRUE,
  main = "lncRNAs Expressed",
  fontsize_col = 14,
  fontsize_row = 14
)

png(file = "isoforms_lncRNA.u01.png", width = 16, height = 8, units = "in", res = 300)
plncRNA.u
dev.off()

##
#mRNA heatmap
#annotation= data.frame(Tissues = pheno.filter$Tissue)
annotation= data.frame(Tissues = pheno.filter$Tissue, System = pheno.filter$System)
rownames(annotation)= colnames(mRNA.expressed.filter.u)

pmRNA.u=pheatmap(
  mRNA.expressed.filter.u, 
  color =heat.cols ,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  scale="row",
  drop_levels = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  #legend_breaks = c(-6, 0, 6),
  #legend_labels = c("-6", "0" ,"6"),
  legend = TRUE,
  main = "mRNAs Expressed",
  fontsize_col = 14,
  fontsize_row = 14
)

png(file = "isoforms_mRNA.u01.png", width = 16, height = 8, units = "in", res = 300)
pmRNA.u
dev.off()

##end
#tissue especific only u
RNA.u=subset(predict, code == "u")
RNA.expressed.u= which(row.names(fpkm.expressed) %in% RNA.u$name)
RNA.expressed.filter.u=fpkm.expressed[RNA.expressed.u, ]
RNA.expressed.filter.u[1:5,1:5]
dim(RNA.expressed.filter.u) #2206   40


##heatmap only U total
annotation= data.frame(Tissues = pheno.filter$Tissue, System = pheno.filter$System)
rownames(annotation)= colnames(RNA.expressed.filter.u)

RNA.u=pheatmap(
  RNA.expressed.filter.u, 
  color =heat.cols ,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  scale="row",
  drop_levels = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  #legend_breaks = c(-6, 0, 6),
  #legend_labels = c("-6", "0" ,"6"),
  legend = TRUE,
  main = "Transcript Expressed",
  fontsize_col = 14,
  fontsize_row = 14
)

png(file = "isoforms_totalU.png", width = 16, height = 8, units = "in", res = 300)
RNA.u
dev.off()


#top5
fpkm.expressed[1:5,1:5]
head(top)

#selected ony the top

mRNA.expressed.top= which(row.names(fpkm.expressed) %in% top$transcript)
mRNA.expressed.filter.top=fpkm.expressed[mRNA.expressed.top, ]
mRNA.expressed.filter.top[1:5,1:5]
dim(mRNA.expressed.filter.top) #302  40
##
mRNA.expressed.filter.top[1:5,1:5]
topU_expressed=data.frame(row.names(mRNA.expressed.filter.top))
colnames(topU_expressed)= "transcript"
head(topU_expressed)

write.table(topU_expressed, file= "topU_expressed.txt", row.names = F, sep = "\t",
            na= "NA", quote = F)

##heatmap top
library(pheatmap)
dim(mRNA.expressed.filter.top) #302  40


# Annotation colors
annotation_colors= list(
  Tissues= c("BCell" ="#0DFAFA", "Bursa"="#004C99","Proximal.Cecum"="#CCFF99",
             "Macrophage(D0)" = "#F02B6D" , "Macrophage(D3)" = "#FF0091",
             "Macrophage(D6)" = "#F572BC", "Iliotibial.major"= "#A1122A",
             "Ileum" = "#98D55A", "Jejunum" ="#4C9900", "Kidney" ="#43009A",
             "Macrophage(lung)" = "#43009A", "Monocyte(blood)" = "#91155B", 
             "Magnum" ="#DE5100", "Isthmus"= "#FF8000", "Shell.Gland" = "#FFC78E",
             "Ovary"= "#CCCC00", "TCell(spleen)" = "#13B7B7", "Thymus" ="#2685E4",
             "Trachea"="#990099", "Pectoralis.major" = "#FFCCCC"),
  System = c("Immune" = "#31E1F7", "Respiratory" = "#400D51", "Excretory"="#FF7777",
             "Muscular" = "#D800A6","Intestine" ="#6499E9", "Reproductive" = "#836FFF")
)

##
heat.cols=colorRampPalette(hodor_pal("black_m")(7))(250)

#heatmap 
annotation= data.frame(Tissues = pheno.filter$Tissue, System = pheno.filter$System)
rownames(annotation)= colnames(mRNA.expressed.filter.top)

topU=pheatmap(
  mRNA.expressed.filter.top, 
  color =heat.cols ,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  scale="row",
  drop_levels = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  #legend_breaks = c(-6, 0, 6),
  #legend_labels = c("-6", "0" ,"6"),
  legend = TRUE,
  main = "Top 5 Expressed",
  fontsize_col = 14,
  fontsize_row = 14
)

png(file = "isoforms_topRNAU.png", width = 16, height = 8, units = "in", res = 300)
topU
dev.off()

##
head(domain)
data.iso.top=as.data.frame(unique(domain$name))
dim(data.iso.top) ##transcript 309 and 87 genes
write.table(data.iso.top, file= "data.iso.top.txt", row.names = F, sep = "\t",
            na= "NA", quote = F)

##
##domain
head(domain)
length(unique(domain$Signature_domain))
table(domain$Signature_domain)

##data screenshot

transcript = data.frame(c("MSTRG.20945.2", "MSTRG.20945.5")) 
colnames(transcript) = "transcript"

mRNA.mstrg= which(row.names(fpkm) %in% transcript$transcript)
mRNA.expressed.mstrg=fpkm[mRNA.mstrg, ]
##

#vector tissues 
tissue_groups= list(
  BCell = c("Bcell.01","BCell.02"),
  Bursa = c("Bursa.01","Bursa.02"),
  "Proximal.Cecum" = c("P.cecum.01","P.cecum.02"),
  "Macrophage(D0)" = c("Macrophage.D0.01","Macrophage.D0.02"),
  "Macrophage(D3)" = c("Macrophage.D3.01","Macrophage.D3.02"),
  "Macrophage(D6)" = c("Macrophage.D6.01","Macrophage.D6.02"),
  Iliotibial.major = c("T.muscle.01", "T.muscle.02"),
  Ileum = c("Ileum.01","Ileum.02"),
  Jejunum = c("Jejunum.01","Jejunum.02"),
  Kidney = c("Kidney.01","Kidney.02"),
  "Macrophage(lung)" = c("Macrophage.lung.01","Macrophage.lung.02"),
  "Monocyte(blood)" = c("Monocyte.01","Monocyte.02"),
  Magnum = c("Magnum.od1.01","Magnum.od1.01"),
  Isthmus = c("Isthmus.od2.01","Isthmus.od2.02"),
  "Shell.Gland" = c("Shell.Gland.od3.01","Shell.Gland.od3.02"),
  Ovary = c("Ovary.01","Ovary.02"),
  "TCell(spleen)" = c("TCell.01","TCell.02"),
  Thymus = c("Thymus.01","Thymus.02"),
  Trachea = c("Trachea.01","Trachea.02"),
  "Pectoralis.major" = c("B.muscle.01","B.muscle.02")
)

##cal mean
data_mean= list()
for (tissue in names(tissue_groups)) {
  columns=tissue_groups[[tissue]]
  col=as.data.frame(columns)
  expr= (colnames(mRNA.expressed.mstrg) %in% col$columns)
  dt =as.matrix(mRNA.expressed.mstrg[,expr])
  mean_tissue=rowMeans(dt)
  data_mean[[tissue]]= mean_tissue
}

data_mean1=as.data.frame(t(bind_cols(data_mean)))
colnames(data_mean1)=c("MSTRG.20945.2", "MSTRG.20945.5")
data_mean1=tibble::rownames_to_column(data_mean1, "Tissues")
write.table(data_mean1, file= "data_mean_iso_screen.txt", row.names = F, sep = "\t",
            na= "NA", quote = F)

##
#color iso
melted_dt=reshape2::melt(data_mean1, id.vars = "Tissues", variable.name = "Transcript_Isoform", 
                         value.name = "Mean.Expression")
colors1= c("#B50698", "#8E94E1")

p6= ggplot(melted_dt, aes(x = Transcript_Isoform, y = Mean.Expression, fill = Transcript_Isoform)) +
  #geom_boxplot() +
  geom_violin(trim=FALSE) +
  ggtitle("") +
  scale_fill_manual(values = colors1) +
  ggthemes::theme_few() +
  theme01 +
  ylab("Expression (FPKM)") +
  xlab("") 

ggsave(here::here(out,"screnshot_mRNA.svg"), plot = p6, width = 12, height = 10,
       units = "in", dpi = 300)
##



