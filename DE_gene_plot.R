##***************************************************************************
##*DEGs Genes Oviduct 
##***************************************************************************
##colors tissues vector
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
  `Ovary`="#CCCC00")

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

##*****************************************************************************
##*libraries
library(ggplot2)
library(dplyr)
library(plotly)
library(hrbrthemes)
library(reshape2)
library(wesanderson)
##******************************************************************************
##*##data path
input=file.path("D:", "RNA-seq", "Reproduction_data", "DEG", "DEG",
                fsep="/")
out=file.path("D:", "RNA-seq", "Reproduction_data", "DEG", "DEGv2_gene",
              fsep="/")
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##DE genes
#DE1 - od12
od12= read.csv(here::here(input,"DEdataod12.csv"),
               header =T, row.names = 1)
#threshold column
od12$abs=abs(od12$log2FoldChange)
od12$threshold <- od12$padj <= 0.05 & od12$abs >= 1.5
#DE2 - od13
od13= read.csv(here::here(input,"DEdataod13.csv"),
               header =T, row.names = 1)
#threshold column
od13$abs=abs(od13$log2FoldChange)
od13$threshold <- od13$padj <= 0.05 & od13$abs >= 1.5
#DE3 - od23
od23= read.csv(here::here(input,"DEdataod23.csv"),
               header =T, row.names = 1)
#threshold column
od23$abs=abs(od23$log2FoldChange)
od23$threshold <- od23$padj <= 0.05 & od23$abs >= 1.5
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##*
library(limma)
library(edgeR)

cnt= read.table(here::here(input,"all_cnt_0625.txt"),header=T)
cnt[1:5,1:5]
##
pheno=read.table(here::here(input,"pheno_2_edit.txt"),header=T)

#transform rownames
row.names(cnt) = cnt[,1] 
head(cnt)
##
cnt = cnt[,-1]          
dim(cnt) 
cnt[1:5,1:5]
##

##filter zero
#remove genes non-expressed and low expressed
Totalcounts = rowSums(cnt)
## genes with zero count?
table(Totalcounts==0)
## filter genes with 0 counts
rm = rowMeans(cnt)==0
expr.filter = cnt[!rm,]
dim(expr.filter) 
##********************************************************************
##normalization using TMM
d = DGEList(expr.filter)
TMM = calcNormFactors(d, method = "TMM")
head(TMM$samples)
#get the normalized counts:
cpm = cpm(TMM, log=FALSE)
cpm=log(cpm+1, base=2)  #log scale
cpm[1:5,1:5]
##**********************************************************************
##volcano plot ods
#Adding a -log10(p-value) column for each df 
od12$neglog10FDR = -log10(od12$padj)
od13$neglog10FDR = -log10(od13$padj)
od23$neglog10FDR = -log10(od23$padj)

#PLOTTING VOLCANO PLOTS 
#identifying top -log10(p-values)
top_positiveod12 <-
  od12[od12$log2FoldChange > 0, ][order(-od12$neglog10FDR), ][1:10, ]
top_negativeod12 <-
  od12[od12$log2FoldChange < 0, ][order(-od12$neglog10FDR), ][1:10, ]
top_positiveod13 <- 
  od13[od13$log2FoldChange > 0, ][order(-od13$neglog10FDR), ][1:10, ]
top_negativeod13 <-
  od13[od13$log2FoldChange < 0, ][order(-od13$neglog10FDR), ][1:10, ]
top_positiveod23 <- 
  od23[od23$log2FoldChange > 0, ][order(-od23$neglog10FDR), ][1:10, ]
top_negativeod23 <- 
  od23[od23$log2FoldChange < 0, ][order(-od23$neglog10FDR), ][1:10, ]

##saving the top 20
data.top=rbind(top_positiveod12,top_negativeod12,top_positiveod13,top_negativeod13,
               top_positiveod23, top_negativeod23)

data_top.gene=as.data.frame(unique(data.top[,c("Gene")]))
colnames(data_top.gene)= "Gene"

#plotting 
od12.color=mycolors1[c("Magnum.od1","Isthmus.od2")] #order up and down 
od12volcano=ggplot(od12, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(padj >= 0.05 & abs <= 1.5 , "grey",
                                ifelse(log2FoldChange > 0, 
                                       od12.color["Magnum.od1"], 
                                       od12.color["Isthmus.od2"]))), 
             size = 1,na.rm=TRUE) +
  scale_color_identity(guide = "legend", "Tissues", labels=c("Up-Regulated",
                                                             "Down-Regulated",
                                                             "Non-Significant","") 
                       ) +
  labs(x = "log2 Fold Change", y = "-log10(p-value)",
       title = "DEGs Magnum vs. Isthmus") +
  ggthemes::theme_few() + 
  #theme_classic() +
  #theme_ipsum() +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(color = 'white', fill = 'white'),
    panel.background = element_rect(color = 'white', fill = 'white'),
    text = element_text(size = 12, family = "Sans"),  # Set the base text size
    axis.title = element_text(size = 12, family = "Sans"),  # Set axis title size
    axis.text = element_text(size = 12,family = "Sans", hjust = 1),
    legend.title = element_text(size = 12, family = "Sans", face= "bold"), # Adjust legend title size
    legend.text = element_text(size = 12, family = "Sans")  # Adjust legend text size
  ) +
  # Add labels for the top 10 positive log2FC values
  geom_text(data = top_positiveod12, aes(label = Gene), vjust = -0.5, 
            hjust = 1, color = "black", size = 3,check_overlap = TRUE) +
  # Add labels for the top 10 negative log2FC values
  geom_text(data = top_negativeod12, aes(label = Gene), vjust = -0.5, 
            hjust = -0.5, color = "black", size = 3,check_overlap = TRUE)
##plot2
od13.color=mycolors1[c("Magnum.od1","Shell.Gland.od3")] #order up and down 
od13volcano=ggplot(od13, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(padj >= 0.05 & abs <= 1.5, "grey",
                                ifelse(log2FoldChange > 0, 
                                       mycolors1["Magnum.od1"], 
                                       mycolors1["Shell.Gland.od3"]))), 
             size = 1,na.rm =TRUE) +
scale_color_identity(guide = "legend", "Tissues", labels=c("Up-Regulated",
                                                             "Down-Regulated",
                                                             "Non-Significant","")) +
  labs(x = "log2 Fold Change", y = "-log10(p-value)", 
       title = "DEGs Magnum vs. Shell Gland") +
  ggthemes::theme_few() + 
  #theme_classic() +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(color = 'white', fill = 'white'),
    panel.background = element_rect(color = 'white', fill = 'white'),
    text = element_text(size = 12, family = "Sans"),  # Set the base text size
    axis.title = element_text(size = 12, family = "Sans"),  # Set axis title size
    axis.text = element_text(size = 12,family = "Sans", hjust = 1),
    legend.title = element_text(size = 12, family = "Sans", face= "bold"), # Adjust legend title size
    legend.text = element_text(size = 12, family = "Sans")  # Adjust legend text size
  ) + 
  # Add labels for the top 10 positive log2FC values
  geom_text(data = top_positiveod13, aes(label = Gene), 
            vjust = -0.5, hjust = 1, color = "black", size = 3,check_overlap = TRUE) +
  # Add labels for the top 10 negative log2FC values
  geom_text(data = top_negativeod13, aes(label = Gene), 
            vjust = -0.5, hjust = -0.5, color = "black", size = 3,check_overlap = TRUE)

#plot3
od23.color=mycolors1[c("Shell.Gland.od3","Isthmus.od2")] #order up and down 
od23volcano=ggplot(od23, aes(x = log2FoldChange, 
                             y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(padj >= 0.05 & abs <= 1.5, "grey",
                                ifelse(log2FoldChange > 0, 
                                       mycolors1["Shell.Gland.od3"], 
                                       mycolors1["Isthmus.od2"]))), 
             size = 1,na.rm = TRUE) +
  scale_color_identity(guide = "legend", "Tissues", labels=c("Down-Regulated",
                                                             "Up-Regulated",
                                                             "Non-Significant","")) +
  labs(x = "log2 Fold Change", y = "-log10(p-value)",
       title = "DEGs Shell Gland vs. Isthmus") +
  ggthemes::theme_few() +
  #theme_classic() +
  #theme_ipsum() +
  theme(
    legend.position = "bottom",
    plot.background = element_rect(color = 'white', fill = 'white'),
    panel.background = element_rect(color = 'white', fill = 'white'),
    text = element_text(size = 12, family = "Sans"),  # Set the base text size
    axis.title = element_text(size = 12, family = "Sans"),  # Set axis title size
    axis.text = element_text(size = 12,family = "Sans", hjust = 1),
    legend.title = element_text(size = 12, family = "Sans", face= "bold"), # Adjust legend title size
    legend.text = element_text(size = 12, family = "Sans")  # Adjust legend text size
  ) + 
  # Add labels for the top 10 positive log2FC values
  geom_text(data = top_positiveod23, aes(label = Gene), 
            vjust = -0.5, hjust = 1, color = "black", size = 3,check_overlap = TRUE) +
  # Add labels for the top 10 negative log2FC values
  geom_text(data = top_negativeod23, aes(label = Gene), 
            vjust = -0.5, hjust = -0.5, color = "black", size = 3,check_overlap = TRUE)
##
#od12volcano <- volcano_plot(od12, 
 #                           mycolors1[c("Magnum.od1", "Isthmus.od2")], 
#                            "DEGs Magnum vs. Isthmus")
#od13volcano <- volcano_plot(od13, 
#                            mycolors1[c("Magnum.od1", "Shell.Gland.od3")], 
#                            "DEGs Magnum vs. Shell Gland")
#od23volcano <- volcano_plot(od23, 
 #                           mycolors1[c("Shell.Gland.od3","Isthmus.od2")],
#                            "DEGs Isthmus vs. Shell Gland")

#making plot list 
volcanoplots=list(od12volcano,od13volcano,od23volcano)
plot=ggpubr::ggarrange(plotlist = volcanoplots, ncol = 3, nrow = 1)
ggsave(here::here(out,"volcanos_ods2.svg"), plot = plot, width = 18, height = 14, 
       units = "in", dpi = 600)
##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##heatmap
head(data_top.gene)
cpm[1:5,1:5]

ods_data= cpm[, c("Sample_od1_rep4","Sample_od1_rep5","Sample_od2_rep2",
                  "Sample_od2_rep6","Sample_od3_rep1","Sample_od3_rep5")]
colnames(ods_data)=c("Magnum.r1", "Magnum.r2","Isthmus.r1","Isthmus.r2",
                     "Shell_Gland.r1", "Shell_Gland.r2")

dat1=read.table(here::here(input,"pheno.od.txt"), header=T)
##
ods_data[1:5,1:5]
head(data_top.gene)
data.gene=data_top.gene
data.gene$Gene <- paste0("gene-", data.gene$Gene)

ods_data.top= as.data.frame(ods_data[which(rownames(ods_data)  %in% 
                                             data.gene$Gene),])

data.m1=as.matrix(t(ods_data.top))
dim(data.m1)
##heatmap
head(dat1)
dim(dat1)
mycolors <- colorRampPalette(wes_palette("Zissou1"))(100)
annotation= data.frame(Tissue=dat1[,'Tissue'])
row.names(annotation)= colnames(ods_data.top)

##colors tissue
mycol1=list(
       Tissue= c( "Magnum" = "#DE5100",
        "Isthmus" =  "#FF8000",                    
        "Shell_Gland" = "#FFC78E"))

### Run pheatmap
p4=pheatmap::pheatmap(ods_data.top, color = 
              mycolors, cluster_rows = T,
            show_rownames=T,
            #annotation = annotation, 
            border_color=NA, fontsize = 10, 
            scale="row", 
            show_colnames = F,
            #show_rownames = F,
            annotation_col= annotation,
            fontsize_row = 10, height=20,
            cluster_cols = F, 
            #cluster_rows = F,
            annotation_colors = mycol1,
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean", 
            clustering_method = "ward.D",
            legend_breaks = c(-2, 0, 2),
            legend_labels = c("Low", "Medium", "High"),
            annotation_legend = TRUE,
            main = "")

ggsave(here::here(out,"heatmap_ods4.svg"), plot = p4, width = 5, height = 8, 
       units = "in", dpi = 300)
##

##Venn diagram 
library(VennDiagram)
names(od12)
od12.sig= as.data.frame(od12[od12$threshold == "TRUE",])
od12.sig=na.omit(od12.sig)
Magnum_Isthmus= as.data.frame(od12.sig$Gene)
colnames(Magnum_Isthmus) = "Gene"
dim(Magnum_Isthmus)
##
names(od13)
od13.sig= as.data.frame(od13[od13$threshold == "TRUE",])
od13.sig=na.omit(od13.sig)
Magnum_Shell.Gland= as.data.frame(od13.sig$Gene)
colnames(Magnum_Shell.Gland) = "Gene"
dim(Magnum_Shell.Gland)
##
names(od23)
od23.sig= as.data.frame(od23[od23$threshold == "TRUE",])
od23.sig=na.omit(od23.sig)
Shell.Gland_Isthmus= as.data.frame(od23.sig$Gene)
colnames(Shell.Gland_Isthmus) = "Gene"
dim(Shell.Gland_Isthmus)
##

gene_sets= list(
  "od12" = as.vector(Magnum_Isthmus$Gene), #"Magnum-Isthmus"
  "od13" = as.vector(Magnum_Shell.Gland$Gene), #"Magnum-Shell_Gland"
  "od23" = as.vector(Shell.Gland_Isthmus$Gene) #"Shell_Gland-Isthmus"
)
##plot
mycol6=c("#FFD014","#CC6600","#CC0000")
venn.diagram( gene_sets,
    #category.names = c( "Magnum-Isthmus", "Magnum-Shell.Gland" ,
     #                   "Shell.Gland_Isthmus"),
    category.names = c( "DEGs.01", "DEGs.02" ,
                       "DEGs.03"),
    imagetype = "png", 
    filename = here::here(out,"ods5.png"),
    lwd = 2,
    #lty = "blank",
    col=c("#FFD014","#CC6600","#CC0000"),
    fill = c(alpha("#FFD014",0.1), alpha('#CC6600',0.1), alpha('#CC0000',0.1)),
    cat.cex = 0.50,
    height = 720 , 
    width = 720 , 
    resolution = 600,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = mycol6,
  label.col = "black",
  rotation = 1,
)


  
  









