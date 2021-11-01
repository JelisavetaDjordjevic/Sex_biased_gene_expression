library(readr)
library(edgeR)
library(pheatmap)
library(reshape2)

#read in "all_stages_dnds.txt" from:
#"https://github.com/JelisavetaDjordjevic/Sex_biased_gene_expression/tree/main/Data/all_stages_dnds.txt")

all_stages<- read.delim("~/all_stages_dnds.txt")

#get 26 genes that are sex-biased at hatchling stage
hatch_sb<- all_stages[all_stages$stage=="hatch"& all_stages$fact != "Not DE",]
hatch_sb_genes<- as.character(hatch_sb$genes)
hatch_sb_genes
#get the FC vlaues for these 26 genes in other two developmental stages
sub<- all_stages[all_stages$genes%in% hatch_sb_genes, ]

##lets make a heatmap of 26 genes from hatchling in juv and adult stage
sub_heat<-sub[,c(1,2,5)]
sub_heat_mat<-acast(sub_heat, genes~stage,  na.rm = FALSE, value.var ="logFC")

## set breaks ## 
breaksList = c(-4.25, -3.75, -3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25, 0.25,  0.75,  1.25, 1.75, 2.25, 2.75,3.25,3.75,4.25) 
pheatmap(sub_heat_mat, 
         clustering_distance_rows = "euclidean", 
         cluster_cols = FALSE, 
         clustering_method = "complete", 
         color = colorRampPalette(c("#08103A","#08306B","#08417C", "#08519C", "#2171B5", "#4292C6", "#6BAED6" ,"#9DCBE1","#9ECAE1", "#FFFFFF", "#FFFFFF"  ,"#FCBBA1","#FCBBA1" ,"#FC9272" ,"#FB6A4A" ,"#EF3B2C" ,"#CB181D" ,"#A50F15" ,"#67000D","#67000D","#67000D"))(length(breaksList)), 
         breaks = breaksList, 
         show_rownames = F,border_color = NA ,
         main  = paste("GN | overeff - inter top10% |  N =", length(sub_heat_mat[,1])))
