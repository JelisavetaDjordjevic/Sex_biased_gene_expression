###packages
library(readr)
library(edgeR)
library(ggplot2)
library(VennDiagram)
library(SuperExactTest)
getwd()

#read droso_counts.csv
droso_counts<-read_csv("droso_counts.csv")

### Larva ####
###3rd instar Larva###
counts_larvae<- droso_counts[, c(6,8,12,14,17,19,21,23)]
colnames(counts_larvae)
counts_larvae$Gene_name<- droso_counts$Gene_name
sex_l<-c("F","F","F","M","M","F","M","M") 
#make edgeR object for larval stage
x_l_stage<-DGEList(counts = droso_counts[, c(6,8,12,14,17,19,21,23)],
                   genes = droso_counts[, 1],
                   group = sex_l)
x_l_stage$samples

#normalization for library sizes
cpm(10, mean(x_l_stage$samples$lib.size))
x_l_stage
head(str(x_l_stage))
#filter out genes that are not expressed in at least 3 male or 3 female libraries 
keep_l <- rowSums(cpm(x_l_stage[,c(1,2,3,6)])>0.5) >=3 | 
  rowSums(cpm(x_l_stage[,c(4,5,7,8)])>0.5) >=3  

x_l_stage<- x_l_stage[keep_l, , keep.lib.sizes=FALSE]

#show number of genes kept and filtered out
table(keep_l)
#calculate normalization factor
x_l_stage<- calcNormFactors(x_l_stage)

#design a model
design_l <- model.matrix(~ 0 + group, data= x_l_stage$samples)
colnames(design_l)<- levels(x_l_stage$samples$group)
colnames(design_l)

#estimate dispersion
x_l_stage <- estimateDisp(x_l_stage, design_l, robust=TRUE)

####Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.
fit_l <- glmQLFit(x_l_stage, design_l, robust=TRUE)

#make contrast between m and F
FM_l<- makeContrasts(F - M, levels=design_l)

#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
##Conduct genewise statistical tests for a given coefficient or contrast.
res_FM_l <- glmQLFTest(fit_l, contrast=FM_l)

#Identify which genes are significantly differentially expressed from an edgeR fit object
is.de_FM_l <- decideTestsDGE(res_FM_l, p.value=0.05)

#number of Sex-biased genes 
summary(is.de_FM_l)

#for venn
##Up, down and not DE genes at larval stage
DE_FM_l<- cbind(is.de_FM_l, res_FM_l$genes)
DE_FM_l<- as.data.frame(DE_FM_l)
head(DE_FM_l)

DE_l_f_up<-DE_FM_l[DE_FM_l$`1*F -1*M`=="1",]
DE_l_f_down<-DE_FM_l[DE_FM_l$`1*F -1*M`=="-1",]
DE_l<- rbind(DE_l_f_up, DE_l_f_down)
#Matrix for venn diagram, of DE genes at hathcling stage
DE_l<-as.matrix(DE_l$Gene_name)
DE_l

###pupae###
counts_pupae<- droso_counts[, c(2,3,4,13,15,16,24,25)]
colnames(counts_pupae)
counts_pupae$Gene_name<- droso_counts$Gene_name
sex_p<-c("F","F","F","M","M","M","M","F") 

#make edgeR object for pupal stage
x_p_stage<-DGEList(counts = droso_counts[, c(2,3,4,13,15,16,24,25)],
                   genes = droso_counts[, 1],
                   group = sex_p)
x_p_stage$samples

#normalization
cpm(10, mean(x_p_stage$samples$lib.size))

#filter out genes that are not expressed in at least 3 male or 3 female libraries 
keep_p <- rowSums(cpm(x_p_stage[,c(1,2,3,8)])>0.5) >=3 | 
  rowSums(cpm(x_p_stage[,c(4,5,6,7)])>0.5) >=3  
x_p_stage<- x_p_stage[keep_p, , keep.lib.sizes=FALSE]

#show number of genes kept and filtered out
table(keep_p)

#calculate normalization factor
x_p_stage<- calcNormFactors(x_p_stage)
x_p_stage$samples

#design a model
design_p <- model.matrix(~ 0 + group, data= x_p_stage$samples)
colnames(design_p)<- levels(x_p_stage$samples$group)
colnames(design_p)

#estimate dispersion
x_p_stage <- estimateDisp(x_p_stage, design_p, robust=TRUE)

####Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.
fit_p <- glmQLFit(x_p_stage, design_p, robust=TRUE)

#between m and F
FM_p<- makeContrasts(F - M, levels=design_p)

#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
##Conduct genewise statistical tests for a given coefficient or contrast.
res_FM_p <- glmQLFTest(fit_p, contrast=FM_p)

#Identify which genes are significantly differentially expressed from an edgeR fit object
is.de_FM_p <- decideTestsDGE(res_FM_p, p.value=0.05)

#number of Sex-biased genes
summary(is.de_FM_p)

#get fc values and p-values for (n) genes
pupae_top_fm<- topTags(res_FM_p, n=12000)
pupae_top_fm<-pupae_top_fm$table

#for venn
##Up, down and not DE genes at hathcling stage
DE_FM_p<- cbind(is.de_FM_p, res_FM_p$genes)
DE_FM_p<- as.data.frame(DE_FM_p)
head(DE_FM_p)

DE_p_f_up<-DE_FM_p[DE_FM_p$`1*F -1*M`=="1",]
DE_p_f_down<-DE_FM_p[DE_FM_p$`1*F -1*M`=="-1",]
DE_p<- rbind(DE_p_f_up, DE_p_f_down)
#Matrix for venn diagram, of DE genes at hathcling stage

DE_p<-as.matrix(DE_p$Gene_name)
DE_p

### Adult ###
#4 males and 4 females!
counts_adult<- droso_counts[, c(5,7,9,10,11,20,22,26)]
colnames(counts_adult)
counts_adult$Gene_name<- droso_counts$Gene_name
sex_a<-c("F","F","M","M","F","M", "F", "M") 
x_a_stage<-DGEList(counts = droso_counts[, c(5,7,9,10,11,20,22,26)],
                   genes = droso_counts[, 1],
                   group = sex_a)
x_a_stage$samples

#normalization
cpm(10, mean(x_a_stage$samples$lib.size))
x_a_stage
head(str(x_a_stage))
keep_a <- rowSums(cpm(x_a_stage[,c(1,2,5,7)])>0.5) >=3 | 
  rowSums(cpm(x_a_stage[,c(3,4,6,8)])>0.5) >=3  

x_a_stage<- x_a_stage[keep_a, , keep.lib.sizes=FALSE]
table(keep_a)
x_a_stage<- calcNormFactors(x_a_stage)



#design a model
design_a <- model.matrix(~ 0 + group, data= x_a_stage$samples)
colnames(design_a)<- levels(x_a_stage$samples$group)
colnames(design_a)
x_a_stage <- estimateDisp(x_a_stage, design_a, robust=TRUE)


####Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.
fit_a <- glmQLFit(x_a_stage, design_a, robust=TRUE)


#between m and F
FM_a<- makeContrasts(F - M, levels=design_a)
res_FM_a <- glmQLFTest(fit_a, contrast=FM_a)
is.de_FM_a <- decideTestsDGE(res_FM_a, p.value=0.05)
summary(is.de_FM_a)
adult_top_fm<- topTags(res_FM_a, n=14000)
adult_top_fm<-adult_top_fm$table

##Up, down and not DE genes at hathcling stage
DE_FM_a<- cbind(is.de_FM_a, res_FM_a$genes)
DE_FM_a<- as.data.frame(DE_FM_a)
head(DE_FM_a)

DE_a_f_up<-DE_FM_a[DE_FM_a$`1*F -1*M`=="1",]
DE_a_f_down<-DE_FM_a[DE_FM_a$`1*F -1*M`=="-1",]
DE_a<- rbind(DE_a_f_up, DE_a_f_down)
#Matrix for venn diagram, of DE genes at hathcling stage

DE_a<-as.matrix(DE_a$Gene_name)
DE_a

getwd()
venn.diagram(x = list(DE_l,DE_p, DE_a),
             category.names = c("Larva" , "Pupa", "Adult"),
             filename = '#sex_biased_genes_drosophila.png',
             output = TRUE ,
             imagetype="png" ,
             height = 480 , 
             width = 480 , 
             resolution = 300,
             compression = "lzw",
             lwd = 1,
             col=c('#fde725ff',"#440154ff", '#21908dff'),
             fill = c(alpha('#fde725ff',0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.3,
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             cat.col = c('#fde725ff',"#440154ff", '#21908dff'),
             rotation = 1)


x_venn = list(DE_l,DE_p, DE_a)
str(x_venn)

#Perform the super exact test. n=13223 (all expressed genes)
sex<- supertest(x_venn,  n=13223)
sex$overlap.sizes
sex$P.value
sex_sum<-summary(sex)
#with elements (genes from the intersections)
sex_venn_p_value<-sex_sum$Table
