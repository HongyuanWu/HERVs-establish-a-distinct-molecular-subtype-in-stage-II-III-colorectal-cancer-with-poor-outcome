#### CRC-CMS  ###
### dependencies: run if not already installed
### limma has lof of dependencies - takes time
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("Biobase", "limma"))
# install.packages("devtools")
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=rownames(TPM),
  mart=mart)

devtools::install_github("Lothelab/CMScaller")
library(Biobase)
library(CMScaller)
par(mfrow=c(1,2))

### CMS prediction of TCGA primary colorectal cancers
#res <- CMScaller(crcTCGAsubset, RNAseq=TRUE, doPlot=TRUE)
#head(res)
counts<-read.delim('Box/CRC-iMM-nonMet/expression/counts.txt')
rownames(counts)<-rownames(TPM)
res <- CMScaller(counts, RNAseq=TRUE, doPlot=TRUE, rowNames = 'ensg')
head(res)

### Camera Gene Set Analysis with CMS informative gene sets
symbol<-read.csv('Box/CRC-iMM-nonMet/expression/TPM.tsv',stringsAsFactors = FALSE, sep = '\t')$SYMBOL
rownames(counts)[which(genes$ensembl_gene_id%in%rownames(TPM))]<-gsub('X','',make.names(genes$entrezgene_id,unique = TRUE))
#new.rownames<-gsub('X','',make.names(symbol,unique = TRUE))
#new.rownames[which(duplicated(new.rownames))]<-paste0(new.rownames[which(duplicated(new.rownames))],'1')
#new.rownames[which(duplicated(new.rownames))]<-paste0(new.rownames[which(duplicated(new.rownames))],'2')
#rownames(counts)<-new.rownames

cam <- CMSgsa(emat=counts, class=res$prediction, RNAseq=TRUE)

head(cam$CMS4)

### limma differential gene expression analysis and visualization
deg <- subDEG(emat=counts, class=res$prediction, doVoom=TRUE)
subVolcano(deg, geneID="symbol")


##### HERV
CMS1.P<-CMS2.P<-CMS3.P<-CMS4.P<-rep(0,dim(hERV)[1])
M<-data.matrix(hERV)
for (i in seq(1:dim(hERV)[1])){
  CMS1.P[i]<-p.adjust(wilcox.test(M[i,which(res$prediction=='CMS1')],M[i,-which(res$prediction=='CMS1')],alternative = 'greater')$p.value)
  CMS2.P[i]<-p.adjust(wilcox.test(M[i,which(res$prediction=='CMS2')],M[i,-which(res$prediction=='CMS2')],alternative = 'greater')$p.value)
  CMS3.P[i]<-p.adjust(wilcox.test(M[i,which(res$prediction=='CMS3')],M[i,-which(res$prediction=='CMS3')],alternative = 'greater')$p.value)
  CMS4.P[i]<-p.adjust(wilcox.test(M[i,which(res$prediction=='CMS4')],M[i,-which(res$prediction=='CMS4')],alternative = 'greater')$p.value)
}
my_palette <- colorRampPalette(c("yellow","purple"))(n = 1000)
heatmap.2(M[c(which(CMS1.P<0.01),which(CMS2.P<0.01),which(CMS3.P<0.01),which(CMS4.P<0.01)), order(res$prediction)],
        col = my_palette, cexCol = 0.5, cexRow = 0.3, Rowv = FALSE, Colv = FALSE, trace = "none", scale = "row")

##### SNVs #####
SNV<-read.delim('Box/CRC-iMM-nonMet/SNV/summary/summary.tsv', header = FALSE)
rownames(SNV)<-gsub(".somatic.pass.vcf","",SNV$V1)
SNV<-SNV[,-c(1,2)]
mut<-read.delim('Box/CRC-iMM-nonMet/SNV/summary/gene', header = FALSE)
colnames(SNV)<-gsub('.txt','',mut$V1)
SNV[SNV>0]<-1
#### run CRC.R to get combined dataframe first to filter for some samples
SNV<-SNV[rownames(SNV)%in%combined$Sample,]
barplot(100*colMeans(SNV)[order(colMeans(SNV), decreasing = T)], las = 2, ylim = c(0,100), ylab = 'freq (%)')
my_palette <- colorRampPalette(c("yellow","purple"))(n = 1000)
heatmap.2(t(data.matrix(SNV))[,order(res$prediction)], col = my_palette, cexRow = 0.8, Rowv = FALSE, Colv = FALSE, trace = "none")

#### CNV ####
arm.del.event<-ifelse(arm.del>0.2,1,0)
arm.amp.event<-ifelse(arm.amp>0.2,1,0)

total.arm.del<-colSums(arm.del.event)
total.arm.amp<-colSums(arm.amp.event)

CNV.events<-as.data.frame(as.numeric(total.arm.del+total.arm.amp))
colnames(CNV.events)<-'CIN'
CNV.events$arm.amplification<-as.numeric(total.arm.amp)
CNV.events$arm.deletion<-as.numeric(total.arm.del)
rownames(CNV.events)<-rownames(SNV)

my_palette <- colorRampPalette(c("yellow","purple"))(n = 1000)
heatmap.2(t(data.matrix(CNV.events))[,order(res$prediction)], col = my_palette, cexCol = 0.5, cexRow = 1.2, Rowv = FALSE, Colv = FALSE, trace = "none")

##### MSI #####
MSI.num<-ifelse(combined$MSI=='MSH',1,0)
heatmap.2(rbind(MSI.num[order(res$prediction)],MSI.num[order(res$prediction)]), col = my_palette, cexCol = 0.5, cexRow = 1.2, Rowv = FALSE, Colv = FALSE, trace = "none")

###### Lynch Syndrome ###
LS<-read.csv('Box/CRC-iMM-nonMet/LS.csv')
LS<-LS[!LS$Lynch.Syndrome.Status=='',]
LS.filter<-LS[LS$Sample.ID%in%rownames(res),]
for(i in 1:length(LS.filter$Sample.ID)){LS.filter$CMS[i]<-res$prediction[which(rownames(res)==LS.filter$Sample.ID[i])]}
LS.filter$Lynch.Syndrome.Status[order(LS.filter$CMS)]
mean(LS.filter$Lynch.Syndrome.Status[which(LS.filter$CMS==4)]=='Negative')

prop.test(c(sum(LS.filter$Lynch.Syndrome.Status[which(LS.filter$CMS==1)]=='Negative'),
            sum(LS.filter$Lynch.Syndrome.Status[-which(LS.filter$CMS==1)]=='Negative')),
          c(length(LS.filter$Lynch.Syndrome.Status[which(LS.filter$CMS==1)]),
            length(LS.filter$Lynch.Syndrome.Status[-which(LS.filter$CMS==1)])))
y<-ifelse(LS.filter$Lynch.Syndrome.Status=='Negative',0,1)
CMS12<-as.numeric(!LS.filter$Lynch.Syndrome.Status[LS.filter$CMS==1 | LS.filter$CMS==2]=='Negative')
CMS34<-as.numeric(!LS.filter$Lynch.Syndrome.Status[LS.filter$CMS==3 | LS.filter$CMS==4]=='Negative')
wilcox.test(CMS12,CMS34)
LS.CMS<-as.data.frame(t(as.character(res$prediction)))
colnames(LS.CMS)<-rownames(res)
for(i in 1:113){LS.CMS[i]<-LS.filter$Lynch.Syndrome.Status[LS.filter$Sample.ID==colnames(LS.CMS)[i]]}
LS.CMS<-as.numeric(ifelse(LS.CMS=='Negative',0,1))
heatmap.2(rbind(LS.CMS[order(res$prediction)],LS.CMS[order(res$prediction)]), col = my_palette, cexCol = 0.5, cexRow = 1.2, Rowv = FALSE, Colv = FALSE, trace = "none")


#### HERV ######
#HERV.num<-ifelse(clin.fac$median.hERV=='high',1,0)
M<-t(clin.fac[order(res$prediction),c(26,27,30,29)])
M<-ifelse(M=='high',1,0)

heatmap.2(M, col = my_palette, cexCol = 0.5, cexRow = 1.2, Rowv = FALSE, Colv = FALSE, trace = "none")

#### PC ####
require(plot.matrix)
M<-clin.fac
rownames(M)<-combined$Sample
#M[M=='high']<-as.numeric(1); M[M=='low']<-0; M[M=='right']<-1; M[M=='left']<-0; M[M=='other']<-0; M[M=='II']<-0; M[M=='III']<-0; M[M=='MSS']<-0; M[M=='MSH']<-1
#as.factor(M)
#heatmap.2(as.matrix(M[-1,]), col = my_palette, cexCol = 0.5, cexRow = 1.2, Rowv = FALSE, Colv = FALSE, trace = "none")
M$Stage<-ifelse(M$Stage=='II',1,2)
M$Sex<-ifelse(M$Sex=='M',1,2)
M$Matastasis<-ifelse(M$Matastasis==0,1,2)
M$sidedness<-ifelse(M$sidedness=='right',0,ifelse(M$sidedness=='left',1,2))
heatmap.2(data.matrix(M[order(res$prediction),-c(1:10,14,17,20:24)]), na.rm = T, col = my_palette, cexCol = 0.78, cexRow = 0.5, Rowv = FALSE, Colv = FALSE, trace = "none")

#### fusions ###
## fusion processer 
require(dplyr)
file_list <- list.files(path="Box/CRC-iMM-nonMet/Fusions/filter/")
fusions<-list()
fusion.genes<-list()
fusion.summary <- data.frame()
fusion.common <- data.frame()
for (i in seq(length(file_list))) {
  fusions[[i]] <- read.csv(paste0("Box/CRC-iMM-nonMet/Fusions/filter/", file_list[i]), header = FALSE)
  fusions[[i]] <- fusions[[i]][which(fusions[[i]][,20]=="True"),]
  fusion.genes[[i]] <- fusions[[i]][,2:3]
  if (i == 1){fusion.common <- fusion.genes[[1]] }
}
for(i in 1:136){print(i);print(grep('NTRK1',fusion.genes[[i]]))}
file_list[61] ## CR592

for(i in 1:136){print(i);print(grep('ALK',fusion.genes[[i]]))}
file_list[15] ## CR340

for(i in 1:136){print(i);print(grep('BRAF',fusion.genes[[i]]))}
file_list[104] ## CR735

for(i in 1:136){print(i); print(grep('RSPO2',fusion.genes[[i]]))}
file_list[130] ## CR827

for(i in 1:136){print(i); print(grep('RSPO2',fusion.genes[[i]]))}
file_list[130] ## CR827

for(i in 1:136){print(i); print(grep('CCDC7',fusion.genes[[i]]))}
file_list[130] ## CR827


v<-c("ERAS","RSPO2","RSPO3","NCOA2","TCF7L2","TCF7L2","SARAF","GTF3A",
     "NAGLU","RNF121","NFATC3","USP9X","EIF3E","PTPRK","LACTB2","VTI1A","RP11","TMEM66","CDK8","IKZF3","FOLR2","PLA2G15")
for( i in 1:length(v)){print(v[i]); print(grep(v[i],unlist(fusion.genes)))}



#### some stat

prop.test(c(sum(na.rm = T, clin.fac$Age[res$prediction=='CMS1']=='high'),
            sum(na.rm = T, clin.fac$Age[!res$prediction=='CMS1']=='high')),
          c(length(clin.fac$Age[res$prediction=='CMS1']),
            length(clin.fac$Age[!res$prediction=='CMS1'])))

prop.test(c(sum(na.rm = T, clin.fac$IFNG[res$prediction=='CMS1']=='high'),
            sum(na.rm = T, clin.fac$IFNG[!res$prediction=='CMS1']=='high')),
          c(length(clin.fac$IFNG[res$prediction=='CMS1']),
            length(clin.fac$IFNG[!res$prediction=='CMS1'])))

prop.test(c(sum(na.rm = T, clin.fac$MSI[res$prediction=='CMS1']=='MSH'),
            sum(na.rm = T, clin.fac$MSI[!res$prediction=='CMS1']=='MSH')),
          c(length(clin.fac$MSI[res$prediction=='CMS1']),
            length(clin.fac$MSI[!res$prediction=='CMS1'])))

prop.test(c(sum(na.rm = T, clin.fac$CD19[res$prediction=='CMS1']=='high'),
            sum(na.rm = T, clin.fac$CD19[!res$prediction=='CMS1']=='high')),
          c(length(clin.fac$CD19[res$prediction=='CMS1']),
            length(clin.fac$CD19[!res$prediction=='CMS1'])))

prop.test(c(sum(na.rm = T, clin.fac$MATH.score[res$prediction=='CMS2']=='high'),
            sum(na.rm = T, clin.fac$MATH.score[!res$prediction=='CMS2']=='high')),
          c(length(clin.fac$MATH.score[res$prediction=='CMS2']),
            length(clin.fac$MATH.score[!res$prediction=='CMS2'])))

prop.test(c(sum(na.rm = T, clin.fac$ploidy[res$prediction=='CMS2']=='high'),
            sum(na.rm = T, clin.fac$ploidy[!res$prediction=='CMS2']=='high')),
          c(length(clin.fac$ploidy[res$prediction=='CMS2']),
            length(clin.fac$ploidy[!res$prediction=='CMS2'])))

prop.test(c(sum(na.rm = T, clin.fac$MSI[res$prediction=='CMS2']=='MSS'),
            sum(na.rm = T, clin.fac$MSI[!res$prediction=='CMS2']=='MSS')),
          c(length(clin.fac$MSI[res$prediction=='CMS2']),
            length(clin.fac$MSI[!res$prediction=='CMS2'])))

prop.test(c(sum(na.rm = T, clin.fac$CD4[res$prediction=='CMS3']=='high'),
            sum(na.rm = T, clin.fac$CD4[!res$prediction=='CMS3']=='high')),
          c(length(clin.fac$CD4[res$prediction=='CMS3']),
            length(clin.fac$CD4[!res$prediction=='CMS3'])))

prop.test(c(sum(na.rm = T, clin.fac$PD1[res$prediction=='CMS4']=='low'),
            sum(na.rm = T, clin.fac$PD1[!res$prediction=='CMS4']=='low')),
          c(length(clin.fac$PD1[res$prediction=='CMS4']),
            length(clin.fac$PD1[!res$prediction=='CMS4'])))


prop.test(c(sum(na.rm = T, clin.fac$Treg[res$prediction=='CMS4']=='low'),
            sum(na.rm = T, clin.fac$Treg[!res$prediction=='CMS4']=='low')),
          c(length(clin.fac$Treg[res$prediction=='CMS4']),
            length(clin.fac$Treg[!res$prediction=='CMS4'])))

##

prop.test(c(sum(na.rm = T, clin.fac$sidedness[res$prediction=='CMS2']=='left'),
            sum(na.rm = T, clin.fac$sidedness[!res$prediction=='CMS2']=='left')),
          c(length(clin.fac$sidedness[res$prediction=='CMS2']),
            length(clin.fac$sidedness[!res$prediction=='CMS2'])))


M<-clin.fac
M<-(M[,-c(1:10,14,17,19,20:24)])
P<-matrix(1,4,17)
prop<-P
for(i in 1:17){
  k<-1
  for (j in c('CMS1','CMS2','CMS3','CMS4')) {
  prop[k,i]<-prop.test(c(sum(na.rm = T, M[res$prediction==j,i]==labels(table(M[,i]))[[1]][1]),
              sum(na.rm = T, M[!res$prediction==j,i]==labels(table(M[,i]))[[1]][1])),
            c(length(M[res$prediction==j,i]),
              length(M[!res$prediction==j,i])))$estimate[1]
  P[k,i]<-prop.test(c(sum(na.rm = T, M[res$prediction==j,i]==labels(table(M[,i]))[[1]][1]),
                      sum(na.rm = T, M[!res$prediction==j,i]==labels(table(M[,i]))[[1]][1])),
                    c(length(M[res$prediction==j,i]),
                      length(M[!res$prediction==j,i])))$p.value
  k<-k+1
  }
}
rownames(prop)<-c('CMS1','CMS2','CMS3','CMS4')
colnames(P)<-colnames(prop)<-colnames(M[,-18])
for(i in 1:17){P[,i]<-p.adjust(P[,i])}
#n<-ifelse(P<0.1,ifelse(P<0.05,ifelse(P<0.01,ifelse(P<0.001,'****','***'),'**'),'*'),'')
heatmap.2(prop, na.rm = T, col = my_palette, cexCol = 0.78, cexRow = 1, Rowv = FALSE, Colv = FALSE, trace = "none")


###### deletion and CMS ###
gene.del<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.gene.amp.tsv",stringsAsFactors = FALSE)
gene.del<-gene.del[,rownames(res)]

num<-dim(gene.del)[1]
P<-matrix(1,4,num)
prop<-P
M<-t(gene.del)
for(i in 1:num){
  k<-1
  for (j in c('CMS1','CMS2','CMS3','CMS4')) {
    prop[k,i]<-prop.test(c(sum(na.rm = T, M[res$prediction==j,i]==labels(table(M[,i]))[[1]][1]),
                           sum(na.rm = T, M[!res$prediction==j,i]==labels(table(M[,i]))[[1]][1])),
                         c(length(M[res$prediction==j,i]),
                           length(M[!res$prediction==j,i])))$estimate[1]
    P[k,i]<-prop.test(c(sum(na.rm = T, M[res$prediction==j,i]==labels(table(M[,i]))[[1]][1]),
                        sum(na.rm = T, M[!res$prediction==j,i]==labels(table(M[,i]))[[1]][1])),
                      c(length(M[res$prediction==j,i]),
                        length(M[!res$prediction==j,i])))$p.value
    k<-k+1
  }
}

rownames(prop)<-c('CMS1','CMS2','CMS3','CMS4')
heatmap.2(prop[,colMins(P)<0.00001], na.rm = T, col = my_palette, cexCol = 0.78, cexRow = 1, Rowv = FALSE, Colv = FALSE, trace = "none")




##### CMS survival ####
clin.fac$CMS<-res$prediction
clin.fac.1<-clin.fac[clin.fac$median.hERV=='low',]
clin.fac.1 %>% analyse_survival(vars(OS, Dead), CMS) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     risk.table.height = 0.15,
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac.1 %>%
       analyse_survival(vars(RFS, relapse), 
                        by=CMS)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("CMS"),
                                          legend=c('right','right')
                    )) %>%
  print



clin.fac.1 %>% analyse_multivariate(vars(OS, Dead),
                                  covariates = vars(CMS)) -> result
#print(result)

forest_plot(result)



##### CMS5
clin.fac$CMS<-res$prediction
clin.fac.1<-clin.fac
clin.fac.1$CMS<-as.character(clin.fac.1$CMS)
clin.fac.1$CMS[clin.fac$median.hERV=='high']<-'CMS/HERV-high'
clin.fac.1$CMS<-as.factor(clin.fac.1$CMS)
clin.fac.1 %>% analyse_survival(vars(OS, Dead), CMS) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     risk.table.height = 0.15,
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac.1 %>%
       analyse_survival(vars(RFS, relapse), 
                        by=CMS)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("CMS"),
                                          legend=c('right','right')
                    )) %>%
  print




##### LS survival
clin.fac.2<-clin.fac
clin.fac.2$LS<-LS.CMS
clin.fac.2<-clin.fac.2[clin.fac.2$CMS=='CMS3' | clin.fac.2$CMS=='CMS4',]

clin.fac.2$LS<-ifelse(clin.fac.2$LS==1,'positive','negative')

clin.fac.2 %>% analyse_survival(vars(OS, Dead), LS) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     risk.table.height = 0.15,
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac.2 %>%
       analyse_survival(vars(RFS, relapse), 
                        by=LS)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("Lynch Syndrome"),
                                          legend=c('right','right')
                    )) %>%
  print


