library('DESeq2')
library('pheatmap')
library('RColorBrewer')
library('ggplot2')

setwd("~/Desktop/F3-all/")
x<-read.csv('f3_all-apr19.txt', row.names="Geneid", header=T, sep='\t')
counts<-x [which (x$Length>200),] #only select fragments larger than 200bp (there are many small ones )
counts$Length=NULL
#outlier removal
counts$A17=NULL 
counts$sA54=NULL
counts$aAB1=NULL
counts$C63=NULL
colnames(counts)

#load samples based on treatment 
counts=counts[,order(names(counts))]
countdata<-as.matrix(counts)
#countdata<-apply(counts,1,function(x) {all(x>0)})
#conditions for LRT
condition<-factor(c(rep('A',7), rep ('aAB',5), rep('aAC',7), rep('B',8),  rep('C',7), rep('LABt',6), rep('LABtdA',8), rep('LBAt',6), rep('LBAtaB',7), rep('LBCt',8), rep('LCAt',4), rep('LCAtaC',5), rep('sA',7), rep('sAaB',5), rep('sC',5), rep('tSC',8), rep('ttB',5), rep('ttBaC',4))) #each treatment independently 
coldata <-data.frame(row.names=colnames(countdata), condition)

#double check data read counts
totalCounts=colSums(countdata) 
totalCounts
#A108      A109      A110      A117      A126      A136       A18      aAB2      aAB3 
#28757171  32862723  38720434  33192232  33514041  34039815  36038993  45745785  53406682 
#aAB4      aAB5      aAB6      aAC1      aAC2      aAC3      aAC4      aAC5      aAC6 
#55346741  48312959  51715784  42266547  45248288  38515538  47114677  52816995  49194908 
#aAC7      B127      B133      B213      B214      B215      B216      B217      B219 
#46230932  29336568  29719331  33131939  33858863  34920415  39511613  32590437  41023536 
#C104      C105      C107       C15      C156      C159       C62   LABt119   LABt121 
#34103757  36214539  29817356  31395507  30916693  32575268  30991499  28925536  25811459 
#LABt128   LABt130   LABt134   LABt135 LABtdA111 LABtdA124  LABtdA13  LABtdA14 LABtdA141 
#33227005  33545116  27110160  37490761  36706076  34847195  28862646  30413674  34675251 
#LABtdA142 LABtdA143 LABtdA144   LBAt113   LBAt137   LBAt138   LBAt139   LBAt140   LBAt160 
#33182110  31024861  32345754  35800479  33771891  28554585  31623021  38160315  31404330 
#LBAtaB1   LBAtaB2   LBAtaB3   LBAtaB4   LBAtaB5   LBAtaB6   LBAtaB7   LBCt161   LBCt162 
#46059624  50440880  37145384  42515228  40731529  49084239  45678240  44964392  38789834 
#LBCt163   LBCt164   LBCt165   LBCt166   LBCt167   LBCt174   LCAt122   LCAt125    LCAt29 
#48197869  44082223  36641148  37909592  30669340  40300967  35495368  41823153  29960958 
#LCAt30 LCAtaC201 LCAtaC202 LCAtaC203 LCAtaC204 LCAtaC228     sA114     sA115      sA16 
#42939108  34975635  29979066  43315717  34966256  35889197  27587988  26805659  36640183 
#sA21      sA53      sA54      sA97      sA98     sAaB1     sAaB2     sAaB3     sAaB4 
#31749350  34324999  26626823  29258776  43313572  41412856  39527683  42870012  49012565 
#sAaB5     sC131     sC149     sC151     sC153     sC154     tSC19    tSC199    tSC200 
#49950712  33513892  36725826  38976543  39145400  34390740  39605686  37702845  36358274 
#tSC229     tSC23     tSC24     tSC26     tSC88    ttB118    ttB129    ttB175    ttB176 
#36109834  34305870  42090420  36422669  31781123  40094410  37918434  26824672  36496360 
#ttB177    ttBaC1    ttBaC2    ttBaC3    ttBaC4 
#43208915  48182831  41522786  51043620  50399348 

min(totalCounts) #25811459
max(totalCounts) #55346741
mean(totalCounts) #37754112

#DESEQ with all samples just to explor the data 
dds<-DESeqDataSetFromMatrix (countdata, coldata, design= ~condition)
dds<-dds[rowSums(counts(dds))>0,]
summary(dds) # 29681 
dds.c <- DESeq(dds, test="LRT",full=~condition, reduced=~1)
res.c<-results(dds.c)
res.c
res.c<-res.c[order(res.c$padj),]
table(res.c$padj<0.01)
#FALSE  TRUE 
#11065 10507 comparison just out of curiosity, as too many conditions are included to find a pattern 

#to obtain variance stabilized data folow: 
head(res.c)
vsd.c=getVarianceStabilizedData(dds.c)
head(vsd.c)
head(res.c)
vals.c=cbind(res.c$pvalue, res.c$padj)
head(vals.c)
colnames(vals.c)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.c,vals.c)
head(vsdpvals.c)
write.csv(vsdpvals.c, "F3-all_LRT_VSDandPVALS_apr-19.csv", quote=F)

#sample to sample matrix to determine outliers 
sampleDists <- dist(t(vsd.c))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd.c))
colnames(sampleDistMatrix) <- paste(colnames(vsd.c))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(250)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

###PCOA all samples - This won't show much, as it's just a clump of data, it's just exploratory 
library('vegan')
library('rgl')
library('ape')

dat=read.csv("F3-all_LRT_VSDandPVALS_apr-19.csv")
dat=dat[which(dat$padj.c<0.01),]
nrow(dat)#10868
names(dat)
data=dat[,2:114]
row.names(data)=dat$X
head(data)

# extracting experimental conditions to include colors in the PCA 
traits=read.csv("f3-traits-all.csv",header=T)
conditions=c('black','black','black','black','black','black','black','darkorange','darkorange','darkorange','darkorange','darkorange','darkorchid1','darkorchid1','darkorchid1','darkorchid1','darkorchid1','darkorchid1','darkorchid1','dodgerblue','dodgerblue','dodgerblue','dodgerblue','dodgerblue','dodgerblue','dodgerblue','dodgerblue','firebrick1','firebrick1','firebrick1','firebrick1','firebrick1','firebrick1','firebrick1','gray55','gray55','gray55','gray55','gray55','gray55','green4','green4','green4','green4','green4','green4','green4','green4','lightcyan2','lightcyan2','lightcyan2','lightcyan2','lightcyan2','lightcyan2','yellow2','yellow2','yellow2','yellow2','yellow2','yellow2','yellow2','springgreen1','springgreen1','springgreen1','springgreen1','springgreen1','springgreen1','springgreen1','springgreen1','azure3','azure3','azure3','azure3','hotpink','hotpink','hotpink','hotpink','hotpink','coral3','coral3','coral3','coral3','coral3','coral3','coral3','coral3','darkgoldenrod','darkgoldenrod','darkgoldenrod','darkgoldenrod','darkgoldenrod','brown4','brown4','brown4','brown4','brown4','chocolate','chocolate','chocolate','chocolate','chocolate','chocolate','chocolate','chocolate','slategray','slategray','slategray','slategray','slategray','cornflowerblue','cornflowerblue','cornflowerblue','cornflowerblue')
# principal coordinate calculation
dd.pcoa=pcoa(dist(t(data)))
scores=dd.pcoa$vectors
dd.pcoa 

# plotting the first two principal axes
quartz() 
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(scores[,1], scores[,2],col=conditions, pch=19, cex=2, main= "PCoA for all treatments", xlab='PC1=22.74%', ylab='PC2=12.93%')
legend("right", c('A','aAB','aAC','B','C','LABt','LABtdA','LBAt','LABtaB','LBCt','LCAt','LCAtaC','sA','sAaB','sC','tSC','ttB','ttBaC'), col=c('black','darkorange','darkorchid1','dodgerblue','firebrick1','gray55','green4','lightcyan2','yellow2','springgreen1','azure3','hotpink','coral3','darkgoldenrod','brown4','chocolate','slategray','cornflowerblue'), pch=c(19),cex=0.75, horiz=FALSE)
par(mar=c(5, 5, 5, 2) + 0.1)

# plotting the second and third principal axes  
quartz()
plot(scores[,2], scores[,3],col=conditions,pch=symbols)
ordispider(dd.pcoa)

############

###Pairwise comparisons between categories: 
###ttB vs sAaB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=7) #this actually does the analysis of gene expression, replacing count outliers by the average expression value
res<-results(dds, contrast=c('condition', 'ttB', 'sAaB'))  #here is where the two contrasting conditions get defined
write.csv(res, file="ttB-sAaB-res-table4GO.csv", quote=FALSE) #this table was edited and then used for the GO-MWU
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#16190   337    
sigs<-res[which(res$padj<0.01),] #this command is only to write how many DEGs were found on the screen 
write.csv(sigs, file='ttBvssAaB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("ttB118", "ttB175", "ttB176", "ttB177", "ttB129","sAaB1", "sAaB2", "sAaB3", "sAaB4", "sAaB5")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "ttb-saab_VSDandPVALS-feb2018.csv", quote=F)

###########
###ttB vs A
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttB', 'A')) 
write.csv(res, file="ttB-A-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18704    99 
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttBvsA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("ttB118", "ttB175", "ttB176", "ttB177", "ttB129","A108", "A109", "A110", "A117", "A126", "A136", "A18")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "ttb-A_VSDandPVALS-mar2018.csv", quote=F)


#######
###A vs B
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'B'))
write.csv(res, file="A-B-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#17930   333
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='AvsB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","B127", "B213", "B214", "B215", "B216", "B217", "B219","B133")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-B_VSDandPVALS-mar2018.csv", quote=F)

#######
###A vs C
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'C')) 
write.csv(res, file="A-C-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18870   502 
sigs<-res[which(res$padj<0.01),]
write.csv(sigs, file='AvsC-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","C104","C105","C107","C15","C156","C159","C62")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-C_VSDandPVALS-mar2018.csv", quote=F)

#####
###A vs LABtdA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'LABtdA'))
write.csv(res, file="A-LABtdA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#15509    30   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-LABtdA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","LABtdA111", "LABtdA1.24","LABtdA13", "LABtdA141", "LABtdA142", "LABtdA143", "LABtdA144","LABtdA14")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-LABtdA_VSDandPVALS-mar2018.csv", quote=F)

#####
###LBCt vs ttBaC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LBCt', 'ttBaC'))
write.csv(res, file="LBCt-ttBaC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19240   696   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LBCt-ttBaC-sig.csv')

#####
###A vs LABt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'LABt')) 
write.csv(res, file="A-LABt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#False True
#22705    66 
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-LABt-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-LABt_VSDandPVALS-mar2018.csv", quote=F)

#####
###B vs LABt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'B', 'LABt')) 
write.csv(res, file="B-LABt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18585   268   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='B-LABt-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("B127", "B213", "B214", "B215", "B216", "B217", "B219","B133","LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "B-LABt_VSDandPVALS-mar2018.csv", quote=F)

#####
###B vs ttB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'B', 'ttB')) 
write.csv(res, file="B-ttB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  19105   297 
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='B-ttB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("B127", "B213", "B214", "B215", "B216", "B217", "B219","B133","ttB118", "ttB175", "ttB176", "ttB177", "ttB129")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "B-ttB_VSDandPVALS-mar2018.csv", quote=F)

#####
###LABt vs ttB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LABt', 'ttB')) 
write.csv(res, file="LABt-ttB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
# 1667    14 
sigs<-res[which(res$padj<0.01),] 
#write.csv(sigs, file='LABt-ttB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121","ttB118", "ttB175", "ttB176", "ttB177", "ttB129")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LABt-ttB_VSDandPVALS-mar2018.csv", quote=F)

#####
###A vs LBAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'LBAt')) 
write.csv(res, file="A-LBAt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#21563    30  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-LBAt-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-LBAt_VSDandPVALS-mar2018.csv", quote=F)

#####
###A vs aAB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'aAB')) 
write.csv(res, file="A-aAB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18822  1117   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-aAB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","aAB2", "aAB3", "aAB4", "aAB5", "aAB6")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-aAB_VSDandPVALS-mar2018.csv", quote=F)

#####
###B vs aAB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'B', 'aAB')) 
write.csv(res, file="B-aAB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#17813  2692  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='B-aAB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("B127", "B213", "B214", "B215", "B216", "B217", "B219","B133", "aAB3", "aAB4", "aAB5", "aAB6")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "B-aAB_VSDandPVALS-mar2018.csv", quote=F)


#####
###A vs sA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'sA')) 
write.csv(res, file="A-sA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19182   220  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-sA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","sA16","sA21", "sA53", "sA97", "sA98", "sA114", "sA115")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-sA_VSDandPVALS-mar2018.csv", quote=F)

#####
###sA vs LABtdA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LABtdA', 'sA')) 
write.csv(res, file="LABtdA-sA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18407   446 
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LABtda-sA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LABtdA111", "LABtdA1.24","LABtdA13", "LABtdA141", "LABtdA142", "LABtdA143", "LABtdA144","LABtdA14","sA16","sA21", "sA53", "sA97", "sA98", "sA114", "sA115")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LABtda-sA_VSDandPVALS-mar2018.csv", quote=F)

###################################
###LABt vs LBAtaB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LBAtaB', 'LABt')) 
write.csv(res, file="LABt-LBAtaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#17252  2150   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LABt-LBAtaB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LBAtaB1", "LBAtaB2", "LBAtaB3", "LBAtaB4", "LBAtaB5", "LBAtaB6", "LBAtaB7", "LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LABt-LBAtaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###LABtadA vs LBAtaB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LBAtaB', 'LABtdA')) 
write.csv(res, file="LABtadA-LBAtaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#16525   672   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LABtdA-LBAtaB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LBAtaB1", "LBAtaB2", "LBAtaB3", "LBAtaB4", "LBAtaB5", "LBAtaB6", "LBAtaB7","LABtdA111", "LABtdA1.24","LABtdA13", "LABtdA141", "LABtdA142", "LABtdA143", "LABtdA144","LABtdA14")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LABtdA-LBAtaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###LBAt vs LABtdA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LBAt', 'LABtdA')) 
write.csv(res, file="LBAt-LABtdA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#28734     1   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LBAt-LABtdA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LBAt137", "LBAt138", "LBAt139", "LBAt113","LBAt140", "LBAt160","LABtdA111", "LABtdA1.24","LABtdA13", "LABtdA141", "LABtdA142", "LABtdA143", "LABtdA144","LABtdA14")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LBAt-LBAtaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###B vs LABtdA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'B', 'LABtdA')) 
write.csv(res, file="B-LABtdA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#16768  2085     
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='B-LABtdA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("B127", "B213", "B214", "B215", "B216", "B217", "B219","B133","LABtdA111", "LABtdA1.24","LABtdA13", "LABtdA141", "LABtdA142", "LABtdA143", "LABtdA144","LABtdA14")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "B-LBAtaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###B vs sA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'B', 'sA')) 
write.csv(res, file="B-sA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#17321   428    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='B-sA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("B127", "B213", "B214", "B215", "B216", "B217", "B219","B133","sA16","sA21", "sA53", "sA97", "sA98", "sA114", "sA115")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "B-sA_VSDandPVALS-mar2018.csv", quote=F)

###################################
###B vs LBAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'B', 'LBAt')) 
write.csv(res, file="B-LBAt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18500   902     
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='B-LBAt-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("B127", "B213", "B214", "B215", "B216", "B217", "B219","B133","LBAt137", "LBAt138", "LBAt139", "LBAt113","LBAt140", "LBAt160")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "B-LBAt_VSDandPVALS-mar2018.csv", quote=F)

###################################
###B vs LBAtaB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'B', 'LBAtaB')) 
write.csv(res, file="B-LBAtaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#16297  2556   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='B-LBAtaB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("B127", "B213", "B214", "B215", "B216", "B217", "B219","B133","LBAtaB1", "LBAtaB2", "LBAtaB3", "LBAtaB4", "LBAtaB5", "LBAtaB6", "LBAtaB7")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "B-LBAtaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###LABt vs sA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LABt', 'sA')) 
write.csv(res, file="LABt-sA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#22103    35   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LABt-sA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121","sA16","sA21", "sA53", "sA97", "sA98", "sA114", "sA115")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LABt-sA_VSDandPVALS-mar2018.csv", quote=F)

###################################
###LABt vs sAaB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LABt', 'sAaB')) 
write.csv(res, file="LABt-sAaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#17935   368    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LABt-sAaB-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121","sAaB1", "sAaB2", "sAaB3", "sAaB4", "sAaB5")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LABt-sAaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###LABt vs LABtdA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LABt', 'LABtdA')) 
write.csv(res, file="LABt-LABtdA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#14336   271    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LABt-LABtdA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121","LABtdA111", "LABtdA124","LABtdA13", "LABtdA141", "LABtdA142", "LABtdA143", "LABtdA144","LABtdA14")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LABt-LABtdA_VSDandPVALS-mar2018.csv", quote=F)

###################################
###LABt vs LBAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LABt', 'LBAt')) 
write.csv(res, file="LABt-sAaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19369    33   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LABt-LBAt-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LABt128", "LABt130","LABt134", "LABt135", "LABt119", "LABt121","LBAt137", "LBAt138", "LBAt139", "LBAt113","LBAt140", "LBAt160")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LABt-LBAt_VSDandPVALS-mar2018.csv", quote=F)

###################################
###ttB vs LABtdA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttB', 'LABtdA')) 
write.csv(res, file="ttB-LABtdA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#12742   587    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttB-LABtdA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("ttB118", "ttB175", "ttB176", "ttB177", "ttB129", "LABtdA111", "LABtdA1.24","LABtdA13", "LABtdA141", "LABtdA142", "LABtdA143", "LABtdA144","LABtdA14")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "ttB-LABtdA_VSDandPVALS-mar2018.csv", quote=F)

###################################
###sAaB vs LABtdA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sAaB', 'LABtdA')) 
write.csv(res, file="sAaB-LABtdA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#16039   606    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sAaB-LABtdA-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("sAaB1", "sAaB2", "sAaB3", "sAaB4", "sAaB5", "LABtdA111", "LABtdA1.24","LABtdA13", "LABtdA141", "LABtdA142", "LABtdA143", "LABtdA144","LABtdA14")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "sAaB-LABtdA_VSDandPVALS-mar2018.csv", quote=F)

###################################
###ttB vs sA
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttB', 'sA')) 
write.csv(res, file="ttB-sA-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18816    37   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttB-sA-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("ttB118", "ttB175", "ttB176", "ttB177", "ttB129","sA16","sA21", "sA53", "sA97", "sA98", "sA114", "sA115")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "ttB-sA_VSDandPVALS-mar2018.csv", quote=F)

###################################
###ttB vs LBAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttB', 'LBAt')) 
write.csv(res, file="ttB-LBAt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#6672    18    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttB-LBAt-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("ttB118", "ttB175", "ttB176", "ttB177", "ttB129","LBAt137", "LBAt138", "LBAt139", "LBAt113","LBAt140", "LBAt160")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "ttB-LBAt_VSDandPVALS-mar2018.csv", quote=F)

###################################
###ttB vs tsC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttB', 'tSC')) 
write.csv(res, file="ttB-tsC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#16862   233    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttB-tSC-sig.csv') 

###################################
###ttB vs sC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttB', 'sC')) 
write.csv(res, file="ttB-sC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19250   122     
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttB-sC-sig.csv') 

###################################
###ttB vs LBAtaB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttB', 'LBAtaB')) 
write.csv(res, file="ttB-LBAtaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#15209  1988     
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttB-LBAtaB-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("ttB118", "ttB175", "ttB176", "ttB177", "ttB129","LBAtaB1", "LBAtaB2", "LBAtaB3", "LBAtaB4", "LBAtaB5", "LBAtaB6", "LBAtaB7")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "ttB-LBAtaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###sA vs LBAtaB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sA', 'LBAtaB')) 
write.csv(res, file="sA-LBAtaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#15924  1273    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sA-LBAtaB-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("sA16","sA21", "sA53", "sA97", "sA98", "sA114", "sA115","LBAtaB1", "LBAtaB2", "LBAtaB3", "LBAtaB4", "LBAtaB5", "LBAtaB6", "LBAtaB7")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "sA-LBAtaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###sA vs LBAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sA', 'LBAt')) 
write.csv(res, file="sA-LBAt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19085   317    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sA-LBAt-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("sA16","sA21", "sA53", "sA97", "sA98", "sA114", "sA115","LBAt137", "LBAt138", "LBAt139", "LBAt113","LBAt140", "LBAt160")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "sA-LBAt_VSDandPVALS-mar2018.csv", quote=F)

###################################
###A vs LBAtaB
#dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
#dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'LBAtaB')) 
write.csv(res, file="A-LBAtaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18354  1048    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-LBAtaB-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","LBAtaB1", "LBAtaB2", "LBAtaB3", "LBAtaB4", "LBAtaB5", "LBAtaB6", "LBAtaB7")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-LBAtaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###A vs sAaB
#dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
#dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'sAaB')) 
write.csv(res, file="A-sAaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18413   440     
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-sAaB-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("A108", "A109", "A110", "A117", "A126", "A136", "A18","sAaB1", "sAaB2", "sAaB3", "sAaB4", "sAaB5")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "A-sAaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###sAaB vs LBAtaB
#dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
#dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LBAtaB', 'sAaB')) 
write.csv(res, file="LBAtaB-sAaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18845     8      
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LBAtaB-sAaB-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LBAtaB1", "LBAtaB2", "LBAtaB3", "LBAtaB4", "LBAtaB5", "LBAtaB6", "LBAtaB7","sAaB1", "sAaB2", "sAaB3", "sAaB4", "sAaB5")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LBAtaB-sAaB_VSDandPVALS-mar2018.csv", quote=F)

###################################
###sAaB vs LBAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sAaB', 'LBAt')) 
write.csv(res, file="sAaB-LBAt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#17935   368    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sAaB-LBAt-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("sAaB1", "sAaB2", "sAaB3", "sAaB4", "sAaB5", "LBAt137", "LBAt138", "LBAt139", "LBAt113","LBAt140", "LBAt160")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "sAaB-LBAt_VSDandPVALS-mar2018.csv", quote=F)

###################################
###sAaB vs B
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sAaB', 'B')) 
write.csv(res, file="sAaB-B-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#17338   965   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sAaB-B-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("sAaB1", "sAaB2", "sAaB3", "sAaB4", "sAaB5", "B127", "B213", "B214", "B215", "B216", "B217", "B219","B133")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "sAaB-B_VSDandPVALS-mar2018.csv", quote=F)

###################################
###LBAtaB vs LBAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LBAtaB', 'LBAt')) 
write.csv(res, file="LBAtaB-vs-LBAt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#16198   447    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LBAtaB-vs-LBAt-sig.csv')

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("LBAtaB1", "LBAtaB2", "LBAtaB3", "LBAtaB4", "LBAtaB5", "LBAtaB6", "LBAtaB7","LBAt137", "LBAt138", "LBAt139", "LBAt113","LBAt140", "LBAt160")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "LBAtaB-vs-LBAt_VSDandPVALS-mar2018.csv", quote=F)

###################################
###sA vs sAaB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sA', 'sAaB')) 
write.csv(res, file="sA-sAaB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#14983   556      
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sA-sAaB-sig.csv') 

#to get variance stabilized data, and the corresponding table 
vsd=getVarianceStabilizedData(dds)
vsd=as.data.frame(vsd)
vsd.sub=vsd[c("sA16","sA21", "sA53", "sA97", "sA98", "sA114", "sA115","sAaB1", "sAaB2", "sAaB3", "sAaB4", "sAaB5")]
head(vsd.sub)
head(res)
vals=cbind(res$pvalue, res$padj)
head(vals)
colnames(vals)=c("pval.c", "padj.c")
vsdpvals.c=cbind(vsd.sub,vals)
head(vsdpvals.c)
write.csv(vsdpvals.c, "sA-sAaB_VSDandPVALS-mar2018.csv", quote=F)

###PROJECT-2 F3
###C vs tSC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'C', 'tSC')) 
write.csv(res, file="C-tsC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#20179   288  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='CvstSC-sig.csv')

#######
###C vs LBCt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'C', 'LBCt')) 
write.csv(res, file="C-LBCt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19706  1322  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='CvsLBCt-sig.csv')

#######
###tsC vs LBCt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'tSC', 'LBCt')) 
write.csv(res, file="tSC-LBCt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19033  1995   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='tSCvsLBCt-sig.csv')

#######
###aAC vs C
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'aAC', 'C')) 
write.csv(res, file="aAC-C-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19158   185   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='aACvsC-sig.csv')

#######
###B vs C
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'B', 'C')) 
write.csv(res, file="B-C-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#15888    74   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='BvsC-sig.csv')

#######
###aAC vs sC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'aAC', 'sC')) 
write.csv(res, file="aAC-sC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19017   888    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='aAC-SC-sig.csv')

###ttBaC vs sC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttBaC', 'sC')) 
write.csv(res, file="ttBaC-sC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18091   124 (should be 111)
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttBaC-SC-sig.csv')

### aAC vs ttBaC 
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'aAC', 'ttBaC')) 
write.csv(res, file="aAC-ttBaC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18091   124   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='aAC-ttBaC-sig.csv')

#######
###sC vs LBCt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sC', 'LBCt')) 
write.csv(res, file="sC-LBCt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19172  1295    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sC-LBCt-sig.csv')

#######
###sC vs tsC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sC', 'tSC')) 
write.csv(res, file="sC-tSC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#20271   196   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sC-tSC-sig.csv')

#######
###aAC vs LCAtaC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'aAC', 'LCAtaC')) 
write.csv(res, file="aAC-LCAtaC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#20347   120    
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='aAC-LCAtaC-sig.csv')

#######
###A vs LCAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'LCAt')) 
write.csv(res, file="A-LCAt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#22138    24   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-LCAt-sig.csv')

#######
###C vs ttBaC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'C', 'ttBaC')) 
write.csv(res, file="C-ttBaC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='C-ttBaC-sig.csv')

#######
###C vs ttB
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'C', 'ttB')) 
write.csv(res, file="C-ttB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  		600 
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='C-ttB-sig.csv')

#######
###C vs sC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'C', 'sC')) 
write.csv(res, file="C-sC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='C-sC-sig.csv')

#######
###C vs LCAt
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'C', 'LCAt')) 
write.csv(res, file="C-LCAt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='C-LCAt-sig.csv')

#######
###C vs LCAtac
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'C', 'LCAtac')) 
write.csv(res, file="C-LCAtac-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='C-LCAtac-sig.csv')

#######
###A vs C
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'C')) 
write.csv(res, file="A-C-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#18870   502  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-C-sig.csv')

#######
###LBCt vs aAC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LBCt', 'aAC')) 
write.csv(res, file="LBCt-aAC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19206  2994   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LBCt-aAC-sig.csv')

#######
###tSC vs aAC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'tSC', 'aAC')) 
write.csv(res, file="tsC-aAC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19637   302  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='tSC-aAC-sig.csv')

#######
###ttB vs LBCt - this last one is behaving weirdly 
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'ttB', 'LBCt')) 
write.csv(res, file="ttB-LBCt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#17234  3271  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='ttB-LBCt-sig.csv')

#######
###aAB vs B  
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'aAB', 'B')) 
write.csv(res, file="aAB-B-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='aAB-B-sig.csv')

#######
###aAB vs ttB  
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'aAB', 'ttB')) 
write.csv(res, file="aAB-ttB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='aAB-ttB-sig.csv')

#######
###aAB vs ttB  
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'aAB', 'ttB')) 
write.csv(res, file="aAB-ttB-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='aAB-ttB-sig.csv')

#######
###aAC vs LBCt 
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'aAC', 'LBCt')) 
write.csv(res, file="aAC-LBCt-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='aAC-LBCt-sig.csv')

##C vs LCAtaC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'C', 'LCAtaC')) 
write.csv(res, file="C-LCAtaC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#19257   115  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='C-LCAtaC-sig.csv')

##sC vs LCAtaC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'sC', 'LCAtaC')) 
write.csv(res, file="sC-LCAtaC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#20412    93  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='sC-LCAtaC-sig.csv')

##A vs LCAtaC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'A', 'LCAtaC')) 
write.csv(res, file="A-LCAtaC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#20412    32  
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='A-LCAtaC-sig.csv')

##LCAt vs LCAtaC
dds<- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds<-DESeq(dds, minReplicatesForReplace=6) 
res<-results(dds, contrast=c('condition', 'LCAt', 'LCAtaC')) 
write.csv(res, file="LCAt-LCAtaC-res-table4GO.csv", quote=FALSE)
mcols(res,use.names=TRUE)
table(res$padj<0.01) 
#FALSE  TRUE 
#29581    47   
sigs<-res[which(res$padj<0.01),] 
write.csv(sigs, file='LCAt-LCAtaC-sig.csv')
