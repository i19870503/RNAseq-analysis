library("WGCNA", lib.loc="/home/zluna/R/x86_64-pc-linux-gnu-library/3.1")
# 预设置
options(stringsAsFactors=FALSE)  
allowWGCNAThreads()


### 模拟数据


no.obs=30 # 设置样本个数
ESturquoise<-0; ESbrown<-.6; ESgreen<-.6; ESyellow=0 # 设置基因模块的显著性水平
ESvector<-c (ESturquoise, ESbrown, ESgreen,ESyellow)
nGenes1<-3000   # 设置基因个数为 3 000
simulateProportions1 <-c (0.2, 0.15, 0.08, 0.06,0.04) # 模块中基因个数比
set.seed (1)

# 模拟基因网
MEgreen<-rnorm (no.obs);
scaledy<-MEgreen*ESgreen+sqrt (1-ESgreen^2)*rnorm (no.obs);
y<-ifelse(scaledy>median(scaledy), 1, 0)

# 模拟的表型变量

MEturquoise<-ESturquoise*scaledy+sqrt (1-ESturquoise^2)*rnorm (no.obs)
MEblue<-.6*MEturquoise+sqrt (1-.6^2)*rnorm (no.obs)
MEbrown<-ESbrown*scaledy+sqrt (1-ESbrown^2 )*rnorm(no.obs)
MEyellow<-ESyellow*scaledy+sqrt (1-ESyellow^2)*rnorm (no.obs)
MEN1<-data.frame (y, MEturquoise, MEblue, MEbrown, MEgreen, MEyellow)
dat1<-simulateDatExpr5Modules (MEturquoise=MEN1$MEturquoise,MEblue=MEN1$MEblue, MEbrown=MEN1$MEbrown, MEyellow=MEN1$MEyellow,MEgreen=MEN1$MEgreen, nGenes=nGenes1,simulateProportions=simulateProportions1)
datExpr<-dat1$datExpr;
truemodule<-dat1$truemodule;
datME<-dat1$datME;
attach (MEN1)
datExpr<-data.frame (datExpr)

# 模拟的表达谱数据

ArrayName<-paste ("Sample",1:dim (datExpr)[[1]],sep="")
GeneName <-paste ("Gene",1:dim (datExpr) [[2]],sep="")

Dimnames (datExpr)[[1]]=ArrayName  ##不能用
Dimnames (datExpr)[[2]]=GeneName   ##不能用
GeneName
datExpr
str(datExpr)

head(datExpr)
#
#    ...  Gene2523    Gene2524    Gene2525    Gene2526    Gene2527   Gene2528   Gene2529    Gene2530
#Sample1  1.0693947  2.42077209 -1.30880548 -0.82173531 -0.54665674  0.4460768  1.2946842 -0.06434581
#Sample2 -1.7211642  0.04508150 -0.09911778  0.04592958 -0.77328274 -0.6939014 -1.3745110 -1.38650817
#Sample3  0.1728964  1.08583307  1.17529063  1.71279874 -0.49480656  0.7793184  0.6097631 -0.49821162

datExpr[[1]]
rownames(datExpr)
ArrayName
rownames(datExpr)=ArrayName
colnames(datExpr)=GeneName

# 绘制样本聚类图
plotClusterTreeSamples (datExpr=datExpr, y=y)

pdf('wgcna_test.pdf')
plotClusterTreeSamples (datExpr=datExpr, y=y);
dev.off();

######选择构建网络参数
# 设置网络构建参数选择范围，计算无尺度分布拓扑矩阵

powers=c (c (1:10), seq (from=12, to=20, by=2))
sft=pickSoftThreshold (datExpr, powerVector=powers,verbose=5);

# 图形绘制
par(mfrow=c(1,2));
cex1=0.5;

plot(sft$fitIndices [,1], -sign(sft $fitIndices [,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)", ylab=" Scale Free Topology Model Fit, signed R^2",type="n",main=paste("Scale independence"), ylim=c(0,1));
text (sft$fitIndices [,1], -sign (sft$fitIndices[,3])*sft$fitIndices [,2],labels=powers,cex=cex1, col="black");
abline (h=0.90, col="black")

plot (sft $fitIndices [,1], sft $fitIndices [,5], type ="n", main=paste ("Mean connectivity"),ylab =" MeanConnectivity", xlab ="SoftThreshold(power)")
text (sft $fitIndices [,1], sft $fitIndices [,5], labels=powers, cex=cex1, col="black")

#构建网络
softPower<-2; # 设定无尺度参数为 2
adjacency= adjacency (datExpr, power=softPower);
TOM=TOMsimilarity (adjacency);
dissTOM=1-TOM
# 聚类分析
geneTree=flashClust (as.dist (dissTOM), method="average");
minModuleSize=30;
# 设置基因网中至少包含 30 个基因
dynamicMods=cutreeDynamic (dendro=geneTree,distM=dissTOM,deepSplit=2, pamRespectsDendro=FALSE,minClusterSize=minModuleSize);
dynamicColors=labels2colors (dynamicMods);
table (dynamicColors)
# 计算基因网的特征值
MEList=moduleEigengenes (datExpr, colors=dynamicColors)
MEs=MEList$eigengenes
MEDiss=1-cor (MEs);
METree=flashClust(as.dist (MEDiss), method="average");
# 将特征值距离小于 0.2 的基因网合并
MEDissThres=0.2
merge=mergeCloseModules (datExpr,dynamicColors,cutHeight=MEDissThres,verbose=3)
mergedColors=merge$colors;
mergedMEs=merge$newMEs;
## 重新命名合并后的基因网
moduleColors=mergedColors
colorOrder=c ("grey", standardColors (50));
moduleLabels=match(moduleColors,colorOrder)-1;
MEs=mergedMEs;
## 绘制最终基因网络构建图
plotDendroAndColors (geneTree, mergedColors,"Merged dynamic",dendroLabels=FALSE, hang=0.03,addGuide=TRUE, guideHang=0.05)

#相关系数法
#计算每个基因模块的特征值与患病状态(变量 y)的相关系数
p.values=corPvalueStudent (cor (y, datME, use="p"),nSamples=length (y))
signif(cor(y, datME, use="p"), 2);  
#输出带pval代码

#基因网显著性
#计算每个基因模块的 GS 值，并且绘制相应图形：
#计算基因的显著性
GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# 计算基因模块的显著性
ModuleSignificance=tapply (GeneSignificance, moduleColors, mean, na.rm=T)
#图形
plotModuleSignificance (GeneSignificance, moduleColors, ylim=c (0,0.5), main="Module Significance")
