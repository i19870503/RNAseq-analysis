source("http://bioconductor.org/biocLite.R")
library(DESeq2)
library(gplots)
library(gtools)
library(ggplot2)
library(edgeR)
library(heatmap3)
library(pheatmap)
library(GenomicRanges)
library(rtracklayer)
library(pvclust)
library(scatterplot3d)
library(rgl)
library(vegan)

library(org.Mm.eg.db)
library(GSEABase)
library(GOstats)
library(Category)

library(VennDiagram)


contr = list()
# for(i in 1:length(levels(sampleTable[,3]))){
#   #  for(j in 2:length(levels(sampleTable[,3]))){
#   j = i + 1
#   if(j <= length(levels(sampleTable[,3]))) { contr[i] = paste(levels(sampleTable[,3])[i], '_vs_', levels(sampleTable[,3])[j], sep='')}
#   else {  contr[i] = paste(levels(sampleTable[,3])[i], '_vs_', levels(sampleTable[,3])[1], sep='') }
# }

combn = combn(levels(sampleTable[,3]),2)
for(i in 1:dim(combn)[2]) {
  contr[i] = paste(combn[1,i], '_vs_', combn[2,i], sep='')
}

contr = unlist(contr)
deg_num = matrix(nrow = length(contr), ncol = 7)
colnames(deg_num) = c('Group','DEG(pval <= 0.05)','up-regulatory-DEG(pval <= 0.05)','down-regulatory-DEG(pval <= 0.05)','DEG(adj.p <= 0.1)','up-regulatory-DEG(adj.p <= 0.1)','down-regulatory-DEG(adj.p <= 0.1)')
for( i in 1:length(contr)){
  res = results(dds,contrast = c("condition", strsplit(contr, '_vs_')[[i]][1],  strsplit(contr, '_vs_')[[i]][2]))
  #res = results(dds,contrast = c("x", strsplit(contr, '_vs_')[[i]][1],  strsplit(contr, '_vs_')[[i]][2]))
  resSig1 <- res[ which(res$pvalue < 0.05), ]
  resSig1_up <- res[ which(res$pvalue < 0.05 & res$log2FoldChange > 0), ]
  resSig1_down <- res[ which(res$pvalue < 0.05 & res$log2FoldChange < 0), ]
  resSig2 <- res[ which(res$padj < 0.1), ]
  resSig2_up <- res[ which(res$padj < 0.1 & res$log2FoldChange > 0), ]
  resSig2_down <- res[ which(res$padj < 0.1 & res$log2FoldChange < 0), ]
  deg = rownames(resSig1)
  deg_q = rownames(resSig2)
  deg_num[i,1] = contr[i]
  deg_num[i,2] = length(rownames(resSig1))
  deg_num[i,3] = length(rownames(resSig1_up))
  deg_num[i,4] = length(rownames(resSig1_down))
  deg_num[i,5] = length(rownames(resSig2))
  deg_num[i,6] = length(rownames(resSig2_up))
  deg_num[i,7] = length(rownames(resSig2_down))
  ### Output the differentially expression table for each comparision
  write.table(res, paste(contr[i],'.xls',sep = ''), sep = "\t", col.names = NA)
  
  ### Output the vocano plot of DEG
  pdf(paste(contr[i],'.pdf',sep = ''))
  col=rep('grey60',length(res[,2]))
  col[res[,5]< 0.05]='black'
  col[res[,5]< 0.05 & res[,2]>1]='red'
  col[res[,5]< 0.05 & res[,2]<(-1)]='green'
  plot(res[,2],-log10(res[,5]),pch=16,col=col,xlab = '-Log2 Fold Change',ylab='-Log10 P-Value',main=paste('Vocano plot of ', contr[i], sep = ''))
  abline(h=-log10(0.05),lty=2)
  abline(v = 1,lty =2)
  abline(v = -1,lty =2)
  dev.off()
  
  ### Output GO KEGG enrichment result
  Mmu_gokegg(genes = rownames(resSig1), output = contr[i])
  Mmu_gokegg(genes = rownames(resSig1_up), output = paste(contr[i]), 'up-regulated_genes', sep = '-'))
  Mmu_gokegg(genes = rownames(resSig1_down), output = paste(contr[i]), 'down-regulated_genes', sep = '-'))
  
  
  
  ### Output the heatmap of DEG between comparable group, if the count od deg over 100 the heatmap will not print the gene
  pdf(paste('deg_heatmap.',contr[i],'.pdf',sep = ''), onefile=FALSE)
  if(length(deg_q) >= 100){
    pheatmap(rpkm[deg_q,],show_rownames = 0,scale = 'row', col = bluered(55))
    dev.off()
  }
  else{
    pheatmap(rpkm[deg_q,],show_rownames = 1,scale = 'row', col = bluered(55),border_color = 0)
    dev.off()
  }
}
write.table(deg_num, file = "deg_num.xls", col.names = NA, sep ='\t')
