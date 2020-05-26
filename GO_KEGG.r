# GO KEGG
BiocManager::install('org.Mmu.eg.db')
library(org.Mmu.eg.db)
library(GOstats)
library(Category)


Mmu_gokegg = function(genes = deg_gene, output = "result"){
goAnn <- get("org.Mmu.egGO")
universe <- Lkeys(goAnn)
entrezIDs <- mget(deg_gene, org.Mmu.egENSEMBL2EG, ifnotfound=NA)

params <- new("GOHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Mmu.eg.db",
              ontology="BP",
              pvalueCutoff=1,
              conditional=FALSE,
              testDirection="over")
over <- hyperGTest(params)
glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=org.Mmu.egSYMBOL, ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
})
bp <- summary(over)
bp$Symbols <- glist[as.character(bp$GOBPID)]


params <- new("GOHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Mmu.eg.db",
              ontology="MF",
              pvalueCutoff=1,
              conditional=FALSE,
              testDirection="over")
over <- hyperGTest(params)
glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=org.Mmu.egSYMBOL, ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
})
mf <- summary(over)
mf$Symbols <- glist[as.character(mf$GOMFID)]

params <- new("GOHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Mmu.eg.db",
              ontology="CC",
              pvalueCutoff=1,
              conditional=FALSE,
              testDirection="over")
over <- hyperGTest(params)
glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=org.Mmu.egSYMBOL, ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
})
cc <- summary(over)
cc$Symbols <- glist[as.character(cc$GOCCID)]

bp$Category=rep("BP",length(bp[,1]))
mf$Category=rep("MF",length(mf[,1]))
cc$Category=rep("CC",length(cc[,1]))
colnames(bp)[1]='ID'
colnames(mf)[1]='ID'
colnames(cc)[1]='ID'
all_go=rbind(bp,mf,cc)
all_go[,10]=p.adjust(all_go[,2],method = 'BH')
colnames(all_go)[10]='adj_pval'
all_go2 = data.frame(Category = all_go[,9],all_go[,1:7],adj_pval = all_go[,10],Symbols = all_go[,8])

keggAnn <- get("org.Mmu.egPATH")
universe <- Lkeys(keggAnn)
params <- new("KEGGHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Mmu.eg.db",
              categoryName="KEGG",
              pvalueCutoff=1,
              testDirection="over")
over <- hyperGTest(params)
kegg <- summary(over)
library(Category)
glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=org.Mmu.egSYMBOL, ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
})
kegg$Symbols <- glist[as.character(kegg$KEGGID)]

write.table(all_go, paste(output, 'go_enrichment.xls', sep = '_'), sep = '\t', col.names = NA)
write.table(kegg, paste(output, 'kegg_enrichment.xls', sep = '_'), sep = '\t', col.names = NA)
}
