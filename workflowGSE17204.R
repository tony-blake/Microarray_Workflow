# Workflow_for_parkisondata
# Version 3, Wednesday 29 July 2015
# Created by: Tony Blake 

# Description: This  Work Flow analyses microarray data (Human, Mouse and Rat)
# for Gene Expression in Parkinson Disease Cases and considers genes of interest that 
# related to the treatment of Parkinsons

# Genes of interest: 
# For confidentiality purposes I have chosen to use some randomly generated data as the genes of interest in lieu of the actual genes of interest

# Workflow -> Microarray data normilisation, generation of differentially expressed genes, 
# selection of expression data for g.o.i,  prodcution of heatmap, GSEA Pathway Analysis using GAGE and GO

# Install Packages
source('http://www.bioconductor.org/biocLite.R')
biocLite('knitr')
biocLite('GEOquery')
biocLite('limma')
biocLite('marray')
biocLite('simpleaffy')
bioclite('affy')
biocLite('gplots')
#biocLite('biomaRt')
biocLite('hgu133plus2.db')
biocLite('hgu133a2.db')
biocLite('rgu34a.db')
biocLite('mgu74a.db')
biocLite('gage')
biocLite('gageData')
biocLite('pathview')
biocLite('stats')
#biocLite('multtest')
library(bioDist)

# Load librarys
library(knitr)
library(GEOquery)
library(limma)
library(marray)
library(simpleaffy)
library(affy)
library(gplots)
#library(biomaRt)
library(hgu133plus2.db)
library(hgu133a2.db)
library(rgu34a.db)
library(mgu74a.db)
library(gage)
library(gageData)
library(pathview)
library(stats)
library(plyr)
library(GOstats)
library(multtest)
library(xlsx)
library(reshape2)
############################################################################################################################################################################
################################################################  Begin Paramter Block #####################################################################################
runName <- "GSE17204" # string that is applied to all output files.
runName2 <- "mGSE17204"
inputDirectory <- "/Users/tonyblake/Desktop/Bioinformatics/parkinson_project/FinalCorrectedGSE/GSE17204/data_GSE17204" # Path to directory that contains the .CEL files and phenodata text file.
outputDirectory <- "/Users/tonyblake/Desktop/Bioinformatics/parkinson_project/FinalCorrectedGSE/GSE17204/Output" # Path to directory that outputs heatmaps, tables, etc 

baseline <- c("s_control_m", "s_control_h") # string to convey baseline condition from target column in phnodata text file
minFC <- 1.5 # This is the minimum log2 fold change for a gene to be differentially expressed.
ttestPVal <- 0.01 # This is the threshold p value for significance. 
hgCutoff <- 0.01 # This is the GOStats p value threshold.

################################################################ End Parameter Block #######################################################################################

# Download raw data file to a named data directory and uncompress files

getGEOSuppFiles("GSE17204")
untar("GSE17204/GSE17204_RAW.tar", exdir="data_GSE17204") # uncompress folder
cels <- list.files("data_GSE17204/", pattern = "[gz]") # list all files ending in ".gz" as character strings
sapply(paste("data_GSE17204", cels, sep="/"), gunzip) # unzip .CEL files in data directory

################################################# Begin  "phenodata.txt" file creation from commandline #############################################################################

# $ cd /Users/tonyblake/Desktop/Bioinformatics/parkinson_project/GSE17204/data_GSE17204
# $ ls *.CEL > phenodata.txt
# ....open in spreadsheet and create extra columns of data, copy all text and use "paste and match style" command to paste over original text in "phenodata.txt" file
#  so that it looks like 
#
# Name  File_Name  Target
# GSM430339.CEL  GSM430339.CEL  antisense_DJI_B
# GSM430340.CEL  GSM430340.CEL	antisense_DJI_B
# GSM430341.CEL	GSM430341.CEL	antisense_DJI_G
# GSM430342.CEL	GSM430342.CEL	antisense_DJI_G
# GSM430343.CEL	GSM430343.CEL	s_control_m
# GSM430344.CEL	GSM430344.CEL	s_control_m
# GSM430345.CEL	GSM430345.CEL	s_control_h
# GSM430346.CEL	GSM430346.CEL	s_control_h

####################################################### End "phenodata.txt" file creation ################################################################################## 

############################## Begin Normalisation Block ########################################################################################################

celfiles <- read.affy(covdesc="phenodata.txt", path=inputDirectory) # creates affybatch object
expression <- expresso(celfiles, normalize.method='quantiles', bgcorrect.method='rma', pmcorrect.method='pmonly',summary.method='medianpolish')
exprs <- exprs(expression) # extracts logfold change (expression) data
conditions <- pData(expression)$Target
allunique_Conditions <- unique(conditions)
otherconditions <- allunique_Conditions[allunique_Conditions !=baseline]
reversetwo <- rev(otherconditions)
allotherconditions <- c(otherconditions, reversetwo)
nocontrols <- paste(otherconditions[1],"-", otherconditions[2])

############################## End Normalisation Block ###########################################################################################################

################ Selection of Differentially Expressed Genes ################


## Limma DEG calculations are done here:
gc()
f <- factor(conditions, levels=allunique_Conditions)
design <- model.matrix(~0+f)
colnames(design) <- allunique_Conditions
fit <- lmFit(expression, design)
contrast <- paste(allotherconditions,"-", c(baseline,baseline))
contrast.matrix <- makeContrasts(contrast[1],contrast[2],contrast[3],contrast[4], nocontrols,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

##############################################################################################################
### link gene symbols to affymatrix  ######
exprs <- exprs(expression)
affyid=rownames(exprs)
symbol2probe=hgu133a2SYMBOL[affyid]
annots=toTable(symbol2probe)
View(annots)

##### create some random data for gene symbols (note that in proper use you would not "create random data" but use a real dataset
##### such as a list of genes of interest for your particular project)
genes <- annots$symbol
random_genes <- sample(genes, size=36, replace=TRUE)
random_genes <- as.vector(random_genes)
goi_symbols <- data.frame(random_genes)
genesofinterest <- goi_symbols$random_genes
genesofinterest <- as.character(genesofinterest)
mygenelist <- mget(genesofinterest, revmap(hgu133a2SYMBOL), ifnotfound=NA) #use for human datasets only! 
mygenedf <- do.call("rbind", lapply(mygenelist,data.frame)) # creates dataframe from gene and symbol list
which(is.na(mygenelist)) # Gives names of genes of interest  not included in affymatrix
u <- mygenedf$X..i..
mygenedf_nona <- mygenedf[!is.na(u),] # removes NA's from dataframe of genes
mygenedf_nona <- as.character(mygenedf_nona)

# link probeid's to genesymbols
goi_linkedprobes <- select(hgu133a2.db, keys= mygenedf_nona, columns = "SYMBOL", keytype = "PROBEID") #for human dataset only!
affys <- goi_linkedprobes$PROBEID

#F-statistics
fit2table <- fit2
fit2table <- fit2table[affys,]
fstat_table <- merge(fit2table, goi_linkedprobes, by.x=0, by.y="PROBEID")
View(fstat_table)



# loop through limmma coef's to create 3 limma tables for the 3 comparisons
tops <- list()
symbol_linkedtogenes <- list()


for (i in 1:5){
  tops[[i]] <- topTable(fit2, number= Inf,coef=i)[affys,] #adjust="none", p.value=ttestPVal, lfc=log2(minFC))
  symbol_linkedtogenes[[i]] <- merge(tops[[i]], goi_linkedprobes, by.x=0, by.y="PROBEID")
  col_1st_goi <- grep("SYMBOL", names(symbol_linkedtogenes[[i]]))
  symbol_linkedtogenes[[i]] <- symbol_linkedtogenes[[i]][,c(col_1st_goi, (1:ncol(symbol_linkedtogenes[[i]]))[-col_1st_goi])]
  colnames(symbol_linkedtogenes[[i]])[2] <- "ProbesetIds"
  entrezs <- base::unlist(AnnotationDbi::as.list(hgu133a2ENTREZID[affys]))
  ddd <- ldply(entrezs)
  ddd <- ddd[order(ddd$.id),]
  symbol_linkedtogenes[[i]]$ENTREZID<-ddd$V1 # Adds ENTREZID column
  
}

# check that forloop is giving right result

check1 <- topTable(fit2, number= Inf,coef=4)[affys,] #adjust="none", p.value=ttestPVal, lfc=log2(minFC))
linkedtogenes1 <- merge(check1, goi_linkedprobes, by.x=0, by.y="PROBEID")
col_1st_goi1 <- grep("SYMBOL", names(linkedtogenes1))
linkedtogenes1 <- linkedtogenes1[,c(col_1st_goi1, (1:ncol(linkedtogenes1))[-col_1st_goi1])]
colnames(linkedtogenes1)[2] <- "ProbesetIds"
entrezs <- unlist(AnnotationDbi::as.list(hgu133a2ENTREZID[affys]))
ddd <- ldply(entrezs)
ddd <- ddd[order(ddd$.id),]
linkedtogenes1$ENTREZID<-ddd$V1 # Adds ENTREZID column

all.equal(symbol_linkedtogenes[[4]], linkedtogenes1)


# Final versions of Limma tables
tops1 <- symbol_linkedtogenes[[1]]
tops2 <- symbol_linkedtogenes[[2]]
tops3 <- symbol_linkedtogenes[[3]]
tops4 <- symbol_linkedtogenes[[4]]
tops5 <- symbol_linkedtogenes[[5]]




write.xlsx(tops1, file="limma_table_GSE17204_DJI_BvCtrlM.xlsx")
write.xlsx(tops2, file="limma_table_GSE17204_DJI_GvCtrlH.xlsx")
write.xlsx(tops3, file="limma_table_GSE17204_DJI_GvCtrlM.xlsx")
write.xlsx(tops4, file="limma_table_GSE17204_DJI_BvCtrlH.xlsx")
write.xlsx(tops5, file="limma_table_GSE17204_DJI_BvsDJI_G.xlsx")


# create expression tables for heatmap

exprs=exprs(expression)[affys,]
affyid=rownames(exprs)
symbol2probe=hgu133a2SYMBOL[affyid]
annots=toTable(symbol2probe)
str(annots)
#check to see if dataframes are same
df1 <- sort(annots$probe_id, decreasing=TRUE) 
df2 <- sort(goi_linkedprobes$PROBEID, decreasing=TRUE)
all.equal(df1,df2)
#back to main workflow
exprs=exprs[annots$probe_id,]

#if multiple probe sets map to a gene, select the one with maximal IQR
iqrs=apply(exprs, 1, IQR)
sel.rn=tapply(1:nrow(annots), annots$symbol, function(x){
  x[which.max(iqrs[x])]
})
exprs.symb=exprs[sel.rn,] #expression data with single gene names
rownames(exprs.symb)=names(sel.rn)

mexprs.symb <- exprs
rownames(mexprs.symb)=annots$symbol
mexprs.symb <- mexprs.symb[order(rownames(mexprs.symb)),] #expression data where there are more than one of the smae gene

################################################################# End Selection of Differentially Expressed Genes ############

################################################################# Begin Heatmap Block ##########################################



# heatmap for single genes
pdfOutputFile <- paste(outputDirectory, "/", runName, "_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="s_control_h") "#FF0000" else {if (Target=="s_control_m") "#0000FF" else {if (Target=="antisense_DJI_G") "#2dbdfc" else "#660060"}} }
patientcolors <- base::unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(exprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Rowv=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=1, cexCol=.55)
dev.off()

# heatmap for multiple genes
pdfOutputFile <- paste(outputDirectory, "/", runName2, "_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="s_control_h") "#FF0000" else {if (Target=="s_control_m") "#0000FF" else {if (Target=="antisense_DJI_G") "#2dbdfc" else "#660060"}} }
patientcolors <- base::unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(mexprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Rowv=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=.7, cexCol=.55)
dev.off()

# heatmap for single genes with gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName , "genecluster_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="s_control_h") "#FF0000" else {if (Target=="s_control_m") "#0000FF" else {if (Target=="antisense_DJI_G") "#2dbdfc" else "#660060"}} }
patientcolors <- base::unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(exprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=1, cexCol=.55)
dev.off()

# heatmap for multiple genes with gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName2 , "_mgenecluster_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="s_control_h") "#FF0000" else {if (Target=="s_control_m") "#0000FF" else {if (Target=="antisense_DJI_G") "#2dbdfc" else "#660060"}} }
patientcolors <- unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(mexprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Colv=FALSE,density.info="none", trace="none", cexRow=.7, cexCol=.55)
dev.off()

# heatmap for multiple genes with gene clustering
pdfOutputFile <- paste(outputDirectory, "/", runName2 , "_mgenebothcluster_heatmap.pdf", sep="") # Sets file name.
pdf(pdfOutputFile, width=7,height= 7.5)
color.map <- function(Target) { if (Target=="s_control_h") "#FF0000" else {if (Target=="s_control_m") "#0000FF" else {if (Target=="antisense_DJI_G") "#2dbdfc" else "#660060"}} }
patientcolors <- unlist(lapply(pData(celfiles)$Target, color.map))
heatmap.2(mexprs.symb, col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, Colv=TRUE,density.info="none", trace="none", cexRow=.7, cexCol=.55)
dev.off()

################################### End Heatmap Block ##################################################################

################################### Begin Single gene Limma tables ###########################################################

sannots <- annots[sel.rn,]

# Final versions of Limma tables
stops1 <- tops1[rownames(sannots),] 
stops2 <- tops2[rownames(sannots),]
stops3 <- tops3[rownames(sannots),] 
stops4 <- tops4[rownames(sannots),]
stops5 <- tops5[rownames(sannots),]

write.xlsx(stops1, file="slimma_table_GSE17204_DJI_BvCtrlM.xlsx")
write.xlsx(stops2, file="slimma_table_GSE17204_DJI_GvCtrlH.xlsx")
write.xlsx(stops3, file="slimma_table_GSE17204_DJI_GvCtrlM.xlsx")
write.xlsx(stops4, file="slimma_table_GSE17204_DJI_BvCtrlH.xlsx")
write.xlsx(stops5, file="slimma_table_GSE17204_DJI_BvsDJI_G.xlsx")

################################## End Single gene Limma tables #########################################################

################################## Begin GOstats Block ##################################################################


x <- org.Hs.egACCNUM
mapped_genes <- mappedkeys(x)
xx <- AnnotationDbi::as.list(x[mapped_genes])
geneUniverse <- (unique(names(xx)))


# Gene Ontology Categories for condition  that were shown to be relatively Higher (more expressed) 
# for (1) "antisense_DJI_B than s_control_m" 
# for (2) "antisense_DJI_G than s_control_h" 
# for (3) "antisense_DJI_B than s_control_m" 
# for (4) "antisense_DJI_G than s_control_h"
# for (5) "antisense_DJI_B than antisense_DJI_G"

sig_goi <- list()
upreg_goi <- list()
up_gostat <- list()
upgenesCHR <- list()
upgenesLinkedtoEntrezIds <- list()
GOstats_upgenes <- list()
paramsBPH <- list()
paramsCCH <- list()
paramsMFH <- list()
BP_GOdataH <- list()
BP_SummH <- list()
CC_GOdataH <- list()
CC_SummH <- list()
MF_GOdataH <- list()
MF_SummH <- list()

# No significant up regulated genes for 1:3 

for(i in 4:5){
  
  sig_goi[[i]]<-subset(symbol_linkedtogenes[[i]],symbol_linkedtogenes[[i]][,7]<0.05)
  upreg_goi[[i]] <- subset(sig_goi[[i]], sig_goi[[i]][,3]>0)
  up_gostat[[i]] <- upreg_goi[[i]][1]
  upgenesCHR[[i]] <- up_gostat[[i]]$SYMBOL
  upgenesLinkedtoEntrezIds[[i]] <- select(hgu133a2.db, keys= upgenesCHR[[i]], "ENTREZID", "SYMBOL")
  GOstats_upgenes[[i]] <- upgenesLinkedtoEntrezIds[[i]][,2]
  
  paramsBPH[[i]] <- new("GOHyperGParams", geneIds=GOstats_upgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsCCH[[i]] <- new("GOHyperGParams", geneIds=GOstats_upgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="CC", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsMFH[[i]] <- new("GOHyperGParams", geneIds=GOstats_upgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="MF", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  BP_GOdataH[[i]] <- hyperGTest(paramsBPH[[i]]) # Biological Process.
  BP_SummH[[i]] <- summary(BP_GOdataH[[i]])
  
  CC_GOdataH[[i]] <- hyperGTest(paramsCCH[[i]]) # Cellular Component.
  CC_SummH[[i]] <- summary(CC_GOdataH[[i]])
  
  MF_GOdataH[[i]] <- hyperGTest(paramsMFH[[i]]) # Molecular Function.
  MF_SummH[[i]] <- summary(MF_GOdataH[[i]])
  
}

# Gene Ontology Categories for condition  that were shown to be relatively Lower (less expressed) 
# for (1) "antisense_DJI_B than s_control_m" 
# for (2) "antisense_DJI_G than s_control_h" 
# for (3) "antisense_DJI_B than s_control_m" 
# for (4) "antisense_DJI_G than s_control_h"
# for (5) "antisense_DJI_B than antisense_DJI_G"


sig_goi <- list()
downreg_goi <- list()
downgenesLinkedtoEntrezIds <- list()
GOstats_downgenes <- list()
paramsBPL <- list()
paramsCCL <- list()
paramsMFL <- list()
BP_GOdataL <- list()
BP_SummL <- list()
CC_GOdataL <- list()
CC_SummL <- list()
MF_GOdataL <- list()
MF_SummL <- list()

for(i in 1:2){
  
  sig_goi[[i]]<-subset(symbol_linkedtogenes[[i]],symbol_linkedtogenes[[i]][,7]<0.05)
  downreg_goi[[i]]<-subset(sig_goi[[i]], sig_goi[[i]][,3]<0)
  down_gostat <- downreg_goi[[i]][1]
  downgenesCHR <- down_gostat$SYMBOL
  downgenesLinkedtoEntrezIds[[i]] <- select(hgu133a2.db, keys= downgenesCHR, "ENTREZID", "SYMBOL")
  GOstats_downgenes[[i]] <- downgenesLinkedtoEntrezIds[[i]][,2]
  
  paramsBPL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsCCL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="CC", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsMFL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="MF", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  BP_GOdataL[[i]] <- hyperGTest(paramsBPL[[i]]) # Biological Process.
  BP_SummL[[i]] <- summary(BP_GOdataL[[i]])
  
  CC_GOdataL[[i]] <- hyperGTest(paramsCCL[[i]]) # Cellular Component.
  CC_SummL[[i]] <- summary(CC_GOdataL[[i]])
  
  MF_GOdataL[[i]] <- hyperGTest(paramsMFL[[i]]) # Molecular Function.
  MF_SummL[[i]] <- summary(MF_GOdataL[[i]])
  
}

for(i in 4:5){
  
  sig_goi[[i]]<-subset(symbol_linkedtogenes[[i]],symbol_linkedtogenes[[i]][,7]<0.05)
  downreg_goi[[i]]<-subset(sig_goi[[i]], sig_goi[[i]][,3]<0)
  down_gostat <- downreg_goi[[i]][1]
  downgenesCHR <- down_gostat$SYMBOL
  downgenesLinkedtoEntrezIds[[i]] <- select(hgu133a2.db, keys= downgenesCHR, "ENTREZID", "SYMBOL")
  GOstats_downgenes[[i]] <- downgenesLinkedtoEntrezIds[[i]][,2]
  
  paramsBPL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="BP", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsCCL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="CC", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  paramsMFL[[i]] <- new("GOHyperGParams", geneIds=GOstats_downgenes[[i]], universeGeneIds=geneUniverse,
                        annotation="org.Hs.eg.db", ontology="MF", pvalueCutoff=hgCutoff, conditional=TRUE,
                        testDirection="over")
  
  BP_GOdataL[[i]] <- hyperGTest(paramsBPL[[i]]) # Biological Process.
  BP_SummL[[i]] <- summary(BP_GOdataL[[i]])
  
  CC_GOdataL[[i]] <- hyperGTest(paramsCCL[[i]]) # Cellular Component.
  CC_SummL[[i]] <- summary(CC_GOdataL[[i]])
  
  MF_GOdataL[[i]] <- hyperGTest(paramsMFL[[i]]) # Molecular Function.
  MF_SummL[[i]] <- summary(MF_GOdataL[[i]])
  
}


## Write out the results
source_https('https://raw.github.com/bobthecat/codebox/master/write.GOhyper.r')

bph1GO <- write.GOhyper(BP_GOdataH[[1]], filename="hBP_GSE17204_DJI_BvCtrlM.xls")
bph2GO <- write.GOhyper(BP_GOdataH[[2]], filename="hBP_GSE17204_DJI_GvCtrlH.xls")
bph3GO <- write.GOhyper(BP_GOdataH[[3]], filename="hBP_GSE17204_DJI_GvsCtrlM.xls")
bph4GO <- write.GOhyper(BP_GOdataH[[4]], filename="hBP_GSE17204_DJI_BvsCtrlH.xls")
bph5GO <- write.GOhyper(BP_GOdataH[[5]], filename="hBP_GSE17204_DJI_BvsDJI_G.xls") 

cch1GO <- write.GOhyper(CC_GOdataH[[1]], filename="hCC_GSE17204_DJI_BvCtrlM.xls")
cch2GO <- write.GOhyper(CC_GOdataH[[2]], filename="hCC_GSE17204_DJI_GvCtrlH.xls")
cch3GO <- write.GOhyper(CC_GOdataH[[3]], filename="hCC_GSE17204_DJI_GvsCtrlM.xls")
cch4GO <- write.GOhyper(CC_GOdataH[[4]], filename="hCC_GSE17204_DJI_BvsCtrlH.xls")
cch5GO <- write.GOhyper(CC_GOdataH[[5]], filename="hCC_GSE17204_DJI_BvsDJI_G.xls") 

mfh1GO <- write.GOhyper(MF_GOdataH[[1]], filename="hMF_GSE17204_DJI_BvCtrlM.xls")
mfh2GO <- write.GOhyper(MF_GOdataH[[2]], filename="hMF_GSE17204_DJI_GvCtrlH.xls")
mfh3GO <- write.GOhyper(MF_GOdataH[[3]], filename="hMF_GSE17204_DJI_GvsCtrlM.xls")
mfh4GO <- write.GOhyper(MF_GOdataH[[4]], filename="hMF_GSE17204_DJI_BvsCtrlH.xls")
mfh5GO <- write.GOhyper(MF_GOdataH[[5]], filename="hMF_GSE17204_DJI_BvsDJI_G.xls") 

bpl1GO <- write.GOhyper(BP_GOdataL[[1]], filename="lBP_GSE17204_DJI_BvCtrlM.xls")
bpl2GO <- write.GOhyper(BP_GOdataL[[2]], filename="lBP_GSE17204_DJI_GvCtrlH.xls")
bpl3GO <- write.GOhyper(BP_GOdataL[[3]], filename="lBP_GSE17204_DJI_GvsCtrlM.xls")
bpl4GO <- write.GOhyper(BP_GOdataL[[4]], filename="lBP_GSE17204_DJI_BvsCtrlH.xls")
bpl5GO <- write.GOhyper(BP_GOdataL[[5]], filename="lBP_GSE17204_DJI_BvsDJI_G.xls") 

ccl1GO <- write.GOhyper(CC_GOdataL[[1]], filename="lCC_GSE17204_DJI_BvCtrlM.xls")
ccl2GO <- write.GOhyper(CC_GOdataL[[2]], filename="lCC_GSE17204_DJI_GvCtrlH.xls")
ccl3GO <- write.GOhyper(CC_GOdataL[[3]], filename="lCC_GSE17204_DJI_GvsCtrlM.xls")
ccl4GO <- write.GOhyper(CC_GOdataL[[4]], filename="lCC_GSE17204_DJI_BvsCtrlH.xls")
ccl5GO <- write.GOhyper(CC_GOdataL[[5]], filename="lCC_GSE17204_DJI_BvsDJI_G.xls") 

mfl1GO <- write.GOhyper(MF_GOdataL[[1]], filename="lMF_GSE17204_DJI_BvCtrlM.xls")
mfl2GO <- write.GOhyper(MF_GOdataL[[2]], filename="lMF_GSE17204_DJI_GvCtrlH.xls")
mfl3GO <- write.GOhyper(MF_GOdataL[[3]], filename="lMF_GSE17204_DJI_GvsCtrlM.xls")
mfl4GO <- write.GOhyper(MF_GOdataL[[4]], filename="lMF_GSE17204_DJI_BvsCtrlH.xls")
mfl5GO <- write.GOhyper(MF_GOdataL[[5]], filename="lMF_GSE17204_DJI_BvsDJI_G.xls") 

###################################### Colour Map #############################################

# color scale
c1 <- rainbow(32,v=seq(0.5,1,length=32),s=seq(1,0.3,length=32),start=4/6,end=4.0001/6);
pie(rep(1,32),col=c1);   # show off these colors
c2 <- rainbow(32,v=seq(0.5,1,length=32),s=seq(1,0.3,length=32),start=1/6,end=1.0001/6);
pie(rep(1,32),col=c2);   # show off these colors
c3 <- c(c1,rev(c2));  # rev reverses the list
#rm(c1,c2)


######################################## Correlation Analysis ###################################
 
# Here the degree of correlation between different gene ontologies is analysed based on the 
# the differentailly expressed genes that they have in common.

# The following pieces of code take the objects created from the GOstats analysis as input. It then creates a matrix of 
# correlation values for gene ontologies and uses it to create a correlogram showing the percenatge of genes shared 
# between gene ontologies in the form a clustered heatmap. It also creates a table of gene ontologies showing the over
# represented genes in each ontology. Rows are ordered according to the clustering heatmap of correlation values. 
 
 
 
#upregulated pathways

#bph4

bph4square <- makebeta4_pathwaygenetable(bph4GO) 
hr <- hclust(cor.dist(bph4square[[1]]))
hc <- hclust(cor.dist(bph4square[[1]]))
pdf(file='hBP_correlmap_GSE17204_antisense_DJI_Bvs_control_h.pdf', height=13, width=13)
heatmap.2(bph4square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.1, cexCol=0.1)
dev.off()
bph1clustered <- bph4square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(bph4clustered)
clustered <- bph4square[[2]][match(clst, bph4square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
#write.xlsx(my.df, file="hBP_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")
write.xlsx(bph4square[[2]], file="hBP_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")


#cch4
cch4square <- makebeta4_pathwaygenetable(cch4GO)
hr <- hclust(cor.dist(cch4square[[1]]))
hc <- hclust(cor.dist(cch4square[[1]]))
pdf(file='hCC_correlmap_GSE17204_antisense_DJI_Bvs_control_h.pdf', height=13, width=13)
heatmap.2(cch4square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
cch4clustered <- cch4square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(cch4clustered)
clustered <- cch4square[[2]][match(clst, cch4square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
#write.xlsx(my.df, file="hCC_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")
write.xlsx(cch4square[[2]], file="hCC_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")

#mfh4
mfh4square <- makebeta4_pathwaygenetable(mfh4GO)
hr <- hclust(cor.dist(mfh4square[[1]]))
hc <- hclust(cor.dist(mfh4square[[1]]))
pdf(file='hMF_correlmap_GSE17204_antisense_DJI_Bvs_control_h.pdf', height=13, width=13)
heatmap.2(mfh4square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
mfh4clustered <- mfh4square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(mfh4clustered)
clustered <- mfh4square[[2]][match(clst, mfh4square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
#write.xlsx(my.df, file="hMF_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")
write.xlsx(bph4square[[2]], file="hMF_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")

#down regulated

#bpl1

bpl1square <- makebeta1down_pathwaygenetable(bpl1GO)
hr <- hclust(cor.dist(bpl1square[[1]]))
hc <- hclust(cor.dist(bpl1square[[1]]))
pdf(file='lBP_correlmap_GSE17204_antisense_DJI_Bvs_control_m.pdf', height=13, width=13)
heatmap.2(bpl1square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.1, cexCol=0.1)
dev.off()
bpl1clustered <- bpl1square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(bpl1clustered)
clustered <- bpl1square[[2]][match(clst, bpl1square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
#write.xlsx(my.df, file="hMF_genesinpathways_GSE17204_antisense_DJI_Bvs_control_m.xlsx")
#write.xlsx(bpl1square[[2]], file="lBP_genesinpathways_GSE17204_antisense_DJI_Bvs_control_m.xlsx")


bpl2square <- makebeta2down_pathwaygenetable(bpl2GO)
bpl2square[[1]][is.nan(bpl2square[[1]])] = 0
hr <- hclust(cor.dist(bpl2square[[1]]))
hc <- hclust(cor.dist(bpl2square[[1]]))
pdf(file='lBP_correlmap_GSE17204_antisense_DJI_Gvs_control_h.pdf', height=13, width=13)
heatmap.2(bpl2square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.3, cexCol=0.3)
dev.off()
bpl2clustered <- bpl2square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(bpl2clustered)
clustered <- bpl2square[[2]][match(clst, bpl2square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lBP_genesinpathways_GSE17204_antisense_DJI_Gvs_control_h.xlsx")
#write.xlsx(bpl2square[[2]], file="lBP_genesinpathways_GSE17204_antisense_DJI_Gvs_control_h.xlsx")


bpl4square <- makebeta4down_pathwaygenetable(bpl4GO)
#bpl2square[[1]][is.nan(bpl2square[[1]])] = 0
hr <- hclust(cor.dist(bpl4square[[1]]))
hc <- hclust(cor.dist(bpl4square[[1]]))
pdf(file='lBP_correlmap_GSE17204_antisense_DJI_Bvs_control_h.pdf', height=13, width=13)
heatmap.2(bpl4square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.3, cexCol=0.3)
dev.off()
bpl4clustered <- bpl4square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(bpl4clustered)
clustered <- bpl4square[[2]][match(clst, bpl4square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lBP_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")

#ccl2
ccl2square <- makebeta2down_pathwaygenetable(ccl2GO)
hr <- hclust(cor.dist(ccl2square[[1]]))
hc <- hclust(cor.dist(ccl2square[[1]]))
pdf(file='lCC_correlmap_GSE17204_antisense_DJI_Gvs_control_h.pdf', height=13, width=13)
heatmap.2(ccl2square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
ccl2clustered <- ccl2square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(ccl2clustered)
clustered <- ccl2square[[2]][match(clst, ccl2square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lCC_genesinpathways_GSE17204_antisense_DJI_Gvs_control_h.xlsx")

#ccl4
ccl4square <- makebeta4down_pathwaygenetable(ccl4GO)
hr <- hclust(cor.dist(ccl4square[[1]]))
hc <- hclust(cor.dist(ccl4square[[1]]))
pdf(file='lCC_correlmap_GSE17204_antisense_DJI_Bvs_control_h.pdf', height=13, width=13)
heatmap.2(ccl4square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
ccl4clustered <- ccl4square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(ccl4clustered)
clustered <- ccl4square[[2]][match(clst, ccl4square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lCC_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")




#mfl2
mfl2square <- makebeta2down_pathwaygenetable(mfl2GO)
hr <- hclust(cor.dist(mfl2square[[1]]))
hc <- hclust(cor.dist(mfl2square[[1]]))
pdf(file='lMF_correlmap_GSE17204_antisense_DJI_Gvs_control_h.pdf', height=13, width=13)
heatmap.2(mfl2square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
mfl2clustered <- mfl2square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(mfl2clustered)
clustered <- mfl2square[[2]][match(clst, mfl2square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lMF_genesinpathways_GSE17204_antisense_DJI_Gvs_control_h.xlsx")

#mfl4
mfl4square <- makebeta4down_pathwaygenetable(mfl4GO)
hr <- hclust(cor.dist(mfl4square[[1]]))
hc <- hclust(cor.dist(mfl4square[[1]]))
pdf(file='lMF_correlmap_GSE17204_antisense_DJI_Bvs_control_h.pdf', height=13, width=13)
heatmap.2(mfl4square[[1]], col=c3, margins = c(20, 20),  Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), trace='none', density.info='none', cexRow=0.7, cexCol=0.7)
dev.off()
mfl4clustered <- mfl4square[[1]][rev(hr$labels[hr$order]), hc$labels[hc$order]]
clst <- rownames(mfl4clustered)
clustered <- mfl4square[[2]][match(clst, mfl4square[[2]]$pathways),]
my.df <- data.frame(lapply(clustered, as.character), stringsAsFactors=FALSE)
write.xlsx(my.df, file="lMF_genesinpathways_GSE17204_antisense_DJI_Bvs_control_h.xlsx")

############################################# End #############################################################################
