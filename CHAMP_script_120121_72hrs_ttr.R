if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("ChAMP")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DMRcate")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylumi")

library("ChAMP")
#Delete environment 
rm(beta_heatmap, df, df_subt, heat, newdata, newdata2, newdata3, newdata4, newdata5)

#Local lab computer
setwd("C:/Users/bossarr/Desktop/r_aray")

#Local Mac
setwd("~/Analises R/r_array_huvec_120121/IL6_72hrs_vs_PBS")

myLoad <- champ.load(directory = "~/Analises R/r_array_huvec_120121/IL6_72hrs_vs_PBS",arraytype = "EPIC")
myLoad_2 <- champ.load("~/Analises R/r_array_huvec_120121/IL6_72hrs_vs_PBS",method="minfi",arraytype="EPIC") # loaded as minfi way

CpG.GUI(arraytype="EPIC")

#Save specific data
#save(Anno, EPIC.manifest.hg19, multi.hit, myLoad, probe.features, bloodCtl, MatchGeneName, myCNA,
#     myDMP, myNorm, PathwayList, probeInfoALL.lv, file = "champ_huvec_120121_72hrs.RData")

#Save all 

save.image(file = "champ_huvec_120121_72hrs.RData")
load("champ_huvec_120121_72hrs.RData")

#Open the files - save all the idat on the pathway below:

testDir=system.file("extdata",package="ChAMPdata")
myLoad <- champ.load(testDir,arraytype="EPIC")

#Save the .idat and .csv files on the folder below, this will make automatic upload for all the samples
#/Users/ramonbossardi/Library/R/4.0/library/ChAMPdata/extdata/lung_test_set.csv

#QC.GUI(beta=myLoad$beta)
#NOrmalization method of BMIQ
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC")
#NOrmalization method of PBC
myNorm_2 <- champ.norm(method="PBC",arraytype = "EPIC")
#NOrmalization method of SWAN
myNorm_3 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="SWAN",arraytype="EPIC")
#NOrmalization method of FunctionalNormalization
myNorm_4 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="FunctionalNormalization",arraytype="EPIC")


#QC.GUI for BMIQ
QC.GUI(beta=myNorm,arraytype = "EPIC")

#QC.GUI for BMIQ
QC.GUI(beta=myNorm_3,arraytype = "EPIC")

champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)


#Slide is pd file is numeric, which is not correct, so we firstly manually change it into charater
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
myLoad$pd$Array <- as.factor(myLoad$pd$Array)
#champ.SVD() will only check the dimensions between data and pd, instead if checking if Sample_Names are correctly matched (because some user may have no Sample_Names in their pd file),thus please make sure your pd file is in accord with your data sets (beta) and (rgSet).
champ.SVD(beta=myNorm)

table(myLoad$pd$Slide)

myDMP <- champ.DMP(arraytype = "EPIC")
##############################################################################################

#Experiment 2 and 3_just duplicate 1

library("ChAMP")

#Local Mac
setwd("~/Analises R/r_array_huvec_120121/IL6_72hrs_exp2_3")

myLoad <- champ.load(directory = "~/Analises R/r_array_huvec_120121/IL6_72hrs_exp2_3",arraytype = "EPIC")
myLoad_2 <- champ.load("~/Analises R/r_array_huvec_120121/IL6_72hrs_exp2_3",method="minfi",arraytype="EPIC") # loaded as minfi way

CpG.GUI(arraytype="EPIC")

#Save specific data
#save(Anno, EPIC.manifest.hg19, multi.hit, myLoad, probe.features, bloodCtl, MatchGeneName, myCNA,
#     myDMP, myNorm, PathwayList, probeInfoALL.lv, file = "champ_huvec_120121_72hrs.RData")

#Save all 

save.image(file = "champ_huvec_120221_72hrs_exp2_3_1st duplicate.RData")
load("champ_huvec_120121_72hrs.RData")

#QC.GUI(beta=myLoad$beta)
#NOrmalization method of BMIQ
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC")
#NOrmalization method of PBC
myNorm_2 <- champ.norm(method="PBC",arraytype = "EPIC")
#NOrmalization method of SWAN
myNorm_3 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="SWAN",arraytype="EPIC")
#NOrmalization method of FunctionalNormalization
myNorm_4 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="FunctionalNormalization",arraytype="EPIC")


#QC.GUI for BMIQ
QC.GUI(beta=myNorm,arraytype = "EPIC")

#QC.GUI for SWAN
QC.GUI(beta=myNorm_3,arraytype = "EPIC")

champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)


#Slide is pd file is numeric, which is not correct, so we firstly manually change it into charater
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
myLoad$pd$Array <- as.factor(myLoad$pd$Array)
#champ.SVD() will only check the dimensions between data and pd, instead if checking if Sample_Names are correctly matched (because some user may have no Sample_Names in their pd file),thus please make sure your pd file is in accord with your data sets (beta) and (rgSet).
champ.SVD(beta=myNorm)

table(myLoad$pd$Slide)

myDMP <- champ.DMP( adjPVal = 0.1, arraytype = "EPIC")

##############################################################################################

#Experiment 2 and 3 all samples

library("ChAMP")

#Local Mac
setwd("~/Analises R/r_array_huvec_120121/IL6_72hrs_exp2_3_allsamples")

myLoad <- champ.load(directory = "~/Analises R/r_array_huvec_120121/IL6_72hrs_exp2_3_allsamples",arraytype = "EPIC")
myLoad_2 <- champ.load("~/Analises R/r_array_huvec_120121/IL6_72hrs_exp2_3_allsamples",method="minfi",arraytype="EPIC") # loaded as minfi way

CpG.GUI(arraytype="EPIC")

#Save specific data
#save(Anno, EPIC.manifest.hg19, multi.hit, myLoad, probe.features, bloodCtl, MatchGeneName, myCNA,
#     myDMP, myNorm, PathwayList, probeInfoALL.lv, file = "champ_huvec_120121_72hrs.RData")

#Save all 

save.image(file = "champ_huvec_120121_72hrs_exp2_3_allsamples.RData")
load("champ_huvec_120121_72hrs_exp2_3_allsamples.RData")

#QC.GUI(beta=myLoad$beta)
#NOrmalization method of BMIQ
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC")
#NOrmalization method of PBC
myNorm_2 <- champ.norm(method="PBC",arraytype = "EPIC")
#NOrmalization method of SWAN
myNorm_3 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="SWAN",arraytype="EPIC")
#NOrmalization method of FunctionalNormalization
myNorm_4 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="FunctionalNormalization",arraytype="EPIC")


#QC.GUI for BMIQ
QC.GUI(beta=myNorm,arraytype = "EPIC")

#QC.GUI for SWAN
QC.GUI(beta=myNorm_3,arraytype = "EPIC")

champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)
champ.SVD(beta=myNorm, resultsDir="./CHAMP_SVDimages_withoutnormalization/")

############ Batch Effect Correction
#This function implements the ComBat normalization method
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))

#Slide is pd file is numeric, which is not correct, so we firstly manually change it into charater
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
#Transform array as a factor
myLoad$pd$Array <- as.factor(myLoad$pd$Array)
champ.SVD(beta=myCombat, resultsDir="./CHAMP_SVDimages/")

QC.GUI(beta=myCombat,arraytype = "EPIC")

#DMP
# To adjust the p value to identify = , adjPVal = 0.1
myDMP <- champ.DMP(beta=myCombat,adjPVal = 0.1, arraytype = "EPIC")
DMP.GUI()

DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)
## myDMP is a list now, each data frame is stored as myDMP[[1]], myDMP[[2]], myDMP[[3]]
write_xls(myDMP[[1]], "DMPExcelTable.xlsx")
write.csv(myDMP[[1]],file="./DMP_analysis_result.csv",quote=F,row.names = F)


#DMR
myDMR <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="Bumphunter", arraytype = "EPIC")

DMR.GUI(DMR=myDMR,arraytype="EPIC")
#,compare.group=c("PrEC_cells","LNCaP_cells"))

#GSEA
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
# myDMP and myDMR could (not must) be used directly.

head(myGSEA$DMP)
head(myGSEA$DMR)

##Block 
myBlock <- champ.Block(arraytype = "EPIC")
Block.GUI(arraytype="EPIC")

##myEpiMod <- champ.EpiMod(beta=myNorm,arraytype = "EPIC", resultsDir="./CHAMP_EpiMod/", PDFplot=TRUE)


write_xls(myDMP[[1]], "DMPExcelTable.xlsx")
#################################################################################

#Experiment 1, 2 and 3 all samples with batch effect

library("ChAMP")
library(ggplot2)
library(dplyr)
library(tidyr)
#Local Mac
setwd("~/Analises R/r_array_huvec_120121/IL6_72hrs_vs_PBS_batch")

myLoad <- champ.load(directory = "~/Analises R/r_array_huvec_120121/IL6_72hrs_vs_PBS_batch",arraytype = "EPIC")
myLoad_2 <- champ.load("~/Analises R/r_array_huvec_120121/IL6_72hrs_vs_PBS_batch",method="minfi",arraytype="EPIC") # loaded as minfi way


CpG.GUI( arraytype="EPIC")

#Save specific data
#save(Anno, EPIC.manifest.hg19, multi.hit, myLoad, probe.features, bloodCtl, MatchGeneName, myCNA,
#     myDMP, myNorm, PathwayList, probeInfoALL.lv, file = "champ_huvec_120121_72hrs.RData")

#Save all 

save.image(file = "champ_huvec_120121_72hrs_exp1_2_3_allsamples.RData")
load("champ_huvec_120121_72hrs_exp1_2_3_allsamples.RData")

myLoad$pd

#QC.GUI(beta=myLoad$beta)
#NOrmalization method of BMIQ
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC")
#NOrmalization method of PBC
myNorm_2 <- champ.norm(method="PBC",arraytype = "EPIC")
#NOrmalization method of SWAN
myNorm_3 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="SWAN",arraytype="EPIC")
#NOrmalization method of FunctionalNormalization
myNorm_4 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="FunctionalNormalization",arraytype="EPIC")


#QC.GUI for BMIQ
QC.GUI(beta=myNorm,arraytype = "EPIC")

##QC.GUI for PBC
QC.GUI(beta=myNorm_2,arraytype = "EPIC")

#QC.GUI for SWAN
QC.GUI(beta=myNorm_3,arraytype = "EPIC")

champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)
champ.SVD(beta=myNorm, resultsDir="./CHAMP_SVDimages_withoutnormalization/")

############ Batch Effect Correction
#This function implements the ComBat normalization method for BMIQ
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))
#Slide is pd file is numeric, which is not correct, so we firstly manually change it into charater
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
#Transform array as a factor
myLoad$pd$Array <- as.factor(myLoad$pd$Array)
chamSVD_bmiq <- champ.SVD(beta=myCombat, resultsDir=paste(getwd(), "resultsChamp", sep = "/"))
#principal components (PC-1) correlated to the covariates information provided
QC.GUI(beta=myCombat,arraytype = "EPIC")

write.csv(myCombat,file="./raw_result_bmiq.csv")

#This function implements the ComBat normalization method for PBC
myCombat2 <- champ.runCombat(beta=myNorm_2,pd=myLoad$pd,batchname=c("Slide"))
QC.GUI(beta = myCombat2,arraytype="EPIC" )
write.csv2(myCombat2, file="./Normalized_PBC.csv")
champ.SVD(beta = myCombat2, rgSet=NULL,
                  pd=myLoad$pd,
                  RGEffect=FALSE,
                  PDFplot=TRUE,
                  Rplot=TRUE,
               resultsDir="~/Analises R/r_array_huvec_120121/IL6_72hrs_vs_PBS_batch")

write.csv(champSVD,file="./SVD_analysis_result.csv",quote=F,row.names = F)

#This function implements the ComBat normalization method for SWAN
myCombat3 <- champ.runCombat(beta=myNorm_3,pd=myLoad$pd,batchname=c("Slide"))
QC.GUI(beta = myCombat3,arraytype="EPIC" )

#DMP BMIQ
# To adjust the p value to identify = , adjPVal = 0.1
myDMP <- champ.DMP(beta=myCombat, arraytype = "EPIC")
DMP.GUI()

DMP.GUI(DMP=myDMP[[1]],beta=myCombat,pheno=myLoad$pd$Sample_Group)
## myDMP is a list now, each data frame is stored as myDMP[[1]], myDMP[[2]], myDMP[[3]]
write_xls(myDMP[[1]], "DMPExcelTable.xlsx")
write.csv(myDMP[[1]],file="./DMP_analysis_result.csv",quote=F,row.names = F)

write.csv(myDMP[[1]],file="./DMP_analysis_result.csv")


head(myDMP[[1]])

#DMP for the PBC normalization results
myDMP_2 <- champ.DMP(beta = myCombat2,
                   pheno = myLoad$pd$Sample_Group,
                   compare.group = NULL,
                   adjPVal = 0.05,
                   adjust.method = "BH",
                   arraytype = "EPIC")
myDMP_2 <- champ.DMP(beta=myCombat2, arraytype = "EPIC")

head(myDMP_2[[1]])
write.csv2(myDMP[[1]],"/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results/DMP_850_grp1vsgrp2.csv" )
DMP.GUI(DMP=myDMP_2[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)

# DMP for SWAN normalization results
myDMP_3 <- champ.DMP(beta=myCombat3, arraytype = "EPIC")
DMP.GUI(DMP=myDMP_3[[1]],beta=myCombat3,pheno=myLoad$pd$Sample_Group)


#DMR
myDMR <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="Bumphunter", arraytype = "EPIC")
DMR.GUI(DMR=myDMR,arraytype="EPIC")
#,compare.group=c("PrEC_cells","LNCaP_cells"))

head(myDMR$Bumphunter)


myDMR_2 <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="DMRcate", arraytype = "EPIC", cores=1)


#Need to fix
#myDMR_2 <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="DMRcate",arraytype="EPIC")

myDMR_2 <- champ.DMR(arraytype = "EPIC",method="ProbeLasso", adjPvalProbe=0.1, compare.group=c("PBS","IL6"))

#<- champ.DMR(arraytype = "EPIC",method="DMRcate",cores=1)
#DMR.GUI(DMR=myDMR2,arraytype="EPIC")


#GSEA
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
# myDMP and myDMR could (not must) be used directly.

head(myGSEA$DMP)
head(myGSEA$DMR)

##Block 
myBlock <- champ.Block(beta=myCombat,pheno=myLoad$pd$Sample_Group, arraytype = "EPIC")
Block.GUI(arraytype="EPIC")

##myEpiMod <- champ.EpiMod(beta=myNorm,arraytype = "EPIC", resultsDir="./CHAMP_EpiMod/", PDFplot=TRUE)

myEpiMod <- champ.EpiMod(arraytype="EPIC")


write_xls(myDMP[[1]], "DMPExcelTable.xlsx")

phen<-myLoad$pd[match(colnames(myNorm$beta),myLoad$pd$Sample_Name),]$Sample_Group
HeatMap(myCombat$beta,phen=phen,varselect=1000,plot="heatmap.pdf",cexRow = 0.01,cexCol = 1.2,Colv=T,Rowv=T)


myLoad$pd$Sample_Group
######################################################################################################
CpG.GUI(CpG=rownames(myLoad$beta),
        arraytype="EPIC")

champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)

##CNVAnalysis
myCNA <- champ.CNA(intensity=myLoad$beta,pheno=myLoad$pd$Sample_Group,controlGroup='C',arraytype="EPIC")
#Save the differenceCNVResult
write.csv(myCNA$groupResult,file="./CNV_analysis_result.csv",quote=F,row.names = F)




#####################################################################
champ.SVD(beta=myNorm,pd=myLoad$pd)

#Correct with there are two or more factors
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))

#Methylation probe difference analysis
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group, arraytype = "EPIC")
DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)

#myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="ProbeLasso") 


myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")

#### CONTROLS VS. ALL SAMPLES WITHOUT IL6_6hrs_3

analysis_070921_control_vs_alltreatment_without_il6_6hrs_3

setwd("C:/Users/bossarr/Desktop/r_aray/analysis_070921_control_vs_alltreatment_without_il6_6hrs_3")

save.image(file = "My_Object.RData")

#Load samples

testDir=system.file("extdata",package="ChAMPdata")
myLoad <- champ.load(testDir,arraytype="EPIC")

#display the distribution of CpGs on chromosome, CpG island, TSS reagions.
CpG.GUI(CpG=rownames(myLoad$beta),
        arraytype="EPIC")

#Data Quality Control QC
champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)

# provide five interactive plot
QC.GUI(beta=myLoad$beta,arraytype="EPIC")

#Normalization 
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=5)

#plots for normalization data
QC.GUI(beta=myLoad$beta,arraytype="EPIC")

#Differentially Methylated Probes (DMPs)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group, arraytype = "EPIC")

head(myDMP[[1]])

#GUI interface to check the result of myDMP
DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)
#MVP Overview -  left panel indicates the cut off parameters for significant CpGs
#Heat map - rows and columns of this heatmap has been Hierarchical Clustered

#Differential Methylation Regions - extended segments of the genome that show a quantitative alteration in DNA methylation levels between two groups.

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter", arraytype = "EPIC")
#three DMR algorithms implemented within ChAMP: Bumphunter, ProbeLasso and DMRcate.

head(myDMR$BumphunterDMR)
DMR.GUI(DMR=myDMR)

############### CONTROL VS. 72 Hours ##############
library("ChAMP")
library("knitr")

#Local Mac
setwd("~/Analises R/r_aray")

setwd("~/Analises R/r_aray/contros_vs_72hrs")

myLoad <- champ.load(directory = getwd(),arraytype="EPIC")

save.image(file = "control_vs_72hrs.RData")
load(file = "control_vs_72hrs.Rdata")

#display the distribution of CpGs on chromosome, CpG island, TSS reagions.
CpG.GUI(CpG=rownames(myLoad$beta),
        arraytype="EPIC")

#Data Quality Control QC
champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)

# provide five interactive plot
QC.GUI(beta=myLoad$beta,arraytype="EPIC")

#Normalization 
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC", cores=5)

write.csv(myNorm,file="./Normalization Data.csv",quote=F,row.names = T)

#plots for normalization data
QC.GUI(beta=myLoad$beta,arraytype="EPIC")

#Differentially Methylated Probes (DMPs)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group, arraytype = "EPIC")

head(myDMP[[1]])

write.csv(myDMP[[1]],file="./DMP_analysis_result.csv",quote=F,row.names = F)
write.csv(myDMP[[1]],file="./DMP_analysis_result.csv",quote=F,row.names = F)


#GUI interface to check the result of myDMP
DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)

#Differential Methylation Regions - extended segments of the genome that show a quantitative alteration in DNA methylation levels between two groups.

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter", arraytype = "EPIC")
#three DMR algorithms implemented within ChAMP: Bumphunter, ProbeLasso and DMRcate.

head(myDMR1$BumphunterDMR)

DMR.GUI(DMR=myDMR,
        beta=myNorm,
        pheno=myLoad$pd$Sample_Group,
        runDMP=TRUE,
        compare.group=NULL,
        arraytype="EPIC")


#Differential Methylation Blocks
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")

Block.GUI(Block=myBlock,beta=myNorm,pheno=myLoad$pd$Sample_Group,runDMP=TRUE,compare.group=NULL,arraytype="EPIC")

#Gene Set Enrichment Analysis
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05)

myGSEAtest <- champ.GSEA(beta=myNorm,
           DMP=myDMP[[1]],
           DMR=myDMR,
           CpGlist=NULL,
           Genelist=NULL,
           pheno=myLoad$pd$Sample_Group,
           method="fisher",
           arraytype="EPIC",
           Rplot=TRUE,
           adjPval=0.05,
           cores=1)
myGSEAtest

head(myGSEAtest$DMP)

write.csv(myGSEAtest,file="./enrichment_analysis_result.csv",quote=F,row.names = F)

#Differential Methylated Interaction Hotspots - Not working
#myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group, arraytype="EPIC")

myCNA <- champ.CNA(intensity=myLoad$beta,pheno=myLoad$pd$Sample_Group,sampleCNA = FALSE,controlGroup='C',arraytype="EPIC")


write.csv(myCNA$groupResult,file="./CNV_analysis_result.csv",quote=F,row.names = F)


############### CONTROL VS. 6 Hours ##############
library("ChAMP")
library("knitr")

#Local Mac
setwd("~/Analises R/r_aray")

setwd("~/Analises R/r_aray/contros_vs_6hrs")

myLoad <- champ.load(directory = getwd(),arraytype="EPIC")

save.image(file = "control_vs_6hrs.RData")
load(file = "control_vs_6hrs.Rdata")

#display the distribution of CpGs on chromosome, CpG island, TSS reagions.
CpG.GUI(CpG=rownames(myLoad$beta),
        arraytype="EPIC")

#Data Quality Control QC
champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)

# provide five interactive plot
QC.GUI(beta=myLoad$beta,arraytype="EPIC")

#Normalization 
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC", cores=5)

#plots for normalization data
QC.GUI(beta=myLoad$beta,arraytype="EPIC")

#Differentially Methylated Probes (DMPs)
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group, arraytype = "EPIC")

head(myDMP[[1]])

#GUI interface to check the result of myDMP
DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)


#Differential Methylation Regions - extended segments of the genome that show a quantitative alteration in DNA methylation levels between two groups.

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="DMRcate", arraytype = "EPIC")
#three DMR algorithms implemented within ChAMP: Bumphunter, ProbeLasso and DMRcate.

head(myDMR$BumphunterDMR)

DMR.GUI(DMR=myDMR,
        beta=myNorm,
        pheno=myLoad$pd$Sample_Group,
        arraytype="EPIC")

#Differential Methylation Blocks
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="EPIC")

Block.GUI(Block=myBlock,beta=myNorm,pheno=myLoad$pd$Sample_Group,runDMP=TRUE,compare.group=NULL,arraytype="EPIC")

#Gene Set Enrichment Analysis
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")

myGSEAtest <- champ.GSEA(beta=myNorm,
                         DMP=myDMP[[1]],
                         DMR=myDMR,
                         CpGlist=NULL,
                         Genelist=NULL,
                         pheno=myLoad$pd$Sample_Group,
                         method="fisher",
                         arraytype="EPIC",
                         Rplot=TRUE,
                         adjPval=0.05,
                         cores=1)
myGSEAtest

head(myGSEAtest$DMR)

#Differential Methylated Interaction Hotspots - Not working
#myEpiMod <- champ.EpiMod(beta=myNorm,pheno=myLoad$pd$Sample_Group, arraytype="EPIC")

myCNA <- champ.CNA(intensity=myLoad$beta,pheno=myLoad$pd$Sample_Group,controlGroup='C',arraytype="EPIC")
write.csv(myCNA$groupResult,file="./CNV_analysis_result.csv",quote=F,row.names = F)

##END


