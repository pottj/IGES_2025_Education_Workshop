################################################################################
# IGES 2025 Education Workshop: Navigating Biobank Data:
# Advanced Strategies for Genetic Epidemiology Research
# Session 4: Practical Tutorial – A Walkthrough of Cross-Biobank Analysis
# Part 4: Common Variant PRS using different summary statistics 
#   - Part 4.2: Flipping summary statistics and testing if it worked
#   - Date: August 2025
################################################################################

rm(list=ls())

library(bigsnpr)
library(dplyr)

gc()

setwd("/data/williamsjacr/IGES_2025_Education_Workshop/data/")

## create a rds file and bk file with bigsnpr
if(file.exists("chr19.rds")){
  file.remove("chr19.rds")
  file.remove("chr19.bk") 
}
snp_readBed("chr19.bed")

SumStats_CrossBiobank <- read.csv("SumStats_Biobank.csv")

Phenotype_CommonVariant <- read.csv("Phenotype_CommonVariant.csv")

obj.bigSNP <- snp_attach("chr19.rds")
G <- obj.bigSNP$genotypes   # genotype matrix

## Fails as the files are not the same size
## In plink this would fail as the ids are not the same
# Test <- big_prodVec(G, SumStats_CrossBiobank$beta)

## Use the bim file to generate a NULL summary statistics file
chr19_bim <- read.delim("chr19.bim", header=FALSE)
sumstats <- data.frame(SNP = chr19_bim[,2],pos = chr19_bim[,4], ref = chr19_bim[,5], alt = chr19_bim[,6], beta = 0)

## need to join the biobanks summary statistics with our null summary file

## first add chr and pos
SNPinfo <- read.csv("SNPInfo_GrCH38.csv")
SNPinfo <- data.frame(rsid = SNPinfo$rsid,chr = SNPinfo$chr, pos = SNPinfo$pos38)

SumStats_CrossBiobank <- left_join(SumStats_CrossBiobank,SNPinfo, by = c("RSID" = "rsid"))

## Next create a unique id to merge the null file and the biobank sumstats
sumstats$unique_id <- paste0(19,"_",sumstats$pos,"_",sumstats$ref,"_",sumstats$alt)
SumStats_CrossBiobank$unique_id1 <- paste0(19,"_",SumStats_CrossBiobank$pos,"_",SumStats_CrossBiobank$REF,"_",SumStats_CrossBiobank$ALT)
SumStats_CrossBiobank$unique_id2 <- paste0(19,"_",SumStats_CrossBiobank$pos,"_",SumStats_CrossBiobank$ALT,"_",SumStats_CrossBiobank$REF)

## Merge by the unique CHR_POS_REF_ALT ID to add the estimated coefficients to our NULL file
sumstats <- left_join(sumstats,SumStats_CrossBiobank[,c("unique_id1","BETA")],by = c("unique_id" = "unique_id1"))
sumstats$beta[!is.na(sumstats$BETA)] <- sumstats$BETA[!is.na(sumstats$BETA)]
sumstats <- subset(sumstats,select = -c(BETA))
sumstats <- left_join(sumstats,SumStats_CrossBiobank[,c("unique_id2","BETA")],by = c("unique_id" = "unique_id2"))
sumstats$beta[!is.na(sumstats$BETA)] <- sumstats$BETA[!is.na(sumstats$BETA)]
sumstats <- subset(sumstats,select = -c(BETA))

## Calculate the PRS
Test <- big_prodVec(G, sumstats$beta)
all.equal(Test,Phenotype_CommonVariant$Y)
cor(Test,Phenotype_CommonVariant$Y)

## However, the PRS above did not account for the mismatched alleles, if the reference allele is different (target data vs training data) we also have to flip beta's
sumstats$beta[sumstats$unique_id %in% SumStats_CrossBiobank$unique_id2] <- (-1)*sumstats$beta[sumstats$unique_id %in% SumStats_CrossBiobank$unique_id2]

## This improves are accuracy with our response
Test <- big_prodVec(G, sumstats$beta)
all.equal(Test,Phenotype_CommonVariant$Y)
cor(Test,Phenotype_CommonVariant$Y)

write.csv(sumstats,file = "SumStats_Aligned.csv", row.names = FALSE)
