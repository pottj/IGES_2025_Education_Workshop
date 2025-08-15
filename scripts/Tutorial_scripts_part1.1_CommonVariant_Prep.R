################################################################################
# IGES 2025 Education Workshop: Navigating Biobank Data:
# Advanced Strategies for Genetic Epidemiology Research
# Session 4: Practical Tutorial – A Walkthrough of Cross-Biobank Analysis
# Part 4: Common Variant PRS using summary statistics 
#   - Part 4.1: Obtaining common variants and setting up summary statistics
#   - Date: August 2025
################################################################################

rm(list=ls())

library(bigsnpr)
library(dplyr)
library(stringr)

gc()

setwd("/data/williamsjacr/IGES_2025_Education_Workshop/data/")

## extract phenotype sample.id
pheno_cov <- read.table("integrated_call_samples_v3.20130502.ALL.panel",header=TRUE)
phenotype.id <- as.vector(pheno_cov$sample)
## write a keep file to filter vcf
write.table(phenotype.id,file = paste0("Keep.txt"),row.names = FALSE,col.names = FALSE,quote = FALSE)
## filter vcf to the 2504 participants and MAF >= 0.05
system("/data/williamsjacr/software/plink2 --vcf 1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --keep Keep.txt --maf 0.05 --make-bed --out chr19")

## create a rds file and bk file with bigsnpr
if(file.exists("chr19.rds")){
  file.remove("chr19.rds")
  file.remove("chr19.bk") 
}
snp_readBed("chr19.bed")

## using the bim file generate create a simulated response and a summary statistics file
chr19_bim <- read.delim("chr19.bim", header=FALSE)

obj.bigSNP <- snp_attach("chr19.rds")
G <- obj.bigSNP$genotypes   # genotype matrix

## create a simulated response
set.seed(1440)
simulated_beta <-  rep(0, length.out = nrow(chr19_bim))
simulated_beta[sample(1:nrow(chr19_bim),5000)] <- rnorm(5000)
simulated_beta <- data.frame(SNP = chr19_bim[,2],pos = chr19_bim[,4], ref = chr19_bim[,5], alt = chr19_bim[,6], beta = simulated_beta,se = runif(nrow(chr19_bim),0.5,5),pval = runif(nrow(chr19_bim),0,1), n = 2504)

Trait_of_Interest <- big_prodVec(G, simulated_beta$beta)
chr19_fam <- read.delim("chr19.fam", header=FALSE)
Phenotype_Data <- data.frame(IID = chr19_fam[,2],Y = Trait_of_Interest)

write.csv(Phenotype_Data,file = "Phenotype_CommonVariant.csv",row.names = FALSE)

## create a fake external biobank summary statistics 

SNPinfo <- read.csv("SNPInfo_GrCH38.csv")

SNPinfo_38 <- data.frame(rsid = SNPinfo$rsid,unique_id1 = paste0(SNPinfo$chr,"_",SNPinfo$pos38,"_",SNPinfo$allele1_38,"_",SNPinfo$allele2_38),
                         unique_id2 = paste0(SNPinfo$chr,"_",SNPinfo$pos38,"_",SNPinfo$allele2_38,"_",SNPinfo$allele1_38))
rm(SNPinfo)

gc()

## use rsids instead of GrCH38

simulated_beta$unique_id <- paste0(19,"_",simulated_beta$pos,"_",simulated_beta$ref,"_",simulated_beta$alt)
sumstats_modified <- left_join(simulated_beta,SNPinfo_38,by = c("unique_id" = "unique_id1"))
sumstats_modified$SNP[!is.na(sumstats_modified$rsid)] <- sumstats_modified$rsid[!is.na(sumstats_modified$rsid)] 
sumstats_modified <- subset(sumstats_modified,select = -c(rsid))
sumstats_modified <- left_join(sumstats_modified,SNPinfo_38,by = c("unique_id" = "unique_id2"))
sumstats_modified$SNP[!is.na(sumstats_modified$rsid)] <- sumstats_modified$rsid[!is.na(sumstats_modified$rsid)] 
sumstats_modified <- sumstats_modified[str_detect(sumstats_modified$SNP,"rs"),]
sumstats_modified <- sumstats_modified[,c("SNP","ref","alt","beta","se","pval","n")]
colnames(sumstats_modified) <- c("RSID","REF","ALT","BETA","SE","PVAL","N")
sumstats_modified <- sumstats_modified[!duplicated(sumstats_modified$RSID),]

## swap some alleles for illustrative purposes

swap_indexes <- sample(1:nrow(sumstats_modified),5000)
sumstats_crossbiobank <- sumstats_modified
sumstats_crossbiobank$ALT[swap_indexes] <- sumstats_modified$REF[swap_indexes]
sumstats_crossbiobank$REF[swap_indexes] <- sumstats_modified$ALT[swap_indexes]
sumstats_crossbiobank$BETA[swap_indexes] <- (-1)*sumstats_modified$BETA[swap_indexes]

write.csv(sumstats_crossbiobank,file = "SumStats_Biobank.csv",row.names = FALSE)

