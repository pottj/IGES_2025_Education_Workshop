################################################################################
# IGES 2025 Education Workshop: Navigating Biobank Data:
# Advanced Strategies for Genetic Epidemiology Research
# Session 4: Practical Tutorial – A Walkthrough of Cross-Biobank Analysis
# Part 4: Common Variant Meta-Analysis using Metal
#   - Part 4.3: Creating summary statistics files to use with METAL
#   - Date: August 2025
################################################################################

rm(list=ls())

library(bigsnpr)
library(dplyr)

gc()

setwd("/data/williamsjacr/IGES_2025_Education_Workshop/data/")

SumStats_Biobank1 <- read.csv("SumStats_Biobank.csv")

SNPinfo <- read.csv("SNPInfo_GrCH38.csv")
SNPinfo <- data.frame(rsid = SNPinfo$rsid,chr = SNPinfo$chr, pos = SNPinfo$pos38)
SumStats_Biobank1 <- left_join(SumStats_Biobank1,SNPinfo, by = c("RSID" = "rsid"))

SumStats_Biobank1 <- data.frame(SNP = SumStats_Biobank1$RSID, CHR = SumStats_Biobank1$chr, BP = SumStats_Biobank1$pos, 
                                EA = SumStats_Biobank1$ALT,NEA = SumStats_Biobank1$REF, BETA = SumStats_Biobank1$BETA,
                                SE = SumStats_Biobank1$SE, P = SumStats_Biobank1$PVAL, N = SumStats_Biobank1$N)

SumStats_Biobank2 <- SumStats_Biobank1[sample(1:nrow(SumStats_Biobank1),round(0.8*nrow(SumStats_Biobank1))),]

write.table(SumStats_Biobank1,file = "SumStats_Biobank1.txt",sep = "\t",quote = FALSE,row.names = FALSE)
write.table(SumStats_Biobank2,file = "SumStats_Biobank2.txt",sep = "\t",quote = FALSE,row.names = FALSE)

metal_script <- c(
  "SCHEME STDERR",
  "MARKER SNP",
  "ALLELE EA NEA",
  "EFFECT BETA",
  "STDERR SE",
  "PVAL P",
  "WEIGHT N",
  "",
  "PROCESS SumStats_Biobank1.txt",
  "PROCESS SumStats_Biobank2.txt",
  "",
  "OUTFILE meta_results_ .tbl",
  "ANALYZE"
)

writeLines(metal_script, "metal_script.txt")

# system("./metal metal_script.txt")