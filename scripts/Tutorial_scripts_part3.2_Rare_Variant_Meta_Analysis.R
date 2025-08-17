################################################################################
# IGES 2025 Education Workshop: Navigating Biobank Data:
# Advanced Strategies for Genetic Epidemiology Research
# Session 4: Practical Tutorial – A Walkthrough of Cross-Biobank Analysis
# Part 3: Meta-analysis of sequencing data using MetaSTAARlite
#   - Part 3.2: Rare variant meta-analysis using MetaSTAARlite
#   - More tutorial: https://github.com/li-lab-genetics/MetaSTAARlite-tutorial
#   - Date: August 2025
################################################################################

rm(list=ls())
gc()

### Setup repository
setwd("/data/williamsjacr/IGES_2025_Education_Workshop/data/")

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(MetaSTAAR)
library(MetaSTAARlite)
library(Matrix)

#####################################################
#           User Input
#####################################################
## Sample sizes of participating studies
sample.sizes <- c(661,347,504,503,489)

## variant_type
variant_type <- "SNV"
## cov_maf_cutoff
cov_maf_cutoff <- c(0.05,0.05,0.05,0.05,0.05)

## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name (for variant weighting using the MetaSTAAR method)
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")


#####################################################
#   Rare variant meta-analysis using MetaSTAARlite
#   Gene-Centric Coding: LDLR
#####################################################
chr <- 19

### run coding mask of LDLR
gene_name <- "LDLR"

## Directories of the study-specific summary statistics file folders
file.dir <- c("","","","","")
file.prefix <- c("AFR_LDLR_coding",
                 "AMR_LDLR_coding",
                 "EAS_LDLR_coding",
                 "EUR_LDLR_coding",
                 "SAS_LDLR_coding")

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat.Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov.Rdata")
coding_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
coding_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)

genes <- genes_info

coding_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
  coding_sumstat_list[[paste0(x,".coding_sumstat")]][[gene_name]]
})
coding_cov_gene_list <- lapply(cov.file.path, function(x) {
  coding_cov_list[[paste0(x,".coding_cov")]][[gene_name]]
})
results_coding_meta <- coding_MetaSTAARlite(chr=chr,gene_name=gene_name,genes=genes,
                                            sample.sizes=sample.sizes,coding_sumstat_gene_list=coding_sumstat_gene_list,
                                            coding_cov_gene_list=coding_cov_gene_list,
                                            cov_maf_cutoff=cov_maf_cutoff,
                                            rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                            variant_type=variant_type,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
results_coding_meta

#####################################################
#   Rare variant meta-analysis using MetaSTAARlite
#   Gene-Centric Noncoding: APOE 
#####################################################
chr <- 19

### run noncoding mask of APOE
gene_name <- "APOE"

## Directories of the study-specific summary statistics file folders
file.dir <- c("","","","","")
file.prefix <- c("AFR_APOE_noncoding",
                 "AMR_APOE_noncoding",
                 "EAS_APOE_noncoding",
                 "EUR_APOE_noncoding",
                 "SAS_APOE_noncoding")

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat.Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov.Rdata")
noncoding_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
noncoding_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)

noncoding_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
  noncoding_sumstat_list[[paste0(x,".noncoding_sumstat")]][[gene_name]]
})
noncoding_cov_gene_list <- lapply(cov.file.path, function(x) {
  noncoding_cov_list[[paste0(x,".noncoding_cov")]][[gene_name]]
})
results_noncoding_meta <- noncoding_MetaSTAARlite(chr=chr,gene_name=gene_name,
                                                  sample.sizes=sample.sizes,noncoding_sumstat_gene_list=noncoding_sumstat_gene_list,
                                                  noncoding_cov_gene_list=noncoding_cov_gene_list,
                                                  cov_maf_cutoff=cov_maf_cutoff,
                                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                  variant_type=variant_type,
                                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
results_noncoding_meta

