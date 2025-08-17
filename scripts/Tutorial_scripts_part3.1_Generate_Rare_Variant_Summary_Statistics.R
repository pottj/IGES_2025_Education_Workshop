################################################################################
# IGES 2025 Education Workshop: Navigating Biobank Data:
# Advanced Strategies for Genetic Epidemiology Research
# Session 4: Practical Tutorial – A Walkthrough of Cross-Biobank Analysis
# Part 3: Meta-analysis of sequencing data using MetaSTAARlite
#   - Part 3.1: Generate rare variant summary statistics for each study
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

## agds file
agds_dir <- "1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel_annotated.gds"
## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- read.csv(url("https://raw.githubusercontent.com/xihaoli/STAARpipeline-Tutorial/refs/heads/main/FAVORannotator_csv/Annotation_name_catalog.csv"))

## channel name of the QC label in the GDS/aGDS file
QC_label <- "annotation/info/QC_label"
## variant type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## genes info
genes <- genes_info
## extract phenotype sample.id
pheno_cov <- read.table("integrated_call_samples_v3.20130502.ALL.panel",header=TRUE)
phenotype.id <- as.vector(pheno_cov$sample)
length(phenotype.id)
# [1] 2504

## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name (for variant weighting using the MetaSTAAR method)
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

##########################################
#       load phenotype and ancestry PCs
##########################################
### phenotype
pheno <- read.csv("phenotype_LDLR_coding_APOE_noncoding.csv")
### PCs
PCs <- read.csv("1000G_PCA.csv")

pheno <- dplyr::left_join(pheno,PCs[,2:12],by=c("sample"="id"))

##########################################
#       fit null model for each "study"
##########################################
obj_nullmodel_AFR <- fit_nullmodel(Y~gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   data = pheno[pheno$super_pop=="AFR",], kins = NULL, id = "sample",
                                   family = gaussian(link = "identity"), verbose=T)
obj_nullmodel_AMR <- fit_nullmodel(Y~gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   data = pheno[pheno$super_pop=="AMR",], kins = NULL, id = "sample",
                                   family = gaussian(link = "identity"), verbose=T)
obj_nullmodel_EAS <- fit_nullmodel(Y~gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   data = pheno[pheno$super_pop=="EAS",], kins = NULL, id = "sample",
                                   family = gaussian(link = "identity"), verbose=T)
obj_nullmodel_EUR <- fit_nullmodel(Y~gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   data = pheno[pheno$super_pop=="EUR",], kins = NULL, id = "sample",
                                   family = gaussian(link = "identity"), verbose=T)
obj_nullmodel_SAS <- fit_nullmodel(Y~gender+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                                   data = pheno[pheno$super_pop=="SAS",], kins = NULL, id = "sample",
                                   family = gaussian(link = "identity"), verbose=T)


#####################################################
#        For each study, generate rare variant
#        summary statistics using MetaSTAARlite worker
#        Gene-Centric Coding: LDLR
#####################################################

### run coding mask of LDLR
gene_name <- "LDLR"

## genotype: chr 
chr <- 19
gds.path <- agds_dir
genofile <- seqOpen(gds.path)

genes <- genes_info

## AFR
coding_sumstat <- list()
coding_cov <- list()
results_coding_AFR <- coding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_AFR,gene=genes,
                                                  cov_maf_cutoff=0.05,
                                                  QC_label=QC_label,variant_type=variant_type,
                                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

coding_sumstat[[gene_name]] <- results_coding_AFR$summary_stat_list
coding_cov[[gene_name]] <- results_coding_AFR$GTSinvG_rare_list

save(coding_sumstat,file="AFR_LDLR_coding_sumstat.Rdata",compress = "xz")
save(coding_cov,file="AFR_LDLR_coding_cov.Rdata",compress = "xz")

## AMR
coding_sumstat <- list()
coding_cov <- list()
results_coding_AMR <- coding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_AMR,gene=genes,
                                                  cov_maf_cutoff=0.05,
                                                  QC_label=QC_label,variant_type=variant_type,
                                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

coding_sumstat[[gene_name]] <- results_coding_AMR$summary_stat_list
coding_cov[[gene_name]] <- results_coding_AMR$GTSinvG_rare_list

save(coding_sumstat,file="AMR_LDLR_coding_sumstat.Rdata",compress = "xz")
save(coding_cov,file="AMR_LDLR_coding_cov.Rdata",compress = "xz")

## EAS
coding_sumstat <- list()
coding_cov <- list()
results_coding_EAS <- coding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_EAS,gene=genes,
                                                  cov_maf_cutoff=0.05,
                                                  QC_label=QC_label,variant_type=variant_type,
                                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

coding_sumstat[[gene_name]] <- results_coding_EAS$summary_stat_list
coding_cov[[gene_name]] <- results_coding_EAS$GTSinvG_rare_list

save(coding_sumstat,file="EAS_LDLR_coding_sumstat.Rdata",compress = "xz")
save(coding_cov,file="EAS_LDLR_coding_cov.Rdata",compress = "xz")

## EUR
coding_sumstat <- list()
coding_cov <- list()
results_coding_EUR <- coding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_EUR,gene=genes,
                                                  cov_maf_cutoff=0.05,
                                                  QC_label=QC_label,variant_type=variant_type,
                                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

coding_sumstat[[gene_name]] <- results_coding_EUR$summary_stat_list
coding_cov[[gene_name]] <- results_coding_EUR$GTSinvG_rare_list

save(coding_sumstat,file="EUR_LDLR_coding_sumstat.Rdata",compress = "xz")
save(coding_cov,file="EUR_LDLR_coding_cov.Rdata",compress = "xz")

## SAS
coding_sumstat <- list()
coding_cov <- list()
results_coding_SAS <- coding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_SAS,gene=genes,
                                                  cov_maf_cutoff=0.05,
                                                  QC_label=QC_label,variant_type=variant_type,
                                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

coding_sumstat[[gene_name]] <- results_coding_SAS$summary_stat_list
coding_cov[[gene_name]] <- results_coding_SAS$GTSinvG_rare_list

save(coding_sumstat,file="SAS_LDLR_coding_sumstat.Rdata",compress = "xz")
save(coding_cov,file="SAS_LDLR_coding_cov.Rdata",compress = "xz")

seqClose(genofile)

########################################################
#.       For each study, generate rare variant
#        summary statistics using MetaSTAARlite worker
#        Gene-Centric Noncoding: APOE 
########################################################

### run noncoding mask of APOE
gene_name <- "APOE"

## genotype: chr 
chr <- 19
gds.path <- agds_dir
genofile <- seqOpen(gds.path)

## AFR
noncoding_sumstat <- list()
noncoding_cov <- list()
results_noncoding_AFR <- noncoding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_AFR,
                                                        cov_maf_cutoff=0.05,
                                                        QC_label=QC_label,variant_type=variant_type,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

noncoding_sumstat[[gene_name]] <- results_noncoding_AFR$summary_stat_list
noncoding_cov[[gene_name]] <- results_noncoding_AFR$GTSinvG_rare_list

save(noncoding_sumstat,file="AFR_APOE_noncoding_sumstat.Rdata",compress = "xz")
save(noncoding_cov,file="AFR_APOE_noncoding_cov.Rdata",compress = "xz")

## AMR
noncoding_sumstat <- list()
noncoding_cov <- list()
results_noncoding_AMR <- noncoding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_AMR,
                                                        cov_maf_cutoff=0.05,
                                                        QC_label=QC_label,variant_type=variant_type,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

noncoding_sumstat[[gene_name]] <- results_noncoding_AMR$summary_stat_list
noncoding_cov[[gene_name]] <- results_noncoding_AMR$GTSinvG_rare_list

save(noncoding_sumstat,file="AMR_APOE_noncoding_sumstat.Rdata",compress = "xz")
save(noncoding_cov,file="AMR_APOE_noncoding_cov.Rdata",compress = "xz")

## EAS
noncoding_sumstat <- list()
noncoding_cov <- list()
results_noncoding_EAS <- noncoding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_EAS,
                                                        cov_maf_cutoff=0.05,
                                                        QC_label=QC_label,variant_type=variant_type,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

noncoding_sumstat[[gene_name]] <- results_noncoding_EAS$summary_stat_list
noncoding_cov[[gene_name]] <- results_noncoding_EAS$GTSinvG_rare_list

save(noncoding_sumstat,file="EAS_APOE_noncoding_sumstat.Rdata",compress = "xz")
save(noncoding_cov,file="EAS_APOE_noncoding_cov.Rdata",compress = "xz")

## EUR
noncoding_sumstat <- list()
noncoding_cov <- list()
results_noncoding_EUR <- noncoding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_EUR,
                                                        cov_maf_cutoff=0.05,
                                                        QC_label=QC_label,variant_type=variant_type,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

noncoding_sumstat[[gene_name]] <- results_noncoding_EUR$summary_stat_list
noncoding_cov[[gene_name]] <- results_noncoding_EUR$GTSinvG_rare_list

save(noncoding_sumstat,file="EUR_APOE_noncoding_sumstat.Rdata",compress = "xz")
save(noncoding_cov,file="EUR_APOE_noncoding_cov.Rdata",compress = "xz")

## SAS
noncoding_sumstat <- list()
noncoding_cov <- list()
results_noncoding_SAS <- noncoding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel_SAS,
                                                        cov_maf_cutoff=0.05,
                                                        QC_label=QC_label,variant_type=variant_type,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

noncoding_sumstat[[gene_name]] <- results_noncoding_SAS$summary_stat_list
noncoding_cov[[gene_name]] <- results_noncoding_SAS$GTSinvG_rare_list

save(noncoding_sumstat,file="SAS_APOE_noncoding_sumstat.Rdata",compress = "xz")
save(noncoding_cov,file="SAS_APOE_noncoding_cov.Rdata",compress = "xz")

seqClose(genofile)

