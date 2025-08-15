################################################################################
# IGES 2025 Education Workshop: Navigating Biobank Data:
# Advanced Strategies for Genetic Epidemiology Research
# Session 4: Practical Tutorial – A Walkthrough of Cross-Biobank Analysis
# Part 2: Rare variant analysis of sequencing data using STAARpipeline
#   - Part 2.3: WGS rare variant association analysis using STAARpipeline
#   - More tutorial: https://github.com/xihaoli/STAARpipeline-tutorial
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
library(STAARpipelineSummary)

## gds file
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
## Annotation name (for variant weighting using the STAAR method)
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
## Annotation name (for variants info in summary)
Annotation_name_info <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category",
                          "MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF",
                          "aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## known loci: rs7412 (APOE), rs429358 (APOE), rs35136575 (APOE), rs12151108 (LDLR), rs688 (LDLR), rs6511720 (LDLR)
known_loci <- read.csv("known_loci_info.csv")

##########################################
#       load phenotypes and ancestry PCs
##########################################
### phenotype
pheno <- read.csv("phenotype_LDLR_coding_APOE_noncoding.csv")
### PCs
PCs <- read.csv("1000G_PCA.csv")

pheno <- dplyr::left_join(pheno,PCs[,2:12],by=c("sample"="id"))

##########################################
#       fit null model
##########################################
obj_nullmodel <- fit_nullmodel(Y~gender+super_pop+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                               data = pheno, kins = NULL, id = "sample",
                               family = gaussian(link = "identity"), verbose=T)


#####################################################
#        Gene-Centric Coding: LDLR
#####################################################

### run coding mask of LDLR
gene_name <- "LDLR"

## genotype: chr 
chr <- 19
gds.path <- agds_dir
genofile <- seqOpen(gds.path)

results_coding <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                      rare_maf_cutoff=0.01,rv_num_cutoff=2,
								                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

results_coding

## Conditional Analysis
category <- "plof_ds"								
results_coding_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,category=category,
                                                known_loci=known_loci,rare_maf_cutoff=0.01,rv_num_cutoff=2,
								                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
								
results_coding_cond

## variants info for summary
results_coding_info <- Gene_Centric_Coding_Info(category=category,chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,
										QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
										Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name_info)
		
dim(results_coding_info)
# [1] 17 33						
head(results_coding_info)

seqClose(genofile)
								
########################################################
#        Gene-Centric Noncoding: APOE
########################################################

### run noncoding mask of APOE
gene_name <- "APOE"

## genotype: chr 
chr <- 19
gds.path <- agds_dir
genofile <- seqOpen(gds.path)

results_noncoding <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                            rare_maf_cutoff=0.01,rv_num_cutoff=2,
								                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

results_noncoding

## Conditional Analysis
category <- "enhancer_DHS"								
results_noncoding_cond <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,category=category,
                                                      known_loci=known_loci,rare_maf_cutoff=0.01,rv_num_cutoff=2,
								                                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
								                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								                                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

results_noncoding_cond
								
## variants info for summary
results_noncoding_info <- Gene_Centric_Noncoding_Info(category=category,chr=chr,genofile=genofile,obj_nullmodel=obj_nullmodel,gene_name=gene_name,known_loci=known_loci,
										                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
										                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Annotation_name=Annotation_name_info)

dim(results_noncoding_info)
# [1] 50 33						
head(results_noncoding_info)
								
seqClose(genofile)

