[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# IGES 2025 Cross Biobank Practical Tutorials
This repository contains code used for the tutorial (`IGES_CrossBioBank_PracticalTutorial.Rmd`) as well as scripts used to create the data and test analyses used in the tutorial (`scripts/`).

- All the source data are publicly available;
- All the source code are open-access;
- All the analysis results are reproducible.

## Prerequisites
The packages used in the tutorial include: 

```R
install.packages("devtools")
install.packages("Matrix")
install.packages("BiocManager")
install.packages('restfulr')
install.packages('logistf')
install.packages('GMMAT')
install.packages("remotes")
install.packages("dplyr")

BiocManager::install('S4Vectors')
BiocManager::install('gdsfmt')
BiocManager::install('SeqArray')
BiocManager::install('SeqVarTools')
BiocManager::install('GENESIS')
BiocManager::install('GenomicFeatures')
BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')

devtools::install_github('xihaoli/STAAR')
devtools::install_github('zilinli1988/SCANG')
devtools::install_github('xihaoli/MultiSTAAR',ref='main')
devtools::install_github('xihaoli/STAARpipeline',ref='main')
devtools::install_github('xihaoli/STAARpipelineSummary',ref='main')
devtools::install_github('xihaoli/MetaSTAAR',ref='main')
devtools::install_github("li-lab-genetics/MetaSTAARlite",ref="main")

remotes::install_github("privefl/bigsnpr")
```

The data used for this tutorial is available at https://dataverse.harvard.edu/api/access/datafile/11945376. 	
This contains a tar.gz file with data used in the IGES 2025 Cross Biobank Practical Tutorial. The vcf used in this tutorial to create common variants and agds files is from 1000G (https://doi.org/10.1016/j.cell.2022.08.004) and downloadable at https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/. Unrelated individuals were identified using the panel downloadable from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/. All other files were created or simulated using the code in the `scripts/` folder.

```R
download.file(url = "https://dataverse.harvard.edu/api/access/datafile/11945376",
              destfile = "")
```

## Tutorial Vignette
Please see the <a href="https://htmlpreview.github.io/?https://github.com/li-lab-genetics/IGES_2025_Education_Workshop/blob/main/IGES_CrossBioBank_PracticalTutorial.html">**IGES Cross Biobank Practical Tutorial**</a> for a step-by-step tutorial of the <a href="https://www.geneticepi.org/2025-education-workshop">2025 IGES Education Workshop</a> in Cologne, Germany.

## Session Info
```R
> sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 15.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] MetaSTAARlite_0.9.7   bigsnpr_1.12.21       bigstatsr_1.6.2       STAARpipeline_0.9.7.2

loaded via a namespace (and not attached):
  [1] rstudioapi_0.17.1                        jsonlite_1.8.9                          
  [3] shape_1.4.6.1                            magrittr_2.0.3                          
  [5] jomo_2.7-6                               GenomicFeatures_1.58.0                  
  [7] logistf_1.26.0                           nloptr_2.1.1                            
  [9] rmarkdown_2.29                           BiocIO_1.16.0                           
 [11] zlibbioc_1.52.0                          vctrs_0.6.5                             
 [13] memoise_2.0.1                            minqa_1.2.8                             
 [15] Rsamtools_2.22.0                         SNPRelate_1.40.0                        
 [17] RCurl_1.98-1.16                          quantsmooth_1.72.0                      
 [19] CompQuadForm_1.4.3                       htmltools_0.5.8.1                       
 [21] S4Arrays_1.6.0                           curl_6.0.1                              
 [23] broom_1.0.7                              SparseArray_1.6.0                       
 [25] mitml_0.4-5                              desc_1.4.3                              
 [27] sandwich_3.1-1                           zoo_1.8-12                              
 [29] cachem_1.1.0                             GenomicAlignments_1.42.0                
 [31] lifecycle_1.0.4                          iterators_1.0.14                        
 [33] pkgconfig_2.0.3                          Matrix_1.7-1                            
 [35] R6_2.6.1                                 fastmap_1.2.0                           
 [37] GenomeInfoDbData_1.2.13                  MatrixGenerics_1.18.0                   
 [39] digest_0.6.37                            colorspace_2.1-1                        
 [41] SCANG_1.0.4                              GWASTools_1.52.0                        
 [43] AnnotationDbi_1.68.0                     S4Vectors_0.44.0                        
 [45] ps_1.8.1                                 rmio_0.4.0                              
 [47] GenomicRanges_1.58.0                     RSQLite_2.3.9                           
 [49] httr_1.4.7                               abind_1.4-8                             
 [51] mgcv_1.9-1                               compiler_4.4.1                          
 [53] rngtools_1.5.2                           remotes_2.5.0                           
 [55] doParallel_1.0.17                        bit64_4.5.2                             
 [57] backports_1.5.0                          BiocParallel_1.40.0                     
 [59] DBI_1.2.3                                MultiSTAAR_0.9.7.1                      
 [61] bigassertr_0.1.7                         pkgbuild_1.4.5                          
 [63] pan_1.9                                  MASS_7.3-61                             
 [65] quantreg_5.99.1                          STAAR_0.9.7.2                           
 [67] DelayedArray_0.32.0                      rjson_0.2.23                            
 [69] DNAcopy_1.80.0                           gdsfmt_1.42.0                           
 [71] tools_4.4.1                              lmtest_0.9-40                           
 [73] GMMAT_1.4.2                              nnet_7.3-19                             
 [75] glue_1.8.0                               quadprog_1.5-8                          
 [77] restfulr_0.0.15                          callr_3.7.6                             
 [79] nlme_3.1-166                             bigparallelr_0.3.2                      
 [81] grid_4.4.1                               GENESIS_2.36.0                          
 [83] generics_0.1.4                           operator.tools_1.6.3                    
 [85] gtable_0.3.6                             bigsparser_0.7.3                        
 [87] flock_0.7                                formula.tools_1.7.1                     
 [89] tidyr_1.3.1                              data.table_1.16.4                       
 [91] XVector_0.46.0                           BiocGenerics_0.52.0                     
 [93] GWASExactHW_1.2                          foreach_1.5.2                           
 [95] pillar_1.11.0                            splines_4.4.1                           
 [97] dplyr_1.1.4                              lattice_0.22-6                          
 [99] survival_3.7-0                           rtracklayer_1.66.0                      
[101] bit_4.5.0.1                              SparseM_1.84-2                          
[103] tidyselect_1.2.1                         SeqVarTools_1.44.0                      
[105] kinship2_1.9.6.1                         Biostrings_2.74.0                       
[107] knitr_1.49                               IRanges_2.40.1                          
[109] SummarizedExperiment_1.36.0              stats4_4.4.1                            
[111] xfun_0.49                                expm_1.0-0                              
[113] Biobase_2.66.0                           matrixStats_1.4.1                       
[115] UCSC.utils_1.2.0                         yaml_2.3.10                             
[117] boot_1.3-31                              TxDb.Hsapiens.UCSC.hg38.knownGene_3.20.0
[119] evaluate_1.0.1                           codetools_0.2-20                        
[121] tibble_3.3.0                             cli_3.6.5                               
[123] rpart_4.1.23                             munsell_0.5.1                           
[125] processx_3.8.4                           MetaSTAAR_0.9.6.3                       
[127] Rcpp_1.1.0                               GenomeInfoDb_1.42.1                     
[129] png_0.1-8                                XML_3.99-0.17                           
[131] parallel_4.4.1                           MatrixModels_0.5-3                      
[133] ggplot2_3.5.2                            blob_1.2.4                              
[135] STAARpipelineSummary_0.9.7.1             doRNG_1.8.6.2                           
[137] bitops_1.0-9                             lme4_1.1-35.5                           
[139] glmnet_4.1-8                             SeqArray_1.46.0                         
[141] scales_1.3.0                             purrr_1.0.2                             
[143] crayon_1.5.3                             rlang_1.1.6                             
[145] cowplot_1.2.0                            KEGGREST_1.46.0                         
[147] mice_3.17.0
```

## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
