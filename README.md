[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# IGES 2025 Cross Biobank Practical Tutorials
This repository contains code used for the tutorial (`IGES_CrossBioBank_PracticalTutorial.Rmd`) as well as scripts used to create the data and test analyses used in the tutorial (`scripts/``).

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

The data used for this tutorial is available at https://dataverse.harvard.edu/api/access/datafile/11912614. 	
This contains a tar.gz file with data used in the IGES 2025 Cross Biobank Practical Tutorial. The vcf used in this tutorial to create common variants and agds files is from 1000G (https://doi.org/10.1016/j.cell.2022.08.004) and downloadable at https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/. Unrelated individuals were identified using the panel downloadable from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/. All other files were created or simulated using the code in the scripts folder.

```R
download.file(url = "https://dataverse.harvard.edu/api/access/datafile/11912614",
              destfile = "")
```

## Tutorial Vignette
Please see the <a href="https://htmlpreview.github.io/?https://github.com/li-lab-genetics/IGES_2025_Education_Workshop/blob/main/IGES_CrossBioBank_PracticalTutorial.html">**IGES Cross Biobank Practical Tutorial**</a> for a step-by-step tutorial of the 2025 IGES Education Workshop.

## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
