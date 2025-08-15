[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# IGES 2025 Cross Biobank Practical Tutorials
This repository contains code used the IGES 2025 cross biobank practical tutorial. 

## Prerequisites
This tutorial is presented in the IGES_CrossBioBank_PracticalTutorial.Rmd R markdown script. The packages used in the vignette include: 

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

The data used for this tutorial is available at https://dataverse.harvard.edu/api/access/datafile/11912216 

```R
download.file(url = "https://dataverse.harvard.edu/api/access/datafile/11912216",
              destfile = "")
```

## License
This software is licensed under GPLv3.

![GPLv3](http://www.gnu.org/graphics/gplv3-127x51.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)
