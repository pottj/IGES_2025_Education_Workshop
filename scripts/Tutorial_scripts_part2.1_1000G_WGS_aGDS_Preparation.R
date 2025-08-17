################################################################################
# IGES 2025 Education Workshop: Navigating Biobank Data:
# Advanced Strategies for Genetic Epidemiology Research
# Session 4: Practical Tutorial – A Walkthrough of Cross-Biobank Analysis
# Part 2: Rare variant analysis of sequencing data using STAARpipeline
#   - Part 2.1: Preparing 1000 Genome WGS annotated GDS files
#   - Date: August 2025
################################################################################

### Setup repository
setwd("/data/williamsjacr/IGES_2025_Education_Workshop/data/")
getOption('timeout')
options(timeout=200)

### Download 1000 Genomes Project (1kGP) high-coverage Illumina integrated phased panel
# ETA: 1 minute
download.file(url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
              destfile = "1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

### Convert 1000G WGS Data VCF files to Genomic Data Structure (GDS) files
library(gdsfmt)
library(SeqArray)

# input_dir
vcf.fn <- "1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

# output_dir
out.fn <- "1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds"

# modify the header 
h <- seqVCF_Header(vcf.fn)
# h$info
h$info$Number[h$info$ID=="SOURCE"] <- "."

# ETA: 10 minutes
seqVCF2GDS(vcf.fn, out.fn, header = h, genotype.var.name = "GT", info.import=NULL, fmt.import=NULL, 
           ignore.chr.prefix="chr", raise.error=TRUE, verbose=TRUE)

### Functionally annotate 1000G WGS Data GDS files to annotated GDS (aGDS) files
### using FAVOR database (https://favor.genohub.org/) and FAVORannotator
### The following scripts are used for the favorannotator app in DNAnexus
### There are also offline scripts for favorannotator available for use:
### (https://github.com/xihaoli/STAARpipeline-Tutorial#generate-annotated-gds-agds-file-using-favorannotator)

# Upload 1000G WGS GDS files to DNAnexus RAP Cloud (Drag the files)

#-------------------------------------------------------------------
# The following are dx-toolkit command line scripts, not R scripts
# Install DNAnexus Platform SDK (MacOS)
tar -xzf dx-toolkit-0.398.0.tar.gz
source dx-toolkit-0.398.0/environment

##### Usage
cd dx-toolkit-0.398.0/bin
#dx upgrade

##### Log in
dx login
#dx ls (here we assume to have a project named ukbb_lilab)

# Clone this github repo to some directory:
git clone https://github.com/li-lab-genetics/favorannotator-rap.git

# Navigate to a relevant directory within the project directory on the DNAnexus platform
dx cd UKB_PRS:/
  
# Compile the source code:
dx build -f favorannotator-rap

# Create a new folder (UKB_PRS/IGES_2025_Education_Workshop/1000G/aGDS)
# on DNAnexus under the project directory
# Run favorannotator on DNAnexus (ETA: 26 minutes)
dx run UKB_PRS:/favorannotator \
-igds_file=UKB_PRS:/IGES_2025_Education_Workshop/1000G/GDS/1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.gds \
-ichromosome=19 \
-iuse_compression=YES \
-ioutfile=1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel_annotated \
--destination=UKB_PRS:/IGES_2025_Education_Workshop/1000G/aGDS --yes
#-------------------------------------------------------------------

### Adds QC_label with all "PASS" to a post-QC aGDS file (R scripts)
gds.path <- "1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel_annotated.gds"
genofile <- seqOpen(gds.path, readonly = FALSE)
#genofile
position <- as.integer(seqGetData(genofile, "position"))
length(position)
Anno.folder <- index.gdsn(genofile, "annotation/info")
add.gdsn(Anno.folder, "QC_label", val=factor(rep("PASS", length(position))), compress="LZMA_ra", closezip=TRUE)
#genofile
seqClose(genofile)





