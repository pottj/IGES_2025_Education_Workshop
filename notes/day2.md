# Day 2 (01/09/2025)

## Opening Session

- slido questionaires
- IGES action plan for LLM and IGES

## Session 1: Williams Awards (9am - 10am)

### Yixin Zhang: Biological Group-Guided Mediation Analysis (abstract 97)

- BMI on HbA1c via metabolites
- Classic approach: test all metabolites
- Here: group metabolites to increase power & biological interpretation
- Framingham Heart Study data: BMI at exam 1, metabolites at exam 2, HbA1c at exam 3 to prevent inverse causation

### Alec McKinlay: Influence of Heritable Covariates on Genetic Studies of the Human Gut Microbiome (abstract 18) 

- ... and Consequences for downstream MR
- lack of replication of microbiome QTL & complex mechanisms of association
- example about lactose: if you are intolerant, because _LCT_ is not working, but continue to drink milk, your microbiome will start digesting lactose for you! environmental effect
- Flemish Gut Flora Project (FGFP)
- conclusion: do not do MR with microbiome as exposure (maybe as an outcome -> do always reverse MR) 

### Jakob Woerner: Polygenic Scores in Proteomic Risk Prediction (abstract 43) 

- UKB proteomics 2027 release for all participants & up to 5000 proteins (pharma release 2026)
- prevalent cases (before baseline), incident cases (after baseline), controls
- train set: prevalent + controls, test set: incident + controls
- in general: protein risk scores perform better than the simple polygenetic risk scores

### Fiona Pantrin: Characterisation of Diverse Global Ancestries Within Participants of the UK Biobank (abstract 58) 

- Most researcher focus on "White British" participants of the UKB - the rest is understudied
- get ancestry label for each individual
- try fine-scale ancestry group by IBD
- outlier: high IBD - Ashkenazi Jewish reference

## Session 2: Neel Awards (10.20am - 11:50am)

### Keynote Talk - Claudia Langenberg: Multiomics Insights into Diverse Human Diseases

- genetic architecture of plasma proteins
- start with some old example
- Sex differences in genetic regulation? --> mainly no, effect difference statistical significance, but not clinical relevance, only a handful of proteins with specific effects ([publication](https://www.nature.com/articles/s41467-025-59034-4) --> check these few for causal effect of hormones?)
- Olink & UKB: scale up will be relevant for disease prediction & prognosis! Focus on trans, because that is what will increase with sample size, not cis (already saturated)
- find the pleiotropic common genes
- Genes & Health Study (>55k, David van Heel)
- Try to meet and talk about position in Berlin?

### Robin Hofmeister: Parent-of-Origin Effects on Complex Traits: Evidence from 265,000 Individuals (abstract 26)

- relevant for heterozygous cases, where it matters if the effect allele comes from mom or dad (dad wants to optimize survival of this baby, mom wants to optimize for future babies)
- assign relative parental side by chromosome X
- problem: siblings, solution leverage sex-specific recombination events
- UKB and Estonian Biobank
- 30 POEs, with different effect direction 

### Siru Wang: Heterogeneity Due to Ancestry and Environment Improves the Resolution of Multi-Ancestry Fine-Mapping (abstract 69) 

- MR-MEGA recap
- similar to talk she gave at the CM theme meeting a while back
- too long, I think no one controls the speakers time here...

### Hyunkyu Lee: covImpute: Leveraging Genetic Correlations to Impute Missing EHR Phenotypes (abstract 21) 

- check publication - cannot follow all those equations on the slide (not yet published, maybe some preprint in the next month?)
- new imputation software, better than the others (what all developers say...)

### Cathal Ormond: Bayesian Inference Model to Prioritise Rare Variants from Family-Based Whole Genome Sequencing Data (abstract 28) 

- BICEP: Bayesian Inference for Causality Evaluation in Pedigrees: pathogenicity + co-segregation
- too tired & hungry - hard to follow

## Session 3: Diverstiy (1:25pm - 2:55pm)

### Keynote Talk - Eleftheria Zeggini: Translational Genomics of Complex Disease

- distracted because of CM away day planning chaos

### Brian L. Browning: Estimating Gene Conversion Rates (abstract 15) 

- IBD again (haven't heard about IBD so many times for long)
- data used: TopMed (more samples, higher resolution) & deCODE (possible to stratify & identify the direction of conversion) & UKB

### Peyton McClelland: WGS-Based HLA Allele Imputation Quantifies Pleiotropic Associations with Disease-Related Traits (abstract 110) 

- French Canadian, CARTaGENE (n=30,000 with genotypes)
- IMGT HLA allele database

### Mohamad Saad: High Consanguinity, Underrepresented Populations, and Genome Sequencing for Rare Variant Discoveries (abstract 109) 

- distracted by mails
- everyone under lipid-lowering medication

### Valentina Rukin: Multi-Ancestral GWAS of GDM and Glycemic Traits During Pregnancy (abstract 57) 

- missed beginning, don't know which data was used

## Poster Session 1 (2:55pm - 4:15pm)

## Session 4: Statistical Modeling (4:15pm - 5:45pm)

### Gillian King: Leveraging a Multi-Population Likelihood Framework for Bayesian Model Uncertainty in PRS Construction (abstract 113) 

- Bayesian!
- Simulation with 500 SNPs, 10 causal SNPs, and h2=0.1

### Sebastian Sendel: TrACES of Time: A Targeted mRNA Sequencing Approach for Estimating Time-of-Day of Bloodstain Deposition in Forensic Casework (abstract 86) 

- panel of 69 markers + prediction model
- time as continuous trait: use combination of sin and cos of the 24h clock
- works, but high SE, so at the moment not good to indicate time, but maybe to exclude time

### Chris Wallace: Why Meta-Analysis Fine-Mapping Can Be Confidently Wrong and How to Fix It (abstract 27) 

- nifty - also kind of the CM talk from a month ago

### Chuan Fu Yap: Framework for Allelic Effect Heterogeneity Assessment in Genome-Wide Association Study Meta-Analyses (abstract 45) 

- many equations - still not sure what the goal is, GxE with E being the PCs?
- PC scatter plots - with gradients - PC1-3 are important for continental differences

### Merli Koitmäe: Uncovering Genetic Pathways in Type 2 Diabetes Using Genomic Structural Equation Modeling (abstract 85) 

- Genomic SEM - been a while ...
- SEM needs same ancestry
- many covariables - merge into latent factors

### Julien St-Pierre: Penalized Generalized Linear Mixed Models for Longitudinal Outcomes in Genetic Association Studies (abstract 60) 

- GMMAT as starting point, but it is only quasi-penalized - problem when want to predict stuff
- so lets add random effects per measurement by sample
- "two random slopes" - did not really understand - hope he has a paper/preprint
