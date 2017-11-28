SNV selection intensity
================
Vincent L. Cannataro

This [R Markdown](http://rmarkdown.rstudio.com) Notebook, output as a `github_document`, describes the process by which we estimate the selection intensity for single nucleotide variants (SNVs) in whole-exome sequencing data.

Here, as an example, we estimate the selection intensity of SNVs in lung adenocarcinoma (LUAD) data.

Import and preprocess the data
==============================

The data for this pipeline needs to be in [MAF format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification), hg19 coordinates, and have the following headers:

Import the NCI data and the non-NCI data
----------------------------------------

Obtain the NCI data from <https://portal.gdc.cancer.gov/> and store it within the `input_data/` folder. The `YG.data` dataset may be found along with Cannataro *et al.* (2017)[1]

``` r
# import the LUAD data from NCI 
# UUID: 81ccaef3-4550-494d-882c-895fb5a3de3b
NCI.data.38 <- read.csv(file = "input_data/TCGA.LUAD.mutect.81ccaef3-4550-494d-882c-895fb5a3de3b.DR-7.0.somatic.maf",
                        header = T,
                        skip = 5,
                        sep = "\t",
                        stringsAsFactors = F)

# import data from a non-NCI source
YG.data <- read.csv(file = "input_data/adc_inc_counts.txt",
                    header = T,
                    sep = "\t",
                    stringsAsFactors = F)
```

Preprocess the data
-------------------

``` r
# This downloads the file hg38ToHg19.over.chain.gz from 
# http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/ 
# and then loads in a function that does the converting
source("R/hg39_to_hg19_converter.R") 

source("R/unique_tumor_addition.R") # Adds tumor name to dataframe
source("R/flip_function.R") # Finds nucleotide complement
source("R/DNP_remover.R") # removes possible DNP, and recurrent (non-primary) tumors from TCGA data
source("R/tumor_allele_adder.R") # adds column with the tumor allele

NCI.data.19 <- hg38.to.hg19.converter(chain='input_data/hg38Tohg19.chain',hg38_maf=NCI.data.38)
```

    ## Loading in specified MAF...

    ## Warning: closing unused connection 5 (input_data/hg38Tohg19.chain)

    ## Number of rows in the MAF that failed to convert:  11

``` r
NCI.data.19 <- DNP.remover(MAF = NCI.data.19)
```

    ## Removing possible DNP

    ## Total count of potential DNP removed:  6930

    ## DNP removal complete

    ## Deleting any mutations detected in TCGA recurrent tumors

``` r
NCI.data.19$Tumor_Seq_Allele2 <- toupper(NCI.data.19$Tumor_Seq_Allele2)
NCI.data.19$Reference_Allele <- toupper(NCI.data.19$Reference_Allele)
NCI.data.19 <- tumor.allele.adder(MAF = NCI.data.19)
```

[1] Cannataro, V. L., Gaffney, S. G., Stender, C., Zhao, Z., Philips, M., Greenstein, A. E., Townsend, J. P. (2017) “Heterogeneity and mutation in KRAS and associated oncogenes: evaluating the potential for the evolution of resistance to targeting of KRAS G12C" Oncogene, in press
