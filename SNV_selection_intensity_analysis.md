SNV selection intensity
================
Vincent L. Cannataro

This [R Markdown](http://rmarkdown.rstudio.com) Notebook, output as a `github_document`, describes the process by which we estimate the selection intensity for single nucleotide variants (SNVs) in whole-exome sequencing data.

Here, as an example, we estimate the selection intensity of SNVs in lung adenocarcinoma (LUAD) data.

Import and preprocess the data
==============================

The data for this pipeline needs to be in MAF format

Import NCI data and non-NCI data
--------------------------------

Obtain the NCI data from <https://portal.gdc.cancer.gov/> and store it within the `input_data` folder.

``` r
# import the LUAD data from NCI 
# UUID: 81ccaef3-4550-494d-882c-895fb5a3de3b
NCI.data.38 <- read.csv(file = "input_data/TCGA.LUAD.mutect.81ccaef3-4550-494d-882c-895fb5a3de3b.DR-7.0.somatic.maf",
                        header = T,
                        skip = 5,
                        sep = "\t",
                        stringsAsFactors = F)

# import data from a non-NCI source
YG.data <- read.csv(file = "input_data/adc_inc_counts.txt",header = T,sep = "\t",stringsAsFactors = F)
```
