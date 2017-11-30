SNV selection intensity
================
Vincent L. Cannataro

This document describes our pipeline to estimate the selection intensity for single nucleotide variants (SNVs) in whole-exome sequencing data.

Here, as an example, we estimate the selection intensity of SNVs in lung adenocarcinoma (LUAD) data.

Import and preprocess the data
==============================

The data for this pipeline needs to be in [MAF format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification), hg19 coordinates, and have the following headers:

-   Hugo\_Symbol
-   Chromosome
-   Start\_Position
-   Variant\_Type
-   Tumor\_Sample\_Barcode

Import the NCI data and any non-NCI data
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
YG.data <- read.csv(file = "input_data/YG_LUAD_data.txt",
                    header = T,
                    sep = "\t",
                    stringsAsFactors = F)
```

Preprocess the data
-------------------

Convert the NCI dataset to hg19 coordinates so that the coordinates work with `MutSigCV`. Then, add the unique tumor names based off of `Tumor_Sample_Barcode`, remove possible DNP mistakenly reported as SNP, and add a column to the dataframe that has the "tumor" allele.

``` r
## First, need to source functions

# This downloads and unzips the file hg38ToHg19.over.chain.gz from 
# http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/  (if it is not already there)
# and then loads in a function that does the converting
source("R/hg39_to_hg19_converter.R") 

source("R/unique_tumor_addition.R") # Adds tumor name to dataframe
source("R/flip_function.R") # Finds nucleotide complement
source("R/DNP_remover.R") # removes possible DNP, and recurrent (non-primary) tumors from TCGA data
source("R/tumor_allele_adder.R") # adds column with the tumor allele

## Then, use functions. 
NCI.data.19 <- hg38.to.hg19.converter(chain='input_data/hg38Tohg19.chain',hg38_maf=NCI.data.38)
```

    ## Loading in specified MAF...

    ## Warning in length(x): closing unused connection 5 (input_data/
    ## hg38Tohg19.chain)

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

#Save the NCI dataframe 
save(NCI.data.19, file = "output_data/NCI_data_LUAD.RData")
```

Combine the NCI data with any other dataset along common headers.

``` r
source("R/merging_NCI_and_local_MAF_files.R")

MAF_for_analysis <- merging_TCGA_and_local_MAFdata_function(NCI_data = NCI.data.19,
                                                            Local_data = YG.data)
```

    ## These are the important headers that need to be contained in both files:

    ##  [1] "Hugo_Symbol"            "Chromosome"            
    ##  [3] "Tumor_Seq_Allele2"      "Variant_Classification"
    ##  [5] "Variant_Type"           "trv_type"              
    ##  [7] "transcript_error"       "Reference_Allele"      
    ##  [9] "Start_Position"         "strand"                
    ## [11] "Tumor_Sample_Barcode"   "t_ref_count"           
    ## [13] "t_alt_count"

    ## Important headers not in NCI_data:

    ## [1] "trv_type"         "transcript_error"

    ## Important headers not in Local_data:

    ## [1] "Hugo_Symbol"      "trv_type"         "transcript_error"
    ## [4] "strand"           "t_ref_count"      "t_alt_count"

    ## Making sure all the essential column headers are the same so they can be properly merged...

    ## [1] "Local_data is missing column name header strand"
    ## [1] "Local_data had column name header STRAND2 . That was automatically changed to strand"
    ## [1] "Local_data is missing column name header t_ref_count"
    ## [1] "Local_data had column name header Tumor_ref_cov . That was automatically changed to t_ref_count"
    ## [1] "Local_data is missing column name header t_alt_count"
    ## [1] "Local_data had column name header Tumor_nonref_cov . That was automatically changed to t_alt_count"
    ## [1] "Local_data is missing column name header Hugo_Symbol"
    ## [1] "Local_data had column name header hugo_symbol . That was automatically changed to Hugo_Symbol"

    ## Still a problem and need to be fixed manually:

    ## Important headers not in NCI_data:

    ## [1] "trv_type"         "transcript_error"

    ## Important headers not in Local_data:

    ## [1] "trv_type"         "transcript_error"

    ## Merging the data frames along their common headers...

    ## Merging Completed

[1] Cannataro, V. L., Gaffney, S. G., Stender, C., Zhao, Z., Philips, M., Greenstein, A. E., Townsend, J. P. (2017) “Heterogeneity and mutation in KRAS and associated oncogenes: evaluating the potential for the evolution of resistance to targeting of KRAS G12C" Oncogene, in press
