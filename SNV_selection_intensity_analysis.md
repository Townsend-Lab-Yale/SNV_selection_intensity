SNV selection intensity
================
Vincent L. Cannataro

-   [Import and preprocess the data](#import-and-preprocess-the-data)
    -   [Import the NCI data and any non-NCI data](#import-the-nci-data-and-any-non-nci-data)
    -   [Preprocess the data](#preprocess-the-data)
-   [Calculate mutation rate with MutSigCV](#calculate-mutation-rate-with-mutsigcv)
-   [Calculate the trinucleotide mutation profile of this tissue](#calculate-the-trinucleotide-mutation-profile-of-this-tissue)
-   [Calculating selection intensity](#calculating-selection-intensity)

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
# First, need to source functions

# This downloads and unzips the file hg38ToHg19.over.chain.gz from 
# http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/  (if it is not already there)
# and then loads in a function that does the converting
source("R/hg39_to_hg19_converter.R") 

source("R/unique_tumor_addition.R") # Adds tumor name to dataframe
source("R/flip_function.R") # Finds nucleotide complement
source("R/DNP_remover.R") # removes possible DNP, and recurrent (non-primary) tumors from TCGA data
source("R/tumor_allele_adder.R") # adds column with the tumor allele

## Then, use the functions. 
NCI.data.19 <- hg38.to.hg19.converter(chain='input_data/hg38Tohg19.chain',hg38_maf=NCI.data.38)
```

    ## Loading in specified MAF...

    ## Warning in base.OK: closing unused connection 5 (input_data/
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

If necessary, combine the NCI data with any other dataset along common headers, add custom headers again.

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

``` r
MAF_for_analysis <- unique.tumor.addition.function(MAF.file = MAF_for_analysis,non.TCGA.characters.to.keep = 'all',figures=F)
```

    ## Summary statistics of the number of mutations per unique tumor:

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     1.0    94.0   217.0   336.1   454.0  2681.0

``` r
MAF_for_analysis <- DNP.remover(MAF = MAF_for_analysis)
```

    ## Removing possible DNP

    ## Total count of potential DNP removed:  0

    ## DNP removal complete

    ## Deleting any mutations detected in TCGA recurrent tumors

``` r
MAF_for_analysis$Tumor_Seq_Allele2 <- toupper(MAF_for_analysis$Tumor_Seq_Allele2)
MAF_for_analysis$Reference_Allele <- toupper(MAF_for_analysis$Reference_Allele)
MAF_for_analysis <- tumor.allele.adder(MAF = MAF_for_analysis)

save(MAF_for_analysis, file = "output_data/MAF_for_analysis_LUAD.RData")
```

Match the gene names in the `Hugo_Symbol` column in the MAF file and `geneName` column in the isoform file with gene names in the `MutSigCV` covariates file. The isoform file was downloaded from the [genome.ucsc.edu](https://genome.ucsc.edu/cgi-bin/hgTables) table browser, specifying `group: Genes and Gene Predictions`, `track: RefSeq Genes`, and `table: refFlat` in the human hg19 genome. Gene expression data are median values from lung data sets in `The Cancer Cell Line Encyclopedia`[2] matched to the genes and data in the [MutSigCV covariate input file](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/MutSigCV) and saved as `input_data/LUNG_expression.txt`.

``` r
source("R/synonym_matcher.R")

MAF_for_analysis <- synonym.matcher.function(input_to_be_changed = MAF_for_analysis,synonym.table = read.csv(file='input_data/matching_table.txt',stringsAsFactors = F,sep='\t'),covariates.file = read.table(file='input_data/LUNG_expression.txt',sep='\t',header = T,stringsAsFactors = F),isoforms_or_MAF = 'MAF')
```

    ## [1] "More than 1 match in the synonym table for KMT2D in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for DXO in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for PTPN20 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for TENM1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for GALNT15 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for AOC1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for KMT2B in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for DNAAF3 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CFAP47 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for B4GAT1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for HID1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for AZIN2 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CCAR2 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for ZNF106 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CMC4 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for NT5C3A in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for ACKR2 in the covariates file! Could not assign."

    ## Number of changed gene names:  581

``` r
isoforms <- synonym.matcher.function(input_to_be_changed = read.csv(file='input_data/UCSC_RefSeqGenes_refFlat.txt',header = T,sep='\t',stringsAsFactors = F),synonym.table = read.csv(file='input_data/matching_table.txt',stringsAsFactors = F,sep='\t'),covariates.file = read.table(file='input_data/LUNG_expression.txt',sep='\t',header = T,stringsAsFactors = F),isoforms_or_MAF = 'isoforms')
```

    ## [1] "More than 1 match in the synonym table for AZIN2 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for NECTIN4 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for LNPK in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for GALNT15 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for ACKR2 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for NECTIN3 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for MELTF in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for SELENOP in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for DXO in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for NT5C3A in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for STAG3L1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for SEM1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for AOC1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CCAR2 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for PLPP6 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CFAP47 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for TENM1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CMC4 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for WASHC2C in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for PTPN20 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for FAM25BP in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for WASHC2A in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for B4GAT1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for KMT2D in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CFAP73 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for PCNX4 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for ZNF106 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CRAMP1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for HID1 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for PLPPR3 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for MYDGF in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for PLPPR2 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for BORCS8-MEF2B in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for KMT2B in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for CGB3 in the covariates file! Could not assign."
    ## [1] "More than 1 match in the synonym table for DNAAF3 in the covariates file! Could not assign."

    ## Number of changed gene names:  850

``` r
save(MAF_for_analysis,file='output_data/MAF_for_analysis_LUAD.RData')
write.table(MAF_for_analysis,file='output_data/MAF_for_analysis_LUAD.txt',quote = F,sep='\t',row.names = F)

save(isoforms,file='output_data/Isoforms_LUAD.RData')
```

Calculate mutation rate with MutSigCV
=====================================

``` r
dir.create("MutSigCV") #Download and unzip MutSigCV 1.41 here 
```

    ## Warning in dir.create("MutSigCV"): 'MutSigCV' already exists

Steps to calculate mutation rates with MutSigCV:

1.  Agree to the MutSigCV license and download MutSigCV 1.41 from <http://archive.broadinstitute.org/cancer/cga/mutsig_download>.
2.  Download MutSigCV 1.41 reference files from <http://archive.broadinstitute.org/cancer/cga/mutsig_run#reference_files>. Save files to the MutSigCV folder and unzip.
3.  Add contents of input\_data/to\_add\_MutSigCV.txt to line 953 of MutSigCV.m.
4.  Run with `input_data/LUNG_expression.txt` as the covariate file and `output_data/MAF_for_analysis_LUAD.txt` as the mutation file

Outputs `*.gene_rates.txt` and `*.overall_rates.txt`, which we use to calculate mutation rate within genes

Calculate the trinucleotide mutation profile of this tissue
===========================================================

Use the `deconstructSigs` package to calculate the average trinucleotide context among tumors with over 50 SNV (as per [3]).

``` r
source("R/trinucleotide_profile.R")

trinuc.mutation_data <- trinuc.profile.function_withweights(input.MAF = MAF_for_analysis,save.figs=F)
```

    ## Loading required package: deconstructSigs

    ## Removing all recurrent mutations...

    ## Finding the number of mutations per tumor

    ## Number of tumors over specified minimum mutation number of 50: 570

    ## Cleaning input to only contain tumors above the minimum...

    ## Calculating trinucleotide mutation counts...

    ## Warning in mut.to.sigs.input(mut.ref = input.MAF, sample.id = "Unique_patient_identifier", : Check ref bases -- not all match context:
    ##   TCGA-17-Z028:chr1:146046444:T:C

    ## Calculating individual tumor mutational signatures...

    ## No id variables; using all as measure variables

    ## Statistical summary of the proportion of the mutational signature in each tumor sample that is 'unknown'

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0000  0.0836  0.1257  0.1330  0.1751  0.3615

``` r
save(trinuc.mutation_data,file='output_data/trinuc_data_LUAD_all.RData')
trinuc.mutation_data.df <- trinuc.mutation_data$trinuc.mutation_data
save(trinuc.mutation_data.df,file='output_data/trinuc_data_LUAD.RData')

p <- ggplot(data=trinuc.mutation_data$trinuc.mutation_data, aes(Downstream, Upstream)) +
  geom_tile(aes(fill = proportion*100), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue", name="Percent")
p <- p + facet_grid(.~section_labels, labeller = label_parsed) 
p <- p +  geom_text(aes(label = round(proportion, 4)*100),size=3)
# p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.ticks = element_blank(),
                            strip.text=element_text(size=15),
                            axis.title.x = element_text(size=15),
                            axis.title.y = element_text(size=15),
                            axis.text.x = element_text(size=12),
                            axis.text.y=element_text(size=12),plot.title = element_text(hjust = 0.5)) 
  # ggtitle(paste("Trinucleotide profile for ",tumor.name,sep=""))
p
```

![](SNV_selection_intensity_analysis_files/figure-markdown_github/trinuc%20profile-1.png)

``` r
ggsave(paste("figures/",tumor.name,"_trinuc_heatmap_noname.pdf",sep=""),height = 2.5,width = 10)
```

Calculating selection intensity
===============================

``` r
source("R/selection_intensity_calc.R")
source("R/flip_function.R")
source("R/lambda_calculation.R")

load("output_data/trinuc_data_LUAD.RData")
load("output_data/Isoforms_LUAD.RData")
load("output_data/MAF_for_analysis_LUAD.RData")

selection.output <- selection.intensity.calculation.function(genes_for_analysis = c("KRAS","TP53","APC","EGFR"),
                                                                 MAF_for_analysis = MAF_for_analysis,
                                                                 this.substitution = c("this_is_not_a_real_gene",34,"T"),
                                                                 trinuc.mutation_data = trinuc.mutation_data.df,
                                                                 LabReference =  isoforms,
                                                                 translations =  read.csv(file = "input_data/translations.csv",header = T,stringsAsFactors = F),
                                                                 mut_rates = read.csv(file="MutSigCV/MutSigCV_1.41/gene_rates.txt",header = T,stringsAsFactors = F,sep="\t"),
                                                                 low.mut = read.csv(file="MutSigCV/MutSigCV_1.41/overall_rates.txt",header = F,stringsAsFactors = F,sep="\t"),tumor.number = length(unique(MAF_for_analysis$Unique_patient_identifier)),mutsig_siggenes = read.csv(file="MutSigCV/MutSigCV_1.41/mutsig_output/.sig_genes.txt",header = T,stringsAsFactors = F,sep="\t"))


save(selection.output,file = "output_data/selection_intensity_output.RData")
```

[1] Cannataro, V. L., Gaffney, S. G., Stender, C., Zhao, Z., Philips, M., Greenstein, A. E., Townsend, J. P. (2017) “Heterogeneity and mutation in KRAS and associated oncogenes: evaluating the potential for the evolution of resistance to targeting of KRAS G12C" Oncogene, in press

[2] Barretina, Jordi, Giordano Caponigro, Nicolas Stransky, Kavitha Venkatesan, Adam A. Margolin, Sungjoon Kim, Christopher J. Wilson, et al. 2012. “The Cancer Cell Line Encyclopedia Enables Predictive Modelling of Anticancer Drug Sensitivity.” Nature 483 (7391): 603–7.

[3] Rosenthal, Rachel, Nicholas McGranahan, Javier Herrero, Barry S. Taylor, and Charles Swanton. 2016. “DeconstructSigs: Delineating Mutational Processes in Single Tumors Distinguishes DNA Repair Deficiencies and Patterns of Carcinoma Evolution.” Genome Biology 17 (February): 31.
