
# Combining all Yale-Gilead tumors into a single file

CESC.maf <- read.csv(file = "~/Documents/Selection_analysis/CESC/Yale-Gilead/mutationsTN_32_Cervical_cancer.maf",header = T,sep = "\t",stringsAsFactors = F)
head(CESC.maf)
CESC.maf$tumor_type <- "CESC"

LIHC.maf <- read.csv(file = "~/Documents/Selection_analysis/LIHC/Yale-Gilead/mutationsTN_29_Hepatocellular_Carcinoma.maf",header = T,sep = "\t",stringsAsFactors = F)
head(LIHC.maf)
LIHC.maf$tumor_type <- "LIHC"

LUAD.maf <- read.csv(file = "~/Documents/manuscripts/SNV_selection_intensity/SNV_selection_intensity_code/input_data/YG_LUAD_data.txt",sep="\t",header = T,stringsAsFactors = F)
head(LUAD.maf)
LUAD.maf$tumor_type <- "LUAD"
colnames(LUAD.maf)[1] <- "Patient_ID"
colnames(LUAD.maf)[2] <- "Hugo_Symbol"

LUSC.maf <- read.csv(file = "~/Documents/Selection_analysis/LUSC/Yale-Gilead/mutationsTN_23_Lung_Squamous_Cell_Carcinoma.maf",sep="\t",header=T,stringsAsFactors = F)
head(LUSC.maf)
LUSC.maf$tumor_type <- "LUSC"

PAAD.maf <- read.csv(file = "~/Documents/Selection_analysis/PAAD/Yale-Gilead/new_mutationsTN_26_Pancreatic_Cancer.maf",sep = "\t", header = T,stringsAsFactors = F)
head(PAAD.maf)
PAAD.maf$tumor_type <- "PAAD"

SKCMP.maf <- read.csv(file = "~/Documents/Selection_analysis/SKCM/Yale-Gilead/Yale_data_primary.txt",sep = "\t",header = T,stringsAsFactors = F)
head(SKCMP.maf)
SKCMP.maf$tumor_type <- "SKCMP"

SKCMM.maf <- read.csv(file = "~/Documents/Selection_analysis/SKCM/Yale-Gilead/Yale_data_met.txt",sep = "\t",header = T,stringsAsFactors = F)
head(SKCMM.maf)
SKCMM.maf$tumor_type <- "SKCMM"

THCA.maf <- read.csv(file= "~/Documents/Selection_analysis/THCA/Yale-Gilead/mutationsTN_36_Anaplastic_thyroid_carcinoma.maf",sep="\t",stringsAsFactors = F,header = T)
head(THCA.maf)
THCA.maf$tumor_type <- "THCA"

UCEC.maf <- read.csv(file = "~/Documents/Selection_analysis/UCEC/Yale-Gilead/mutationsTN_44_Endometrial_sarcoma.maf",sep="\t",header = T,stringsAsFactors = F)
head(UCEC.maf)
UCEC.maf$tumor_type <- "UCEC"

#First, merge LUAD to get rid of the extra columns. Merge along common headers 
merged.df <- rbind(LUAD.maf[,intersect(colnames(LUAD.maf),colnames(UCEC.maf))],UCEC.maf[,intersect(colnames(LUAD.maf),colnames(UCEC.maf))])
common.heads <- Reduce(intersect,list(names(merged.df),names(CESC.maf),names(LIHC.maf),names(LUSC.maf),names(PAAD.maf),names(SKCMP.maf),names(SKCMM.maf),names(THCA.maf)))

#remove info not needed

common.heads
# common.heads <- common.heads[c(1,2,3,4,5,7,8,10,12,13,14,15,16,17,18)]

merged.df <- rbind(merged.df[,common.heads],
                   CESC.maf[,common.heads],
                   LIHC.maf[,common.heads],
                   LUSC.maf[,common.heads],
                   PAAD.maf[,common.heads],
                   SKCMP.maf[,common.heads],
                   SKCMM.maf[,common.heads],
                   THCA.maf[,common.heads])

head(merged.df)
tail(merged.df)
length(unique(merged.df$tumor_type))
merged.df <- merged.df[which(merged.df$Variant_Type=="SNP"),]

# merged.df[which((merged.df$Patient_ID==merged.df$Tumor_Sample_Barcode)==F),]


merged.df <- merged.df[,c(2,3,5,7,8,10,12,13,15,16,18)]

length(unique(merged.df$tumor_type))




write.table(x = merged.df,file = "output_data/YG_data_merged.txt",quote = F,row.names = F,sep="\t")
