
# This function merges different MAF files along common headers. 
# It also checks for and alerts users about common variants in headers.

merging_TCGA_and_local_MAFdata_function <- function(NCI_data,
                                                   Local_data){
  
  
  #headers necessary for analyses
  important.headers <- c("Hugo_Symbol",
                         "Chromosome",
                         "Tumor_Seq_Allele2",
                         "Variant_Classification",
                         "Variant_Type",
                         "trv_type",
                         "transcript_error",
                         "Reference_Allele",
                         "Start_Position",
                         "strand",
                         "Tumor_Sample_Barcode",
                         "t_ref_count",
                         "t_alt_count")
  
  message("These are the important headers that need to be contained in both files: ");print(important.headers,sep="")
  
  message("Important headers not in NCI_data: ");print(important.headers[which(!(important.headers %in% colnames(NCI_data)))])
  message("Important headers not in Local_data: ");print(important.headers[which(!(important.headers %in% colnames(Local_data)))])
  
  
  message("Making sure all the essential column headers are the same so they can be properly merged...")
  
  
  ###Can do this quicker with a function! --- 
  rename <- function(oldname, newname, dataframe){
    if(dataframe=="Local_data"){a <- colnames(Local_data)}
    if(dataframe=="NCI_data"){a <- colnames(NCI_data)}
    if(length(which(a==newname))==0){
      print(paste(dataframe,"is missing column name header", newname))
      if(length(which(a==oldname))>0){
        a[which(a==oldname)] <- newname
        print(paste(dataframe,"had column name header", oldname,". That was automatically changed to", newname)) 
      }else{
        print(paste("Could not find", newname, "header in ",dataframe,". You need to manually find the appropriate header and change it to",newname))
      }  
    } 
    a
  } 
  
  colnames(Local_data) <- rename(oldname = "STRAND2",newname = "strand",dataframe = "Local_data")
  colnames(Local_data) <- rename(oldname = "Tumor_ref_cov",newname = "t_ref_count",dataframe = "Local_data")
  colnames(Local_data) <- rename(oldname = "Tumor_nonref_cov",newname = "t_alt_count",dataframe = "Local_data")
  colnames(Local_data) <- rename(oldname = "hugo_symbol",newname = "Hugo_Symbol",dataframe = "Local_data")
  colnames(Local_data) <- rename(oldname = "Chrom",newname = "Chromosome",dataframe = "Local_data")
  colnames(NCI_data) <- rename(oldname = "Transcript_Strand",newname = "strand",dataframe = "NCI_data")
  colnames(NCI_data) <- rename(oldname = "Start_position",newname = "Start_Position",dataframe = "NCI_data")
  
  
  message("Still a problem and need to be fixed manually: ") 
  message("Important headers not in NCI_data: ");print(important.headers[which(!(important.headers %in% colnames(NCI_data)))])
  message("Important headers not in Local_data: ");print(important.headers[which(!(important.headers %in% colnames(Local_data)))])
  
  
  
  ### 
  # Need to make sure that the local files do not contain duplicate tumors from TCGA 
  ###
  
  message("Checking if any local tumor data is a duplicate of a TCGA tumor...")
  message("Printing tumor pair with largest number of matches...")
  
  Local_data$mutation_name <- paste(Local_data$Hugo_Symbol,Local_data$Start_Position,Local_data$Tumor_Seq_Allele2,sep="")
  NCI_data$mutation_name <- paste(NCI_data$Hugo_Symbol,NCI_data$Start_Position,NCI_data$Tumor_Seq_Allele2,sep="")
  
  percent_same <- matrix(data = NA,nrow=length(unique(Local_data$Tumor_Sample_Barcode)),ncol=3)
  
  rownames(percent_same) <- unique(Local_data$Tumor_Sample_Barcode)
  colnames(percent_same) <- c("number_of_same_mutations","percent","TCGA_tumor")
  
  
  for(local_tumor in 1:length(unique(Local_data$Tumor_Sample_Barcode))){
    percent_for_each_tcga <- rep(NA,length=length(unique(NCI_data$Tumor_Sample_Barcode)))
    names(percent_for_each_tcga) <- unique(NCI_data$Tumor_Sample_Barcode)
    for(tcga_tumor in 1:length(unique(NCI_data$Tumor_Sample_Barcode))){
      percent_for_each_tcga[tcga_tumor] <- length(which(Local_data$mutation_name[which(Local_data$Tumor_Sample_Barcode==unique(Local_data$Tumor_Sample_Barcode)[local_tumor])] %in% NCI_data$mutation_name[which(NCI_data$Tumor_Sample_Barcode==unique(NCI_data$Tumor_Sample_Barcode)[tcga_tumor])])) 
      
    }
    
    if(length(which(percent_for_each_tcga>0))>0){
      percent_same[local_tumor,c("number_of_same_mutations","percent","TCGA_tumor")] <- c(percent_for_each_tcga[which.max(percent_for_each_tcga)],percent_for_each_tcga[which.max(percent_for_each_tcga)]/length(which(Local_data$Tumor_Sample_Barcode==unique(Local_data$Tumor_Sample_Barcode)[local_tumor])),names(percent_for_each_tcga)[which.max(percent_for_each_tcga)])
    }
    print(rownames(percent_same)[local_tumor])
    print(percent_same[local_tumor,])
  }
  
  if(length(which(percent_same[,"percent"]>.1))>0){
    warning("Some of the local files share more than 10% of mutations with NCI files,
            see `potential_duplicates.RData` in the working directory for more information.")
    save(percent_same,file = "potential_duplicates.RData")
  }
  
  message("Merging the data frames along their common headers...")
  
  #What are the common headers among the data frames to be merged? 
  common.headers <- colnames(NCI_data)[which(colnames(NCI_data) %in% colnames(Local_data))]
  
  
  NCI_and_Local_data <- as.data.frame(matrix(nrow=(nrow(NCI_data)+nrow(Local_data)),ncol = length(common.headers),data=NA))
  colnames(NCI_and_Local_data) <- common.headers
  head(NCI_and_Local_data)
  
  
  for(i in 1:length(common.headers)){
    NCI_and_Local_data[1:nrow(NCI_data),common.headers[i]] <- NCI_data[,common.headers[i]] 
  }
  
  for(i in 1:length(common.headers)){
    NCI_and_Local_data[(nrow(NCI_data)+1):nrow(NCI_and_Local_data),common.headers[i]] <- Local_data[,common.headers[i]]
  }
  
  
  message("Merging Completed") 
  
  return(NCI_and_Local_data)
  
}


combining_MAFs_function <- function(MAF_directory,sum.stats=T){
  
  first.12 <- function(string_to_12){
    return(paste(unlist(strsplit(string_to_12,split = ""))[1:12],collapse = ""))
  }
  #Names of the files in the tumor list directory
  MAF_files <- dir(MAF_directory,pattern="*.maf*") #only includes files with .maf in the name
  
  #####Cleaning  ---- 
  ###Need to deduplicate and combine all these MAFs 
  
  ##Data to prove uniqueness 
  #Gene name
  #Start_position
  #End_position (we will only be looking at SNV so this should be the same) 
  #Reference Allele (should be the same for the same start_position anyway)
  #Tumor_Seq_Allele2
  #First 16 characters of tumor_sample_barcode 
  
  
  
  ###
  #Sometimes, the MAFs have different column names, so we cannot just combine directly by rows.
  #COADREAD downloaded from firebrowse.org is an example of this-- some have 294 columns, some have 300. 
  #Need to cut down to the most common rows, then combine. 
  #Solution: cut down MAF to the essential columns. 
  ###
  
  
  
  
  
  #initialize the main MAF dataframe with the first row of the first MAF file
  main_MAF.holder <- read.csv(paste(MAF_directory,MAF_files[1],sep=""),skip = 3,header=T,sep="\t",stringsAsFactors = F,nrows = 1)
  # main_MAF <- main_MAF[-1,] #And then get rid of the data (will replace soon)-- just want the column names 
  
  #Finding column numbers for the analysis (loop is faster when data frame is converted into a matrix, and calls are column numbers instead of names)
  #now including this in the loop, in case some headers are in a different order or if they are a different size. 
  # Hugo <- which(colnames(main_MAF.holder)=="Hugo_Symbol")
  # Start <- which(colnames(main_MAF.holder)=="Start_position")
  # Variant <- which(colnames(main_MAF.holder)=="Variant_Type")
  # Tumor_Allele2 <- which(colnames(main_MAF.holder)=="Tumor_Seq_Allele2")
  # Tumor_Barcode <- which(colnames(main_MAF.holder)=="Tumor_Sample_Barcode")
  # chrom <- which(colnames(main_MAF.holder)=="Chromosome")
  
  important.headers <- c("Hugo_Symbol",
                         "Chromosome",
                         "Tumor_Seq_Allele2",
                         "Variant_Classification",
                         "Variant_Type",
                         "trv_type",
                         "transcript_error",
                         "Reference_Allele",
                         "Start_Position",
                         "strand",
                         "Tumor_Sample_Barcode",
                         "t_ref_count",
                         "t_alt_count")
  
  
  #placeholder matrix to store the data
  #Is initialized with a lot of rows, to be filled in, and then chopped down later, because added rows as needed is very time consuming
  # main_MAF <- matrix(nrow=1e6,ncol=ncol(main_MAF.holder),data=NA)
  # colnames(main_MAF) <- colnames(main_MAF.holder)[c(1:42,341,342)]
  main_MAF <- NULL
  
  # start_list <- NULL
  # tumor_list <- NULL
  
  # rows_to_combine <- matrix(data=NA,nrow=1e6,ncol=2)
  # counter <- 1
  # size <- length(start_list)
  #For each MAF, add their rows, while checking that the same mutation is not already present for the same tumor
  # profvis( #This is for checking speed and optimizing the loop
  for(i in 1:length(MAF_files)){
    
    this_MAF <- read.csv(paste(MAF_directory,MAF_files[i],sep=""),skip = 3,header=T,sep="\t",stringsAsFactors = F)
    
    
    #Finding column numbers for the analysis (loop is faster when data frame is converted into a matrix, and calls are column numbers instead of nmaes)
    Hugo <- which(colnames(this_MAF)=="Hugo_Symbol")
    Start <- which(colnames(this_MAF)=="Start_position")
    End <- which(colnames(this_MAF)=="End_position")
    Variant_Classification <- which(colnames(this_MAF)=="Variant_Classification")
    Reference_Allele <- which(colnames(this_MAF)=="Reference_Allele")
    Variant <- which(colnames(this_MAF)=="Variant_Type")
    Tumor_Allele1 <- which(colnames(this_MAF)=="Tumor_Seq_Allele1")
    Tumor_Allele2 <- which(colnames(this_MAF)=="Tumor_Seq_Allele2")
    Tumor_Barcode <- which(colnames(this_MAF)=="Tumor_Sample_Barcode")
    chrom <- which(colnames(this_MAF)=="Chromosome")
    transcript_strand <- which(colnames(this_MAF)=="Transcript_Strand")
    t_ref_count <- which(colnames(this_MAF)=="t_ref_count")
    t_alt_count <- which(colnames(this_MAF)=="t_alt_count")
    Entrez_ID <- which(colnames(this_MAF)=="Entrez_Gene_Id")
    Center <- which(colnames(this_MAF)=="Center")
    NCBI_Build <- which(colnames(this_MAF)=="NCBI_Build")
    cDNA_Change <- which(colnames(this_MAF)=="cDNA_Change")
    Codon_Change <- which(colnames(this_MAF)=="Codon_Change")
    Protein_Change <- which(colnames(this_MAF)=="Protein_Change")
    
    combining_columns <- c(Hugo,Entrez_ID,
                           Center,
                           NCBI_Build,
                           chrom,
                           Tumor_Barcode,
                           Start,
                           End,
                           Variant_Classification,
                           Variant,
                           Reference_Allele,
                           Tumor_Allele1,
                           Tumor_Allele2,
                           transcript_strand,
                           t_ref_count,
                           t_alt_count,
                           cDNA_Change,
                           Codon_Change,
                           Protein_Change)
    column_names <- colnames(this_MAF[,combining_columns])
    
    this_MAF <- as.matrix(this_MAF)
    # to_delete <- NULL
    
    main_MAF <- rbind(main_MAF,this_MAF[,combining_columns])
    # }
    message(paste("Percent of MAFs merged: ",round(i/length(MAF_files),4)*100,"%",sep=""))
    # print(paste("Number of mutations merged: ",length(start_list),sep=""))
    # size <- length(start_list)
  }
  
  main_MAF[,"Tumor_Seq_Allele2"] <- toupper(main_MAF[,"Tumor_Seq_Allele2"])
  main_MAF[,"Reference_Allele"] <- toupper(main_MAF[,"Reference_Allele"])
  message("Determining unique tumors...")
  main_MAF <- cbind(main_MAF,unlist(lapply(main_MAF[,"Tumor_Sample_Barcode"],first.12)))
  
  colnames(main_MAF) <- c(colnames(main_MAF)[-length(colnames(main_MAF))],"Unique_patient_identifier")
  
  unique.patients <- unique(main_MAF[,"Unique_patient_identifier"])
  
  
  message("Deleting duplicated mutations...")
  start.positions <- as.numeric(main_MAF[,"Start_position"])
  for(i in 1:length(unique.patients)){
    to_delete <- NULL
    this_MAF <- which(main_MAF[,"Unique_patient_identifier"]==unique.patients[i])
    for(j in 1:length(this_MAF)){
      
      if(length(which(start.positions[this_MAF]==start.positions[this_MAF[j]]))>1){
        
        if(length(which(start.positions[this_MAF]==start.positions[this_MAF[j]] & 
                        trimws(main_MAF[this_MAF,"Chromosome"])==trimws(main_MAF[this_MAF[j],"Chromosome"]) &
                        trimws(main_MAF[this_MAF,"Tumor_Seq_Allele2"])==main_MAF[this_MAF[j],"Tumor_Seq_Allele2"]))>1){
          
          matches <- which(start.positions[this_MAF]==start.positions[this_MAF[j]] & 
                             trimws(main_MAF[this_MAF,"Chromosome"])==trimws(main_MAF[this_MAF[j],"Chromosome"]) &
                             trimws(main_MAF[this_MAF,"Tumor_Seq_Allele2"])==main_MAF[this_MAF[j],"Tumor_Seq_Allele2"])
          matches <- matches[which(matches>j)]
          
          if(length(matches)>0){
            for(jj in 1:length(matches)){
              if(!(matches[jj] %in% to_delete)){
                to_delete <- c(to_delete,matches[jj])
              }
            }
          }
          
        }
      }
    }
    
    if(length(to_delete)>0){
      # break
      # print(paste(i," tumor"));print(to_delete)
      main_MAF <- main_MAF[-this_MAF[to_delete],]
      start.positions <- start.positions[-this_MAF[to_delete]]
    }
    
    message(paste("Percent of MAFs deduplicated: ",round(i/length(unique.patients),4)*100,"%",sep=""))
    
  }
  
  
  
  
  
  # )
  
  ##Cut the main MAF matrix down to just where the data is stored, make it a data frame, and save it
  
  # nrow(main_MAF)
  message("Cleaning data...") 
  main_MAF_DF <- data.frame(main_MAF,stringsAsFactors = F)
  
  # colnames(main_MAF_DF) <- column_names
  # main_MAF_DF <- droplevels(main_MAF_DF)
  main_MAF_DF$Start_position <- as.numeric(main_MAF_DF$Start_position)
  
  
  if(sum.stats){
    unique.patients <- unique(main_MAF_DF$Unique_patient_identifier)
    
    tumor.mutation.number <- NULL
    for(i in 1:length(unique.patients)){
      tumor.mutation.number[i] <- nrow(main_MAF_DF[which(main_MAF_DF$Unique_patient_identifier==unique.patients[i]),])
      
    }
    
    
    hist(tumor.mutation.number,breaks=100,xlab="Number of mutations per tumor",main="Histogram of the number of mutations in each tumor")
    
    message("Summary statistics of the number of mutations per unique tumor:")
    print(summary(tumor.mutation.number))
    
  }
  return(main_MAF_DF)
  
}

