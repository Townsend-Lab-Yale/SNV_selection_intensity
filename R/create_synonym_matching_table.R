# creating synonym matcher tables

# load in Homo_sapiens.gene_info synonym list from here: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/ 
synonyms <- read.csv(file="~/Documents/genome_data/cell_line_gene_expression/Homo_sapiens.gene_info",stringsAsFactors = F,header = T,sep="\t")

# create a list of all synonyms
synonym_list <- list()
total_names <- 0
for(i in 1:nrow(synonyms)){
  synonym_list[[i]] <- unlist(strsplit(synonyms$Synonyms[i],split='\\|'))
  total_names <- total_names + 1 + length(synonym_list[[i]])
}


# create a matrix with all synonyms and their corresponding "Symbol"
synonyms_matrix <- NULL
for(i in 1:nrow(synonyms)){
  to.add <- matrix(nrow=(1+length(synonym_list[[i]])),ncol=3)  
  to.add[,2] <- synonyms$Symbol[i]
  to.add[,3] <- synonyms$chromosome[i]
  to.add[1,1] <- synonyms$Symbol[i]
  to.add[2:nrow(to.add),1] <- synonym_list[[i]]
  synonyms_matrix <- rbind(synonyms_matrix,to.add)
}


# saving the matrix 
write.table(synonyms_matrix,file="~/Documents/genome_data/synonyms/matching_table.txt",quote = F,sep="\t",row.names = F)