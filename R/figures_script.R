
# Analyzing data off the cluster

# Figure with all high effect size data ---- 

load("~/Documents/Selection_analysis/combined_selection_output.RData")

supp.table.1 <- combined_all_data[which(combined_all_data$freq>1),]
colnames(supp.table.1)[which(colnames(supp.table.1)=="gamma_epistasis")] <- "selection_intensity"
supp.table.1 <- supp.table.1[,-which(colnames(supp.table.1)=="gamma")]

write.table(x = supp.table.1,file = "output_data/combined_selection_data_SUPP1.txt",quote = F,row.names = F,sep="\t")

source("R/label_nudge.R")
library(ggplot2)
unique(combined_all_data$tumor_type)

recurrent.data <- combined_all_data[which(combined_all_data$freq>1),]


recurrent.data$Name <- NA
for(i in 1:nrow(recurrent.data)){
  recurrent.data$Name[i] <- paste(recurrent.data$Gene[i]," ",ifelse(!is.na(recurrent.data$AA_Ref[i]),paste(recurrent.data$AA_Ref[i],recurrent.data$AA_Pos[i],recurrent.data$AA_Change[i],sep=""),"NCSNV"),sep="")
}


recurrent.data$gamma_epistasis <- log10(recurrent.data$gamma_epistasis)




tumor.num.vec <- rep(NA,length=length(unique(recurrent.data$tumor_type)))


tumor.num.matrix <- as.data.frame(matrix(data = NA,nrow=length(unique(recurrent.data$tumor_type)),ncol=3))
colnames(tumor.num.matrix) <- c("tumor","labs","nums")


for(i in 1:length(unique(recurrent.data$tumor_type))){
  
  load(paste("~/Documents/Selection_analysis/",unique(recurrent.data$tumor_type)[i],"/MAF_",unique(recurrent.data$tumor_type)[i],".RData",sep=""))

  tumors <- unique(MAF_for_analysis$Unique_patient_identifier)
  
  TCGA.length <- length(grep(pattern = "TCGA",x = tumors))
  YG.length <- length(tumors) - TCGA.length

  
  tumor.num.matrix$tumor[i] <- unique(recurrent.data$tumor_type)[i]
  
  this.tumor.number <- TCGA.length + YG.length
  tumor.num.vec[i] <- this.tumor.number
  tumor.num.matrix$nums[i] <- this.tumor.number
  tumor.num.matrix$labs[i] <- ifelse(YG.length==0,paste("italic(n)[TCGA]==",TCGA.length,sep=""),paste("italic(n)[TCGA]==",TCGA.length,"~italic(n)[YG]==",YG.length,sep=""))

  print(paste(i,length(tumors)))
}



recurrent.data$tumor_type <- factor(recurrent.data$tumor_type, levels = unique(recurrent.data$tumor_type))

head(recurrent.data)


recurrent.data$tumor_type <- factor(recurrent.data$tumor_type, levels=unique(recurrent.data$tumor_type)[order(tumor.num.vec)])
recurrent.data$tumor_type <- factor(recurrent.data$tumor_type, levels=rev(levels(recurrent.data$tumor_type))[c(1,2,3,4,20,5,6,7,8,9,10,11,23,12,13,14,15,16,17,18,19,21,22)]) #reorder so that related tumors are next to one another. 
recurrent.data$tumor_type <- factor(recurrent.data$tumor_type, levels=rev(levels(recurrent.data$tumor_type)))


tumor.num.matrix <- tumor.num.matrix[order(tumor.num.matrix$nums,decreasing = T),]
tumor.num.matrix <- tumor.num.matrix[c(1,2,3,4,20,5,6,7,8,9,10,11,23,12,13,14,15,16,17,18,19,21,22),]
tumor.num.matrix <- tumor.num.matrix[rev(1:nrow(tumor.num.matrix)),]

to.add <- recurrent.data[which(recurrent.data$tumor_type==unique(recurrent.data$tumor_type)[1]),]
to.add <- to.add[order(to.add$gamma_epistasis,decreasing = T),][1:50,]

top.hits <- to.add
for(i in 2:length(unique(recurrent.data$tumor_type))){
  to.add <- recurrent.data[which(recurrent.data$tumor_type==unique(recurrent.data$tumor_type)[i]),]
  if(nrow(to.add)<50){
    print(paste("Less than 50!",unique(recurrent.data$tumor_type)[i]))
    to.add <- to.add[order(to.add$gamma_epistasis,decreasing = T),][1:nrow(to.add),]
    top.hits <- rbind(top.hits,to.add)
  }else{
    to.add <- to.add[order(to.add$gamma_epistasis,decreasing = T),][1:50,]
    top.hits <- rbind(top.hits,to.add)
  }
}



data.to.push <- data.push.function(data.to.push = top.hits,x.min = 1,x.max = 7.5,touching.distance = 0.03,pushing.distance = 0.03/25,x.data = 'gamma_epistasis',cat.data = 'tumor_type',max.iter = 1e5)


data.to.push$Gene_freq <- NA
for(i in 1:nrow(data.to.push)){
  data.to.push$Gene_freq[i] <- table(data.to.push$Gene)[data.to.push$Gene[i]]
}


data.to.push$Pval <- as.factor((data.to.push$MutSigCV_p<0.05)*1)
data.to.push$Pval_col <- "P > 0.05"
data.to.push$Pval_col[which(data.to.push$Pval==1)] <- "P <= 0.05"

data.to.push$Qval <- as.factor((data.to.push$MutSigCV_q<0.1)*1)
data.to.push$Qval_col <- "Q > 0.1"
data.to.push$Qval_col[which(data.to.push$Qval==1)] <- "Q <= 0.1"

levels(recurrent.data$tumor_type)[which(levels(recurrent.data$tumor_type)=="HNSC_HPVpos")] <- "HPV^{'+'}~HNSC"
recurrent.data$tumor_type[which(recurrent.data$tumor_type=="HNSC_HPVpos")] <- "HPV^{'+'}~HNSC"

levels(recurrent.data$tumor_type)[which(levels(recurrent.data$tumor_type)=="HNSC_HPVneg")] <-"HPV^{'-'}~HNSC"
recurrent.data$tumor_type[which(recurrent.data$tumor_type=="HNSC_HPVneg")] <- "HPV^{'-'}~HNSC"


# which labels need to be lifted higher than others? 
# recurrent.data$label_lifted <- F
# recurrent.data$label_lifted[which((recurrent.data$tumor_type == "SKCMP" & recurrent.data$Name == "BRAF V600E"))] <- T
data.to.push$label_lifted <- F
data.to.push$label_lifted[which((data.to.push$tumor_type == "SKCMP" & data.to.push$Name == "BRAF V600E") |
                                  (data.to.push$tumor_type == "SKCMM" & data.to.push$Name == "BRAF V600E") |
                                  (data.to.push$tumor_type == "THCA" & data.to.push$Name == "BRAF V600E") | 
                                  (data.to.push$tumor_type == "LGG" & data.to.push$Name == "IDH1 R132H") |
                                  (data.to.push$tumor_type == "PAAD" & data.to.push$Name == "KRAS G12R") | 
                                  (data.to.push$tumor_type == "PAAD" & data.to.push$Name == "KRAS G12V") |
                                  (data.to.push$tumor_type == "PAAD" & data.to.push$Name == "KRAS G12D"))] <- T


lolli <- ggplot(data=subset(recurrent.data,gamma_epistasis>3), aes(x=gamma_epistasis,y=tumor_type))
lolli <- lolli + geom_point(shape="I") 

lolli <- lolli + geom_text(data=subset(data.to.push,gamma_epistasis>3),aes(x=new_x,y=ifelse(label_lifted==T, as.numeric(tumor_type)+0.29,as.numeric(tumor_type)+0.2),label=Name),angle=90,size=1.2,hjust = 0,fontface = 'bold')

# 
# lolli <- lolli + geom_text(data=subset(data.to.push,gamma_epistasis>3),aes(x=new_x,y=ifelse(Name %in% c("IDH1 R132H","BRAF V600E","KRAS G12R","KRAS G12V","KRAS G12D"),as.numeric(tumor_type)+0.28,as.numeric(tumor_type)+0.2),label=Name),angle=90,size=1.2,hjust = 0,fontface = 'bold')

# lolli <- lolli + geom_text(data=subset(data.to.push,gamma_epistasis>3),aes(x=new_x,y=as.numeric(tumor_type)+0.2,label=Name),angle=90,size=1.2,hjust = 0,fontface = 'bold')
lolli <- lolli + geom_text(data=subset(data.to.push,Gene_freq>1 & gamma_epistasis>3),aes(x=new_x,y=ifelse(label_lifted==T, as.numeric(tumor_type)+0.29,as.numeric(tumor_type)+0.2),label=Gene,color=Gene),angle=90,size=1.2,hjust = 0,fontface = 'bold') + scale_colour_discrete(guide = FALSE)

# lolli <- lolli + geom_text(data=subset(data.to.push,Gene_freq>1 & gamma_epistasis>3),aes(x=new_x,y=as.numeric(tumor_type)+0.2,label=Gene,color=Gene),angle=90,size=1.2,hjust = 0,fontface = 'bold') + scale_colour_discrete(guide = FALSE)
# data2 <- data.to.push


lolli <-  lolli + geom_text(data = tumor.num.matrix, aes(x=max(recurrent.data$gamma_epistasis)+0.45,y=(1:23)-0.32,label=labs),parse = T,hjust=1,size=3)


pvals <- list(expression(italic("P") <= 0.05), expression(italic("P") > 0.05))
qvals <- list(expression("Q" <= 0.1), expression("Q" > 0.1))

# ifelse(Name=="IDH1 R132H",as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18)

lolli <- lolli + geom_point(data=subset(data.to.push, gamma_epistasis>3),aes(x=new_x,y=as.numeric(tumor_type)+0.1,size=Prop_tumors_with_specific_mut, fill=Qval_col),shape=21,alpha=0.8)+ scale_fill_manual(labels = qvals,values = c("red","black"),breaks=c("Q <= 0.1","Q > 0.1"))

# lolli <- lolli + geom_point(data=subset(data.to.push, gamma_epistasis>3),aes(x=new_x,y=as.numeric(tumor_type)+0.1,size=Prop_tumors_with_specific_mut, fill=Qval_col),shape=21,alpha=0.8)+ scale_fill_manual(labels = qvals,values = c("red","black"),breaks=c("Q <= 0.1","Q > 0.1"))

lolli <- lolli + geom_segment(data=subset(data.to.push, gamma_epistasis>3),aes(x=gamma_epistasis,y=as.numeric(tumor_type)+0.02,xend=new_x,yend=as.numeric(tumor_type)+0.1),size=0.2,alpha=0.8)

lolli <- lolli + theme_bw() + 
  coord_cartesian(ylim=c(1,23.5),xlim=c(3,max(recurrent.data$gamma_epistasis)+0.2)) + 
  scale_size_continuous(breaks = c(0.01,0.05,0.1,0.2,0.4,0.6)) + 
  # labs(size="Prevalence",color="Gene", fill=expression(paste(italic("P")," value"))) + 
  labs(size="Prevalence",color="Gene", fill=expression(paste("Q"," value"))) + 
  guides(size = guide_legend(reverse = TRUE)) + 
  labs(y="Tumor type",x="Selection intensity") + 
  scale_x_reverse(breaks=3:7, labels=expression(10^3, 10^4, 10^5,10^6,10^7)) +#+ scale_y_discrete(labels=as.expression(unlist(sapply(levels(recurrent.data$tumor_type), as.name)))) + 
  scale_y_discrete(labels=parse(text=levels(recurrent.data$tumor_type))) +
  theme(panel.border = element_blank()) + 
  theme(plot.margin = unit(c(1,1,1,1.1),units="cm")) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 1.5, b = 0, l = 0),size=15),
        axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),axis.title.x=element_text(size=15)) +
  theme(legend.position = c(0.1,0.5))

lolli_plot <- ggplot_gtable(ggplot_build(lolli))
lolli_plot$layout$clip[lolli_plot$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)

g <- arrangeGrob(lolli_plot)

ggsave(filename = "figures/full_selection_data.pdf",plot = g,width = 15,height = 14)

# ggsave(filename = "figures/full_selection_data_large.pdf",plot = g,width = 36,height = 24)



# Figure for manuscript----



recurrent.data2 <- recurrent.data[which(recurrent.data$tumor_type=="COAD" |
                                          recurrent.data$tumor_type=="LUAD" |
                                          recurrent.data$tumor_type=="READ" |
                                          recurrent.data$tumor_type=="HPV^{'+'}~HNSC" |
                                          recurrent.data$tumor_type=="HPV^{'-'}~HNSC" |
                                          recurrent.data$tumor_type=="UCEC" |
                                          recurrent.data$tumor_type=="LUSC" |
                                          recurrent.data$tumor_type=="LGG"),]

data.to.push <- data.push.function(data.to.push = top.hits,x.min = 1,x.max = 7.5,touching.distance = 0.05,pushing.distance = 0.03/25,x.data = 'gamma_epistasis',cat.data = 'tumor_type',max.iter = 1e5)




data.to.push2 <- data.to.push[which(data.to.push$tumor_type=="COAD" |
                                      data.to.push$tumor_type=="LUAD" |
                                      data.to.push$tumor_type=="READ" |
                                      data.to.push$tumor_type=="HNSC_HPVpos" |
                                      data.to.push$tumor_type=="HNSC_HPVneg" |
                                      data.to.push$tumor_type=="UCEC" |
                                      data.to.push$tumor_type=="LUSC" |
                                      data.to.push$tumor_type=="LGG"),]

data.to.push2$Gene_freq <- NA
for(i in 1:nrow(data.to.push2)){
  data.to.push2$Gene_freq[i] <- table(data.to.push2$Gene)[data.to.push2$Gene[i]]
}


tumor.num.matrix2 <- tumor.num.matrix[which(tumor.num.matrix$tumor=="COAD" |
                                              tumor.num.matrix$tumor=="LUAD" |
                                              tumor.num.matrix$tumor=="READ" |
                                              tumor.num.matrix$tumor=="HNSC_HPVpos" |
                                              tumor.num.matrix$tumor=="HNSC_HPVneg" |
                                              tumor.num.matrix$tumor=="UCEC" |
                                              tumor.num.matrix$tumor=="LUSC" |
                                              tumor.num.matrix$tumor=="LGG"),]

recurrent.data2$tumor_type <- factor(recurrent.data2$tumor_type, levels=intersect(levels(recurrent.data2$tumor_type),unique(recurrent.data2$tumor_type)))
data.to.push2$tumor_type <- factor(data.to.push2$tumor_type, levels=intersect(levels(data.to.push2$tumor_type),unique(data.to.push2$tumor_type)))
tumor.num.matrix2[3:4,] <- tumor.num.matrix2[4:3,]

data.to.push2$Qval <- as.factor((data.to.push2$MutSigCV_q<0.1)*1)
data.to.push2$Qval_col <- "Q > 0.1"
data.to.push2$Qval_col[which(data.to.push2$Qval==1)] <- "Q <= 0.1"

lolli <- ggplot(data=subset(recurrent.data2,gamma_epistasis>3), aes(x=gamma_epistasis,y=tumor_type))
lolli <- lolli + geom_point(shape="I") 


lolli <- lolli + geom_text(data=subset(data.to.push2,gamma_epistasis>3),aes(x=new_x,y=ifelse(Name=="IDH1 R132H",as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18),label=Name),angle=90,size=2.8,hjust = 0,fontface = 'bold')
lolli <- lolli + geom_text(data=subset(data.to.push2,Gene_freq>1 & gamma_epistasis>3),aes(x=new_x,y=ifelse(Name=="IDH1 R132H",as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18),label=Gene,color=Gene),angle=90,size=2.8,hjust = 0,fontface = 'bold') + scale_colour_discrete(guide = FALSE)
# data2 <- data.to.push

# tumor_labs <- data.frame(x=0,y=1:length(unique(recurrent.data$tumor_label)),lab=unique(recurrent.data$tumor_label))
# tumor_labs$lab <- factor(tumor_labs$lab, levels= levels(recurrent.data$tumor_label))
lolli <-  lolli + geom_text(data = tumor.num.matrix2, aes(x=max(recurrent.data2$gamma_epistasis)+0.8,y=(1:nrow(tumor.num.matrix2))-0.2,label=labs),parse = T,hjust=1,size=4)

# lolli <- lolli + geom_point(data=data.to.push,aes(x=new_x,y=as.numeric(tumor_type)+0.1,size=Prop_tumors_with_specific_mut, fill=factor(col)),shape=21,alpha=0.8)+ scale_fill_manual(values = c("<1%" = "purple", "1–2%" = "blue", "2–3%" = "lightblue", "3–5%" = "green", "5–10%" = "orange", ">10%" = "red"))
pvals <- list(expression(italic("P") <= 0.05), expression(italic("P") > 0.05))
qvals <- list(expression(italic("Q") <= 0.1), expression(italic("Q") > 0.1))

lolli <- lolli + geom_point(data=subset(data.to.push2, gamma_epistasis>3),aes(x=new_x,y=as.numeric(tumor_type)+0.15,size=Prop_tumors_with_specific_mut, fill=Qval_col),shape=21,alpha=0.8)+ scale_fill_manual(labels = qvals,values = c("red","black"),breaks=c("Q <= 0.1","Q > 0.1"))

lolli <- lolli + geom_segment(data=subset(data.to.push2, gamma_epistasis>3),aes(x=gamma_epistasis,y=as.numeric(tumor_type)+0.02,xend=new_x,yend=as.numeric(tumor_type)+0.1),size=0.2,alpha=0.8)
# lolli

lolli <- lolli + theme_bw() + 
  coord_cartesian(ylim=c(1,nrow(tumor.num.matrix2)+0.5),xlim=c(3,max(recurrent.data$gamma_epistasis)+0.1)) + 
  scale_size_continuous(breaks = c(0.01,0.05,0.1,0.2,0.4,0.6)) + 
  labs(size="Prevalence",color="Gene", fill=expression(paste(italic("Q")," value"))) + 
  guides(size = guide_legend(reverse = TRUE)) + 
  labs(y="Tumor type",x="Selection intensity") + 
  scale_x_reverse(breaks=3:7, labels=expression(10^3, 10^4, 10^5,10^6,10^7)) +
  scale_y_discrete(labels=parse(text=levels(recurrent.data2$tumor_type))) +
  theme(panel.border = element_blank()) + 
  theme(plot.margin = unit(c(1,1,1,1.1),units="cm")) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 1.5, b = 0, l = 0),size=15),
        axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),axis.title.x=element_text(size=15)) +
  theme(legend.position = c(0.15,0.8)) + 
  theme(legend.background = element_rect(fill = 'grey95'))

lolli_plot <- ggplot_gtable(ggplot_build(lolli))
lolli_plot$layout$clip[lolli_plot$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)

g <- arrangeGrob(lolli_plot)
# grid.draw(g)
ggsave(filename = "figures/reduced_selection_data.png",width = 12,height = 12,dpi = 300,plot = g)








