
# Analyzing data off the cluster
# 


# Figure with all high effect size data ---- 

load("~/Documents/Selection_analysis/combined_selection_output_withlargeSKCM.RData")

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


tumor.num.matrix <- as.data.frame(matrix(data = NA,nrow=length(unique(recurrent.data$tumor_type)),ncol=4))
colnames(tumor.num.matrix) <- c("tumor","labs","nums","total")


for(i in 1:length(unique(recurrent.data$tumor_type))){
  
  #Load in the MAF files corresponding to the tumor types of interest. 
  load(paste("~/Documents/Selection_analysis/",unique(recurrent.data$tumor_type)[i],"/MAF_",unique(recurrent.data$tumor_type)[i],".RData",sep=""))

  tumors <- unique(MAF_for_analysis$Unique_patient_identifier)
  
  TCGA.length <- length(grep(pattern = "TCGA",x = tumors))
  YG.length <- length(tumors) - TCGA.length

  
  tumor.num.matrix$tumor[i] <- unique(recurrent.data$tumor_type)[i]
  
  this.tumor.number <- TCGA.length + YG.length
  tumor.num.vec[i] <- this.tumor.number
  tumor.num.matrix$nums[i] <- this.tumor.number
  tumor.num.matrix$labs[i] <- ifelse(YG.length==0,paste("italic(n)[TCGA]==",TCGA.length,sep=""),paste("italic(n)[TCGA]==",TCGA.length,"~italic(n)[YG]==",YG.length,sep=""))
  tumor.num.matrix$total[i] <- paste("italic(n)==",this.tumor.number,sep="")

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


# All data, top to bottom, 2e4 and above ---- 



data.to.push <- data.push.function(data.to.push = top.hits,x.min = log10(2e4),x.max = 7.5,touching.distance = 0.03,pushing.distance = 0.03/25,x.data = 'gamma_epistasis',cat.data = 'tumor_type',max.iter = 1e5,gamma.min = log10(2e4))


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


recurrent.data$tumor_type <- factor(recurrent.data$tumor_type, levels=rev(levels(recurrent.data$tumor_type)))
data.to.push$tumor_type <- factor(data.to.push$tumor_type, levels=rev(levels(data.to.push$tumor_type)))





lolli <- ggplot(data=subset(recurrent.data,gamma_epistasis>log10(2e4)), aes(y=gamma_epistasis,x=tumor_type))
lolli <- lolli + geom_point(shape="—") 


lolli <- lolli + geom_text(data=subset(data.to.push,gamma_epistasis>log10(2e4)),aes(y=new_x,x=ifelse(label_lifted==T, as.numeric(tumor_type)+0.24,as.numeric(tumor_type)+0.21),label=Name),angle=0,size=2.8,hjust = 0,fontface = 'bold')
lolli <- lolli + geom_text(data=subset(data.to.push,Gene_freq>1 & gamma_epistasis>log10(2e4)),aes(y=new_x,x=ifelse(label_lifted==T, as.numeric(tumor_type)+0.24,as.numeric(tumor_type)+0.21),label=Gene,color=Gene),angle=0,size=2.8,hjust=0,fontface = 'bold') + scale_colour_discrete(guide = FALSE)



# data2 <- data.to.push

# tumor_labs <- data.frame(x=0,y=1:length(unique(recurrent.data$tumor_label)),lab=unique(recurrent.data$tumor_label))
# tumor_labs$lab <- factor(tumor_labs$lab, levels= levels(recurrent.data$tumor_label))
# tumor.num.matrix[2,"labs"] <- "'italic(n)[TCGA]==104\n'~'italic(n)[YG]==47'"
lolli <-  lolli + geom_text(data = tumor.num.matrix, aes(y=log10(2e4)-0.32,x=(nrow(tumor.num.matrix)):1,label=total),parse = T,hjust=0,size=4)


HPV.status <- as.data.frame(matrix(nrow=2,ncol=1,data=c("HPV^{'+'}","HPV^{'-'}")))
lolli <-  lolli + geom_text(data = HPV.status, aes(y=log10(2e4)-.27,x=c(13,12),label=V1),parse = T,hjust=0,size=4.5)

# lolli <- lolli + geom_point(data=data.to.push,aes(x=new_x,y=as.numeric(tumor_type)+0.1,size=Prop_tumors_with_specific_mut, fill=factor(col)),shape=21,alpha=0.8)+ scale_fill_manual(values = c("<1%" = "purple", "1–2%" = "blue", "2–3%" = "lightblue", "3–5%" = "green", "5–10%" = "orange", ">10%" = "red"))
pvals <- list(expression(italic("P") <= 0.05), expression(italic("P") > 0.05))
qvals <- list(expression(italic("Q") <= 0.1), expression(italic("Q") > 0.1))

lolli <- lolli + geom_point(data=subset(data.to.push, gamma_epistasis>log10(2e4)),aes(y=new_x,x=as.numeric(tumor_type)+0.15,size=Prop_tumors_with_specific_mut, fill=Qval_col),shape=21,alpha=0.8)+ scale_fill_manual(labels = qvals,values = c("red","black"),breaks=c("Q <= 0.1","Q > 0.1"))

lolli <- lolli + geom_segment(data=subset(data.to.push, gamma_epistasis>log10(2e4)),aes(y=gamma_epistasis,x=as.numeric(tumor_type)+0.02,yend=new_x,xend=as.numeric(tumor_type)+0.1),size=0.2,alpha=0.8)
# lolli

xlabels <- levels(recurrent.data$tumor_type)
xlabels[12] <- "HNSC"
xlabels[13] <- "HNSC"
# 
# xlabels[3] <- "'ahhhh\nWe' ~are%=>% bold(here[123])"
# xlabels[3] <- "'HNSC\nHPV' ~^{'+'}"

lolli <- lolli + theme_bw() + 
  coord_cartesian(xlim=c(1,nrow(tumor.num.matrix)+0.5),ylim=c(log10(2e4),max(recurrent.data$gamma_epistasis))) + 
  scale_size_continuous(breaks = c(0.01,0.05,0.1,0.2,0.4,0.6)) + 
  labs(size="Prevalence",color="Gene", fill=expression(paste(italic("Q")," value"))) + 
  guides(size = guide_legend(reverse = TRUE)) + 
  labs(x="Tumor type",y="Selection intensity") + 
  scale_y_continuous(breaks=3:7, labels=expression(10^3, 10^4, 10^5,10^6,10^7)) +
  scale_x_discrete(labels=parse(text=xlabels)) +
  theme(panel.border = element_blank()) + 
  theme(plot.margin = unit(c(t=1,r=1.1,b=1.5,l=1),units="cm")) + 
  theme(axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0),size=15),
        axis.text.x = element_text(size=15,hjust=0,color="black"),axis.text.y = element_text(size=15,color="black"),axis.title.y=element_text(size=15)) +
  theme(legend.position = c(0.15,0.70)) + 
  theme(legend.background = element_rect(fill = 'grey95')) + theme(axis.title.x = element_text(vjust=-12))
# theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

lolli_plot <- ggplot_gtable(ggplot_build(lolli))
lolli_plot$layout$clip[lolli_plot$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)

g <- arrangeGrob(lolli_plot)
# grid.draw(g)
ggsave(filename = "figures/selection_data_2e4_TtoB.png",width = 36,height = 12,dpi = 300,plot = g)
ggsave(filename = "figures/selection_data_2e4_TtoB.pdf",width = 36,height = 12,plot = g,device=cairo_pdf)











# Lollipop figure for manuscript----



recurrent.data2 <- recurrent.data[which(recurrent.data$tumor_type=="SKCMP" |
                                          recurrent.data$tumor_type=="LUAD" |
                                          recurrent.data$tumor_type=="READ" |
                                          recurrent.data$tumor_type=="HPV^{'+'}~HNSC" |
                                          recurrent.data$tumor_type=="HPV^{'-'}~HNSC" |
                                          recurrent.data$tumor_type=="UCEC" |
                                          recurrent.data$tumor_type=="LUSC" |
                                          recurrent.data$tumor_type=="LGG"),]

data.to.push <- data.push.function(data.to.push = top.hits,x.min = 1,x.max = 7.5,touching.distance = 0.05,pushing.distance = 0.03/25,x.data = 'gamma_epistasis',cat.data = 'tumor_type',max.iter = 1e5)




data.to.push2 <- data.to.push[which(data.to.push$tumor_type=="SKCMP" |
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


tumor.num.matrix2 <- tumor.num.matrix[which(tumor.num.matrix$tumor=="SKCMP" |
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


data.to.push2$label_lifted <- F
data.to.push2$label_lifted[which((data.to.push2$tumor_type == "SKCMP" & data.to.push2$Name == "BRAF V600E") | 
                                  (data.to.push2$tumor_type == "LGG" & data.to.push2$Name == "IDH1 R132H"))] <- T


data.to.push2$tumor_type <- factor(data.to.push2$tumor_type, levels=c("READ","SKCMP","HNSC_HPVpos","HNSC_HPVneg","LGG","UCEC","LUSC","LUAD"))

recurrent.data2$tumor_type <- factor(recurrent.data2$tumor_type, levels=c("READ","SKCMP","HPV^{'+'}~HNSC","HPV^{'-'}~HNSC","LGG","UCEC","LUSC","LUAD"))

tumor.num.matrix2 <- tumor.num.matrix2[c(1,5,4,2,3,6,7,8),] #reorder to our custom order



lolli <- ggplot(data=subset(recurrent.data2,gamma_epistasis>3), aes(x=gamma_epistasis,y=tumor_type))
lolli <- lolli + geom_point(shape="I") 


lolli <- lolli + geom_text(data=subset(data.to.push2,gamma_epistasis>3),aes(x=new_x,y=ifelse(label_lifted==T, as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18),label=Name),angle=90,size=2.8,hjust = 0,fontface = 'bold')
lolli <- lolli + geom_text(data=subset(data.to.push2,Gene_freq>1 & gamma_epistasis>3),aes(x=new_x,y=ifelse(label_lifted==T, as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18),label=Gene,color=Gene),angle=90,size=2.8,hjust = 0,fontface = 'bold') + scale_colour_discrete(guide = FALSE)
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
  theme(legend.background = element_rect(fill = 'grey95')) #+ 
  # theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

lolli_plot <- ggplot_gtable(ggplot_build(lolli))
lolli_plot$layout$clip[lolli_plot$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)

g <- arrangeGrob(lolli_plot)
# grid.draw(g)
ggsave(filename = "figures/reduced_selection_data.png",width = 12,height = 12,dpi = 300,plot = g)



# For manuscript, with 2e4 bound ---- 




recurrent.data2 <- recurrent.data[which(recurrent.data$tumor_type=="SKCMP" |
                                          recurrent.data$tumor_type=="LUAD" |
                                          recurrent.data$tumor_type=="READ" |
                                          recurrent.data$tumor_type=="HPV^{'+'}~HNSC" |
                                          recurrent.data$tumor_type=="HPV^{'-'}~HNSC" |
                                          recurrent.data$tumor_type=="UCEC" |
                                          recurrent.data$tumor_type=="LUSC" |
                                          recurrent.data$tumor_type=="LGG"),]

data.to.push <- data.push.function(data.to.push = top.hits,x.min = log10(2e4),x.max = 7.5,touching.distance = 0.04,pushing.distance = 0.03/25,x.data = 'gamma_epistasis',cat.data = 'tumor_type',max.iter = 1e5,gamma.min = log10(2e4))

head(data.to.push[which(data.to.push$tumor_type=="SKCMP"),])


data.to.push2 <- data.to.push[which(data.to.push$tumor_type=="SKCMP" |
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


tumor.num.matrix2 <- tumor.num.matrix[which(tumor.num.matrix$tumor=="SKCMP" |
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


data.to.push2$label_lifted <- F
data.to.push2$label_lifted[which((data.to.push2$tumor_type == "SKCMP" & data.to.push2$Name == "BRAF V600E") | 
                                   (data.to.push2$tumor_type == "LGG" & data.to.push2$Name == "IDH1 R132H"))] <- T


data.to.push2$tumor_type <- factor(data.to.push2$tumor_type, levels=c("SKCMP","HNSC_HPVpos","HNSC_HPVneg","READ","LGG","UCEC","LUSC","LUAD"))

recurrent.data2$tumor_type <- factor(recurrent.data2$tumor_type, levels=c("SKCMP","HPV^{'+'}~HNSC","HPV^{'-'}~HNSC","READ","LGG","UCEC","LUSC","LUAD"))

tumor.num.matrix2 <- tumor.num.matrix2[c(5,4,2,1,3,6,7,8),] #reorder to our custom order



lolli <- ggplot(data=subset(recurrent.data2,gamma_epistasis>log10(2e4)), aes(x=gamma_epistasis,y=tumor_type))
lolli <- lolli + geom_point(shape="I") 


lolli <- lolli + geom_text(data=subset(data.to.push2,gamma_epistasis>log10(2e4)),aes(x=new_x,y=ifelse(label_lifted==T, as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18),label=Name),angle=90,size=2.8,hjust = 0,fontface = 'bold')
lolli <- lolli + geom_text(data=subset(data.to.push2,Gene_freq>1 & gamma_epistasis>log10(2e4)),aes(x=new_x,y=ifelse(label_lifted==T, as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18),label=Gene,color=Gene),angle=90,size=2.8,hjust = 0,fontface = 'bold') + scale_colour_discrete(guide = FALSE)
# data2 <- data.to.push

# tumor_labs <- data.frame(x=0,y=1:length(unique(recurrent.data$tumor_label)),lab=unique(recurrent.data$tumor_label))
# tumor_labs$lab <- factor(tumor_labs$lab, levels= levels(recurrent.data$tumor_label))
lolli <-  lolli + geom_text(data = tumor.num.matrix2, aes(x=max(recurrent.data2$gamma_epistasis)+0.8,y=(1:nrow(tumor.num.matrix2))-0.2,label=labs),parse = T,hjust=1,size=4)

# lolli <- lolli + geom_point(data=data.to.push,aes(x=new_x,y=as.numeric(tumor_type)+0.1,size=Prop_tumors_with_specific_mut, fill=factor(col)),shape=21,alpha=0.8)+ scale_fill_manual(values = c("<1%" = "purple", "1–2%" = "blue", "2–3%" = "lightblue", "3–5%" = "green", "5–10%" = "orange", ">10%" = "red"))
pvals <- list(expression(italic("P") <= 0.05), expression(italic("P") > 0.05))
qvals <- list(expression(italic("Q") <= 0.1), expression(italic("Q") > 0.1))

lolli <- lolli + geom_point(data=subset(data.to.push2, gamma_epistasis>log10(2e4)),aes(x=new_x,y=as.numeric(tumor_type)+0.15,size=Prop_tumors_with_specific_mut, fill=Qval_col),shape=21,alpha=0.8)+ scale_fill_manual(labels = qvals,values = c("red","black"),breaks=c("Q <= 0.1","Q > 0.1"))

lolli <- lolli + geom_segment(data=subset(data.to.push2, gamma_epistasis>log10(2e4)),aes(x=gamma_epistasis,y=as.numeric(tumor_type)+0.02,xend=new_x,yend=as.numeric(tumor_type)+0.1),size=0.2,alpha=0.8)
# lolli

lolli <- lolli + theme_bw() + 
  coord_cartesian(ylim=c(1,nrow(tumor.num.matrix2)+0.5),xlim=c(log10(2e4),max(recurrent.data$gamma_epistasis)+0.1)) + 
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
  theme(legend.background = element_rect(fill = 'grey95')) #+ 
# theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

lolli_plot <- ggplot_gtable(ggplot_build(lolli))
lolli_plot$layout$clip[lolli_plot$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)

g <- arrangeGrob(lolli_plot)
# grid.draw(g)
ggsave(filename = "figures/reduced_selection_data_2e4.png",width = 12,height = 12,dpi = 300,plot = g)




# right to left 

lolli <- ggplot(data=subset(recurrent.data2,gamma_epistasis>log10(2e4)), aes(x=gamma_epistasis,y=tumor_type))
lolli <- lolli + geom_point(shape="I") 


lolli <- lolli + geom_text(data=subset(data.to.push2,gamma_epistasis>log10(2e4)),aes(x=new_x,y=ifelse(label_lifted==T, as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18),label=Name),angle=90,size=2.8,hjust = 0,fontface = 'bold')
lolli <- lolli + geom_text(data=subset(data.to.push2,Gene_freq>1 & gamma_epistasis>log10(2e4)),aes(x=new_x,y=ifelse(label_lifted==T, as.numeric(tumor_type)+0.23,as.numeric(tumor_type)+0.18),label=Gene,color=Gene),angle=90,size=2.8,hjust = 0,fontface = 'bold') + scale_colour_discrete(guide = FALSE)
# data2 <- data.to.push

# tumor_labs <- data.frame(x=0,y=1:length(unique(recurrent.data$tumor_label)),lab=unique(recurrent.data$tumor_label))
# tumor_labs$lab <- factor(tumor_labs$lab, levels= levels(recurrent.data$tumor_label))
lolli <-  lolli + geom_text(data = tumor.num.matrix2, aes(x=log10(2e4)-0.2,y=(1:nrow(tumor.num.matrix2))-0.2,label=labs),parse = T,hjust=1,size=4)

# lolli <- lolli + geom_point(data=data.to.push,aes(x=new_x,y=as.numeric(tumor_type)+0.1,size=Prop_tumors_with_specific_mut, fill=factor(col)),shape=21,alpha=0.8)+ scale_fill_manual(values = c("<1%" = "purple", "1–2%" = "blue", "2–3%" = "lightblue", "3–5%" = "green", "5–10%" = "orange", ">10%" = "red"))
pvals <- list(expression(italic("P") <= 0.05), expression(italic("P") > 0.05))
qvals <- list(expression(italic("Q") <= 0.1), expression(italic("Q") > 0.1))

lolli <- lolli + geom_point(data=subset(data.to.push2, gamma_epistasis>log10(2e4)),aes(x=new_x,y=as.numeric(tumor_type)+0.15,size=Prop_tumors_with_specific_mut, fill=Qval_col),shape=21,alpha=0.8)+ scale_fill_manual(labels = qvals,values = c("red","black"),breaks=c("Q <= 0.1","Q > 0.1"))

lolli <- lolli + geom_segment(data=subset(data.to.push2, gamma_epistasis>log10(2e4)),aes(x=gamma_epistasis,y=as.numeric(tumor_type)+0.02,xend=new_x,yend=as.numeric(tumor_type)+0.1),size=0.2,alpha=0.8)
# lolli

lolli <- lolli + theme_bw() + 
  coord_cartesian(ylim=c(1,nrow(tumor.num.matrix2)+0.5),xlim=c(log10(2e4),max(recurrent.data$gamma_epistasis)+0.1)) + 
  scale_size_continuous(breaks = c(0.01,0.05,0.1,0.2,0.4,0.6)) + 
  labs(size="Prevalence",color="Gene", fill=expression(paste(italic("Q")," value"))) + 
  guides(size = guide_legend(reverse = TRUE)) + 
  labs(y="Tumor type",x="Selection intensity") + 
  scale_x_continuous(breaks=3:7, labels=expression(10^3, 10^4, 10^5,10^6,10^7)) +
  scale_y_discrete(labels=parse(text=levels(recurrent.data2$tumor_type))) +
  theme(panel.border = element_blank()) + 
  theme(plot.margin = unit(c(1,1,1,1.1),units="cm")) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 1.5, b = 0, l = 0),size=15),
        axis.text.y = element_text(size=15),axis.text.x = element_text(size=15),axis.title.x=element_text(size=15)) +
  theme(legend.position = c(0.85,0.8)) + 
  theme(legend.background = element_rect(fill = 'grey95')) #+ 
# theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

lolli_plot <- ggplot_gtable(ggplot_build(lolli))
lolli_plot$layout$clip[lolli_plot$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)

g <- arrangeGrob(lolli_plot)
# grid.draw(g)
ggsave(filename = "figures/reduced_selection_data_2e4_RtL.png",width = 12,height = 12,dpi = 300,plot = g)


# Top to bottom

# need to change ordering of the plot
recurrent.data2$tumor_type <- factor(recurrent.data2$tumor_type, levels=rev(levels(recurrent.data2$tumor_type)))
data.to.push2$tumor_type <- factor(data.to.push2$tumor_type, levels=rev(levels(data.to.push2$tumor_type)))

lolli <- ggplot(data=subset(recurrent.data2,gamma_epistasis>log10(2e4)), aes(y=gamma_epistasis,x=tumor_type))
lolli <- lolli + geom_point(shape="—") 


lolli <- lolli + geom_text(data=subset(data.to.push2,gamma_epistasis>log10(2e4)),aes(y=new_x,x=ifelse(label_lifted==T, as.numeric(tumor_type)+0.24,as.numeric(tumor_type)+0.21),label=Name),angle=0,size=2.8,hjust = 0,fontface = 'bold')
lolli <- lolli + geom_text(data=subset(data.to.push2,Gene_freq>1 & gamma_epistasis>log10(2e4)),aes(y=new_x,x=ifelse(label_lifted==T, as.numeric(tumor_type)+0.24,as.numeric(tumor_type)+0.21),label=Gene,color=Gene),angle=0,size=2.8,hjust=0,fontface = 'bold') + scale_colour_discrete(guide = FALSE)



# data2 <- data.to.push

# tumor_labs <- data.frame(x=0,y=1:length(unique(recurrent.data$tumor_label)),lab=unique(recurrent.data$tumor_label))
# tumor_labs$lab <- factor(tumor_labs$lab, levels= levels(recurrent.data$tumor_label))
# tumor.num.matrix2[2,"labs"] <- "'italic(n)[TCGA]==104\n'~'italic(n)[YG]==47'"
lolli <-  lolli + geom_text(data = tumor.num.matrix2, aes(y=log10(2e4)-0.3,x=(nrow(tumor.num.matrix2)):1,label=total),parse = T,hjust=0,size=4)


HPV.status <- as.data.frame(matrix(nrow=2,ncol=1,data=c("HPV^{'+'}","HPV^{'-'}")))
lolli <-  lolli + geom_text(data = HPV.status, aes(y=log10(2e4)-.25,x=c(7,6),label=V1),parse = T,hjust=0,size=4.5)

# lolli <- lolli + geom_point(data=data.to.push,aes(x=new_x,y=as.numeric(tumor_type)+0.1,size=Prop_tumors_with_specific_mut, fill=factor(col)),shape=21,alpha=0.8)+ scale_fill_manual(values = c("<1%" = "purple", "1–2%" = "blue", "2–3%" = "lightblue", "3–5%" = "green", "5–10%" = "orange", ">10%" = "red"))
pvals <- list(expression(italic("P") <= 0.05), expression(italic("P") > 0.05))
qvals <- list(expression(italic("Q") <= 0.1), expression(italic("Q") > 0.1))

lolli <- lolli + geom_point(data=subset(data.to.push2, gamma_epistasis>log10(2e4)),aes(y=new_x,x=as.numeric(tumor_type)+0.15,size=Prop_tumors_with_specific_mut, fill=Qval_col),shape=21,alpha=0.8)+ scale_fill_manual(labels = qvals,values = c("red","black"),breaks=c("Q <= 0.1","Q > 0.1"))

lolli <- lolli + geom_segment(data=subset(data.to.push2, gamma_epistasis>log10(2e4)),aes(y=gamma_epistasis,x=as.numeric(tumor_type)+0.02,yend=new_x,xend=as.numeric(tumor_type)+0.1),size=0.2,alpha=0.8)
# lolli

xlabels <- levels(recurrent.data2$tumor_type)
xlabels[6] <- "HNSC"
xlabels[7] <- "HNSC"
# 
# xlabels[3] <- "'ahhhh\nWe' ~are%=>% bold(here[123])"
# xlabels[3] <- "'HNSC\nHPV' ~^{'+'}"

lolli <- lolli + theme_bw() + 
  coord_cartesian(xlim=c(1,nrow(tumor.num.matrix2)+0.5),ylim=c(log10(2e4),max(recurrent.data2$gamma_epistasis))) + 
  scale_size_continuous(breaks = c(0.01,0.05,0.1,0.2,0.4,0.6)) + 
  labs(size="Prevalence",color="Gene", fill=expression(paste(italic("Q")," value"))) + 
  guides(size = guide_legend(reverse = TRUE)) + 
  labs(x="Tumor type",y="Selection intensity") + 
  scale_y_continuous(breaks=3:7, labels=expression(10^3, 10^4, 10^5,10^6,10^7)) +
  scale_x_discrete(labels=parse(text=xlabels)) +
  theme(panel.border = element_blank()) + 
  theme(plot.margin = unit(c(t=1,r=1.1,b=1.5,l=1),units="cm")) + 
  theme(axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0),size=15),
        axis.text.x = element_text(size=15,hjust=0,color="black"),axis.text.y = element_text(size=15,color="black"),axis.title.y=element_text(size=15)) +
  theme(legend.position = c(0.15,0.70)) + 
  theme(legend.background = element_rect(fill = 'grey95')) + theme(axis.title.x = element_text(vjust=-12))
# theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

lolli_plot <- ggplot_gtable(ggplot_build(lolli))
lolli_plot$layout$clip[lolli_plot$layout$name == "panel"] <- "off"
library(grid)
library(gridExtra)

g <- arrangeGrob(lolli_plot)
# grid.draw(g)
ggsave(filename = "figures/reduced_selection_data_2e4_TtoB.png",width = 12,height = 12,dpi = 300,plot = g)



# Tornado plot figure for manuscript ---- 
# 
# subset the data to plot
LUAD.data <- subset(combined_all_data, tumor_type == "LUAD")
LUSC.data <- subset(combined_all_data, tumor_type == "LUSC")

#need to figure out common genes among tumor types and asign colors 
# LUAD.and.LUSC.data <- subset(combined_all_data, (tumor_type == "LUSC" | tumor_type == "LUAD") & freq>1 )
# table(LUAD.and.LUSC.data$Gene)



LUAD.data <- LUAD.data[which(LUAD.data$freq > 1),]

LUAD.data <- LUAD.data[which(!is.na(LUAD.data$Gene)),]
LUAD.data <- LUAD.data[order(-LUAD.data$gamma_epistasis),]

if(nrow(LUAD.data)>25){
  LUAD.data.ordered <- LUAD.data[1:25,] #
}else{
  LUAD.data.ordered <- LUAD.data 
}
LUAD.data.ordered <- LUAD.data.ordered[order(LUAD.data.ordered$gamma_epistasis),]
LUAD.data.ordered$Name <- NA

for(i in 1:nrow(LUAD.data.ordered)){
  LUAD.data.ordered$Name[i] <- paste(LUAD.data.ordered$Gene[i]," ",ifelse(!is.na(LUAD.data.ordered$AA_Ref[i]),paste(LUAD.data.ordered$AA_Ref[i],LUAD.data.ordered$AA_Pos[i],LUAD.data.ordered$AA_Change[i],sep=""),"NCSNV"),sep="")
}

#If the name is not unique, need to make it unique
if(length(which(table(LUAD.data.ordered$Name)>1))>0){
  nonunique <- which(table(LUAD.data.ordered$Name)>1)
  for(k in 1:length(nonunique)){
    these.pos <- which(LUAD.data.ordered$Name==names(nonunique[k]))
    
    for(j in 1:length(these.pos)){
      LUAD.data.ordered$Name[these.pos[j]] <- paste(LUAD.data.ordered$Name[these.pos[j]],j)  
    }
    
  }
  
}


# LUAD.data.ordered$Name <- paste(LUAD.data.ordered$Gene," ",LUAD.data.ordered$AA_Ref,LUAD.data.ordered$AA_Pos,LUAD.data.ordered$AA_Change,sep="")
LUAD.data.ordered$Name <- factor(LUAD.data.ordered$Name, levels=unique(LUAD.data.ordered$Name))

LUAD.data.ordered$Name_col <- "black"
LUAD.data.ordered$Name_col[which(LUAD.data.ordered$MutSigCV_q < 0.1)] <- "red"




LUSC.data <- LUSC.data[which(LUSC.data$freq > 1),]

LUSC.data <- LUSC.data[which(!is.na(LUSC.data$Gene)),]
LUSC.data <- LUSC.data[order(-LUSC.data$gamma_epistasis),]

if(nrow(LUSC.data)>25){
  LUSC.data.ordered <- LUSC.data[1:25,] #
}else{
  LUSC.data.ordered <- LUSC.data 
}
LUSC.data.ordered <- LUSC.data.ordered[order(LUSC.data.ordered$gamma_epistasis),]
LUSC.data.ordered$Name <- NA

for(i in 1:nrow(LUSC.data.ordered)){
  LUSC.data.ordered$Name[i] <- paste(LUSC.data.ordered$Gene[i]," ",ifelse(!is.na(LUSC.data.ordered$AA_Ref[i]),paste(LUSC.data.ordered$AA_Ref[i],LUSC.data.ordered$AA_Pos[i],LUSC.data.ordered$AA_Change[i],sep=""),"NCSNV"),sep="")
}

#If the name is not unique, need to make it unique
if(length(which(table(LUSC.data.ordered$Name)>1))>0){
  nonunique <- which(table(LUSC.data.ordered$Name)>1)
  for(k in 1:length(nonunique)){
    these.pos <- which(LUSC.data.ordered$Name==names(nonunique[k]))
    
    for(j in 1:length(these.pos)){
      LUSC.data.ordered$Name[these.pos[j]] <- paste(LUSC.data.ordered$Name[these.pos[j]],j)  
    }
    
  }
  
}


# LUSC.data.ordered$Name <- paste(LUSC.data.ordered$Gene," ",LUSC.data.ordered$AA_Ref,LUSC.data.ordered$AA_Pos,LUSC.data.ordered$AA_Change,sep="")
LUSC.data.ordered$Name <- factor(LUSC.data.ordered$Name, levels=unique(LUSC.data.ordered$Name))

LUSC.data.ordered$Name_col <- "black"
LUSC.data.ordered$Name_col[which(LUSC.data.ordered$MutSigCV_q < 0.1)] <- "red"



# determining colors for the plot 

which(table(c(LUAD.data.ordered$Gene,LUSC.data.ordered$Gene))>1)
# 5 gene names 
# BRAF CTNNB1   KRAS NFE2L2   TP53 
# 3      6     13     19     28 

LUAD.and.LUSC.ord.bound <- rbind(LUAD.data.ordered,LUSC.data.ordered)

LUAD.and.LUSC.ord.bound$bar_color <- "black"

library(RColorBrewer)

colors.to.give <- brewer.pal(n = 5, name = "Set1")

LUAD.and.LUSC.ord.bound$bar_color[which(LUAD.and.LUSC.ord.bound$Gene=="BRAF")] <- colors.to.give[1]
LUAD.and.LUSC.ord.bound$bar_color[which(LUAD.and.LUSC.ord.bound$Gene=="CTNNB1")] <- colors.to.give[2]
LUAD.and.LUSC.ord.bound$bar_color[which(LUAD.and.LUSC.ord.bound$Gene=="KRAS")] <- colors.to.give[3]
LUAD.and.LUSC.ord.bound$bar_color[which(LUAD.and.LUSC.ord.bound$Gene=="NFE2L2")] <- colors.to.give[4]
LUAD.and.LUSC.ord.bound$bar_color[which(LUAD.and.LUSC.ord.bound$Gene=="TP53")] <- colors.to.give[5]

LUAD.data.ordered <- LUAD.and.LUSC.ord.bound[which(LUAD.and.LUSC.ord.bound$tumor_type=="LUAD"),]
LUSC.data.ordered <- LUAD.and.LUSC.ord.bound[which(LUAD.and.LUSC.ord.bound$tumor_type=="LUSC"),]

#adapted from http://stackoverflow.com/questions/18265941/two-horizontal-bar-charts-with-shared-axis-in-ggplot2-similar-to-population-pyr
library('grid')
library('gridExtra')
library(ggplot2)
source("R/fancy_scientific_code.R")


# LUAD tornado ----

tumor.name <- "LUAD"

g.mid <- ggplot(LUAD.data.ordered,aes(x=1,y=Name)) +
  geom_text(aes(label=Name,color=factor(Name_col)),size=6) +
  # geom_segment(aes(x=0.94,xend=0.96,yend=Name)) +
  # geom_segment(aes(x=1.04,xend=1.065,yend=Name)) +
  ggtitle("") +
  ylab(NULL) + 
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065)) + ggtitle(tumor.name) + scale_colour_manual(values=c("black","red")) +
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA,size=12),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"),legend.position = "none",plot.title = element_text(hjust = 0.5))
g.mid

luad.col <- as.character(LUAD.data.ordered$bar_color)
names(luad.col) <- as.character(LUAD.data.ordered$Gene)

g1 <- ggplot(data=LUAD.data.ordered,aes(x=Name,y=mu,fill=Gene)) +
  geom_bar(stat="identity") + scale_fill_manual(values=luad.col) + ggtitle("Mutation rate")  +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,10), "mm")) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  geom_text(aes(label=round(mu*1e6,2)), vjust=-0.5, hjust=0.5, position=position_dodge(width=0.9),size=3,angle=90) +
  geom_text(aes(label=freq,y=-max(mu)/26), position=position_dodge(width=0.9),size=5,angle=0,color="black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none',axis.text.x = element_text(size=12)) +
  scale_y_reverse(labels=fancy_scientific) + coord_flip() 
g1

g2 <- ggplot(data=LUAD.data.ordered, aes(x=Name,y=gamma_epistasis,fill=Gene)) +
  geom_bar(stat="identity") + ggtitle("Selection intensity") + scale_fill_manual(values=luad.col) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,5,1,-1), "mm")) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=12)) +
  coord_flip() + theme(legend.position = c(0.85, .5),legend.text = element_text(size=10)) + theme(legend.position="none")
g2

gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4.25/10,2/10,3.5/10))


gg.combined.luad <- arrangeGrob(gg1,gg.mid,gg2,ncol=3,widths=c(5/10,3/10,5/10))
ggsave(gg.combined.luad, filename = paste("figures/selection_tornado_plot_",tumor.name,".png",sep=""),units = "in",height=7,width = 10)







# LUSC tornado ----


#from http://stackoverflow.com/questions/18265941/two-horizontal-bar-charts-with-shared-axis-in-ggplot2-similar-to-population-pyr
library('grid')
library('gridExtra')
library(ggplot2)
source("R/fancy_scientific_code.R")

tumor.name <- "LUSC"

g.mid <- ggplot(LUSC.data.ordered,aes(x=1,y=Name)) +
  geom_text(aes(label=Name,color=factor(Name_col)),size=6) +
  # geom_segment(aes(x=0.94,xend=0.96,yend=Name)) +
  # geom_segment(aes(x=1.04,xend=1.065,yend=Name)) +
  ggtitle("") +
  ylab(NULL) + 
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065)) + ggtitle(tumor.name) + scale_colour_manual(values=c("black","red")) +
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA,size=12),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"),legend.position = "none",plot.title = element_text(hjust = 0.5))
g.mid

LUSC.col <- as.character(LUSC.data.ordered$bar_color)
names(LUSC.col) <- as.character(LUSC.data.ordered$Gene)

g1 <- ggplot(data=LUSC.data.ordered,aes(x=Name,y=mu,fill=Gene)) +
  geom_bar(stat="identity") + scale_fill_manual(values=LUSC.col) + ggtitle("Mutation rate")  +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,10), "mm")) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  geom_text(aes(label=round(mu*1e6,2)), vjust=-0.5, hjust=0.5, position=position_dodge(width=0.9),size=3,angle=90) +
  geom_text(aes(label=freq,y=-max(mu)/26), position=position_dodge(width=0.9),size=5,angle=0,color="black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = 'none',axis.text.x = element_text(size=12)) +
  scale_y_reverse(labels=fancy_scientific) + coord_flip() 
g1

g2 <- ggplot(data=LUSC.data.ordered, aes(x=Name,y=gamma_epistasis,fill=Gene)) +
  geom_bar(stat="identity") + ggtitle("Selection intensity") + scale_fill_manual(values=LUSC.col) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,5,1,-1), "mm")) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(size=12)) +
  coord_flip() + theme(legend.position = c(0.85, .5),legend.text = element_text(size=10)) + theme(legend.position="none")
g2

gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4.25/10,2/10,3.5/10))


gg.combined.LUSC <- arrangeGrob(gg1,gg.mid,gg2,ncol=3,widths=c(5/10,3/10,5/10))
ggsave(gg.combined.LUSC, filename = paste("figures/selection_tornado_plot_",tumor.name,".png",sep=""),units = "in",height=7,width = 10)


# Loading in whole data set, then looking at trinucleotide mutations responsible for specific TP53 mutations. 

# load("~/Documents/Selection_analysis/combined_selection_output_full_data.RData")

head(combined_full_data)
combined_full_data[which(combined_full_data$tumor_type=="LUSC" & combined_full_data$Gene=="TP53" & combined_full_data$Amino_acid_position==298),]

combined_full_data[which(combined_full_data$tumor_type=="LUSC" & combined_full_data$Gene=="TP53" & combined_full_data$Amino_acid_position==158),]

combined_full_data[which(combined_full_data$tumor_type=="LUSC" & combined_full_data$Gene=="TP53" & combined_full_data$Amino_acid_position==245 & combined_full_data$Amino_acid_alternative=="C"),]

