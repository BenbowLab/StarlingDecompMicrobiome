# Package load
library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)
library(ggpubr)





#####
#Functions
#####


table_for_plotting<-function(Physeq, Level, Variable1,Variable2=NULL){
  #####
  # Inputs
  #####
  # Required
  #     Physeq: A phyloseq object of abundance combined at the desired level (using tax_glom)
  #     Level: The taxonomic level at which you want the rel. abu. summarized 
  #         This is the same level used in tax_glom()
  #     Variable1: The first variable to summarize by 
  # Optional
  #     Variable2: A second variable to summarize by
  # Example syntax: table_for_plotting(PhylumLevel,"Phylum","Region","sample_Species")
  # Note: If you are using 'Species' as a metadata variable, phyloseq will
  # rename it sample_Species to avoid conflicts
  
  
  #####
  # Required packages
  #####
  library(plyr)
  library(phyloseq)
  
  df <- psmelt(Physeq) # Convert phyloseq object to a dataframe
  df$Abundance=df$Abundance*100 # Multiply by 100 so scale is 0-100%
  # The section below will run if only one variable to summarize by is present
  if(is.null(Variable2)==T) {
    # Summarize the abundance based on desired level (phylum,family,etc) and variable of interest (variable 1)
    Trtdata <- ddply(df, c(Level, Variable1), summarise,  
                     N    = length(Abundance),
                     mean = mean(Abundance),
                     sd   = sd(Abundance),
                     se   = sd / sqrt(N))
    # Create a new table with the sum relative of all taxa not filtered out
    Summarized<-ddply(Trtdata,c(Variable1),summarize,sum= sum(mean),N=mean(N)) 
    # Subtract sum of all other taxa to get the 'Other' category rel abu 
    Summarized$mean<-100-Summarized$sum
    
    # Modified the summarized object to match the Trtdata format
    TableForCombining<-data.frame(Variable1=Summarized[,1],N=Summarized$N,Level="Other",mean=Summarized$mean,sd="NA",se="NA")
    # Rename columns based on variable name and level
    colnames(TableForCombining)[1]<-Variable1
    colnames(TableForCombining)[3]<-Level
    # Combine the other category with the already summarized data
    
    TableWithOther<-rbind(Trtdata,TableForCombining)
    TableWithOther[Level] = factor(TableWithOther[,Level], levels = unique(TableWithOther[,Level])) # Moves other to the bottom of taxa list
  }
  
  # The below section will run if two variables to summarize by are present
  else {Trtdata <- ddply(df, c(Level, Variable1,Variable2), summarise,  
                         N    = length(Abundance),
                         mean = mean(Abundance),
                         sd   = sd(Abundance),
                         se   = sd / sqrt(N))
  
  Summarized<-ddply(Trtdata,c(Variable1,Variable2),summarize,sum= sum(mean),N=mean(N))
  Summarized$mean<-100-Summarized$sum
  
  
  TableForCombining<-data.frame(Variable1=Summarized[,1],Variable2=Summarized[,2],N=Summarized$N,Level="Other",mean=Summarized$mean,sd="NA",se="NA")
  colnames(TableForCombining)[1]<-Variable1
  colnames(TableForCombining)[2]<-Variable2
  colnames(TableForCombining)[4]<-Level
  TableWithOther<-rbind(Trtdata,TableForCombining)}
  TableWithOther[Level] = factor(TableWithOther[,Level], levels = unique(TableWithOther[,Level])) # Moves other to the bottom of taxa list
  
  return(TableWithOther)
}
#####
#Import
#####

setwd("~/Documents/Microbiome/2022/Starling/starlingdecomp")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
biom=import_biom("table-with-taxonomy500.biom",parseFunction= parse_taxonomy_greengenes)
meta=read.csv("StarlingMeta5.13.22.csv",header=TRUE)
meta=sample_data(meta)
sample_names(meta)=meta$SampleID
physeq=merge_phyloseq(biom,meta) 
physeq
physeq=subset_samples(physeq, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
physeq


#######
#Taxonomic Filtering
########
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-3, TRUE) #filter out any taxa lower tha 0.1%
#PhylumLevel  = transform_sample_counts(PhylumLevel, function(x) x / sum(x) )
GPrFamily=tax_glom(GPr,"Family")
FamilyLevel = filter_taxa(GPrFamily, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
#FamilyLevel  = transform_sample_counts(FamilyLevel, function(x) x / sum(x) )
GPrGenus=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GPrGenus, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%


######
#Figure 1A
######
Antimortem=subset_samples(physeq, Anti_Post== "Perimortem")

GPr  = transform_sample_counts(Antimortem, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%

TrtdataWOther<-table_for_plotting(PhylumLevel,"Phylum","Sample_Type")

TrtdataWOther$Phylum = factor(TrtdataWOther$Phylum, levels = c("Acidobacteria","Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria","Tenericutes","Other"))
TrtdataWOther$Sample_Type = factor(TrtdataWOther$Sample_Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

PhylumAntePlot2A=ggplot(TrtdataWOther, aes(x=Sample_Type,y=mean,color=Phylum))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+xlab("Sample Type")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 0,size=12),legend.position = "bottom",legend.title= element_blank())+ scale_fill_manual(values=cbPalette)+ scale_x_discrete(labels=c("Ceca" = "Ceca", "Cloaca" = "Cloaca","Large Intestine"="Large\n Intestine","Small Intestine"="Small\n Intestine"))

PhylumAntePlot2A

#####
# Figure 1B
#####
set.seed(546)
sample_data(physeq)$Sample_Type = factor(sample_data(physeq)$Sample_Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))

Antimortem=subset_samples(physeq, Anti_Post== "Perimortem")
Antimortem
sample_data(Antimortem)
ord=ordinate(Antimortem,"PCoA", "jaccard") #other common methods besides for weighted unifrac include jaccard, bray, jsd
ordplot=plot_ordination(Antimortem, ord,"samples", color="Sample_Type")+geom_point(size=4)+scale_colour_manual(name="Sample Type",values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot2B<-ordplot+ theme(legend.position=c(0.7,0.75),legend.title = element_blank(),legend.background = element_rect(color="black"))
#Remove legend title+facet_wrap(~Anti_Post)+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Anti_Post))
ordplot2B


AntimortemOrd=phyloseq::distance(Antimortem, "jaccard") #Tree will not be a submultiple as some of the taxa are subset out
adonis2(AntimortemOrd ~ Sample_Type, as(sample_data(Antimortem), "data.frame"))


AntimortemNoCloaca<-subset_samples(Antimortem,Sample_Type!="Cloaca")
GPdist=phyloseq::distance(AntimortemNoCloaca, "jaccard")
adonis2(GPdist ~ Sample_Type, as(sample_data(AntimortemNoCloaca), "data.frame"))

########
# Figure 3A
########
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%

TrtdataOther<-table_for_plotting(PhylumLevel,"Phylum","Sample_Type","Hours")


TrtdataOther$Sample_Type = factor(TrtdataOther$Sample_Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))
TrtdataOther$Hours <- as.character(TrtdataOther$Hours)
Figure3A=ggplot(TrtdataOther, aes(x=Hours,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+ylab("Relative Abundance (> 1%)") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size=12), legend.position = "bottom",legend.title= element_blank(),strip.text.x = element_text(size = 15))+
  scale_fill_manual(values=cbPalette)+facet_wrap(~Sample_Type)+xlab("Time (hrs)")+scale_x_discrete(labels=c("0_hrs" = "0","24_hrs"="24","48_hrs"="48","72_hrs"="72"))

Figure3A<-Figure3A+scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
)


#Figure3B####
#Small Intestine
SmallIntestine=subset_samples(physeq,Sample_Type=="Small Intestine")
sample_data(SmallIntestine)
SmallIntestineOrd=phyloseq::distance(SmallIntestine, "jaccard") #Tree will not be a submultiple as some of the taxa are subset out
adonis2(SmallIntestineOrd ~ Anti_Post, as(sample_data(SmallIntestine), "data.frame"))

#F1,12 = 1.87, P = 0.007

#Ceca
Ceca=subset_samples(physeq,Sample_Type=="Ceca")
sample_data(Ceca)
CecaOrd=phyloseq::distance(Ceca, "jaccard")
adonis2(CecaOrd ~ Anti_Post, as(sample_data(Ceca), "data.frame"))
#F1,12 =1.79, P = 0.016
#Large Intestine
LargeIntestine=subset_samples(physeq,Sample_Type=="Large Intestine")
sample_data(LargeIntestine)
LargeIntestineOrd=phyloseq::distance(LargeIntestine, "jaccard")
adonis2(LargeIntestineOrd ~ Anti_Post, as(sample_data(LargeIntestine), "data.frame"))
#F1,15 = 2.49, P < 0.001


#Cloaca Not run as Cloaca only has 2 antemortem samples
# Cloaca=subset_samples(physeq,Sample.Type=="Cloaca")
# sample_data(Cloaca)
# CloacaOrd=phyloseq::distance(Cloaca, "jaccard")


PERMANOVAText<-c("F = 6.06,P = 0.033","F = 2.43, P = 0.007","F = 7.70, P < 0.001"," ")
PERMANOVAText
sample_data(physeq)$Sample.Type = factor(sample_data(physeq)$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))

theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

ord=ordinate(physeq,"NMDS", "jaccard") #other common methods besides for weighted unifrac include jaccard, bray, jsd
ordplot=plot_ordination(physeq, ord,"samples", color="Anti_Post",shape="Anti_Post")+geom_point(size=2)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)+
  facet_wrap(~Sample_Type,scales = "free")+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Anti_Post))+
  theme(legend.title = element_blank(),legend.text = element_text(size = 10),legend.position = "bottom",strip.text.x = element_text(size = 15))#+  annotate("text", label = PERMANOVAText, size = 2, x = 0, y = c(-0.15))

#legend.position=c(0.73,0.14)
ordplot
Figure3B <-ordplot
Figure3B
dev.off()
tiff("Figures/Figure3B.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
Figure3B
dev.off()



######
# Phylum level plot with error bars
#####
Perimortem=subset_samples(physeq, Anti_Post =="Perimortem")
GPr  = transform_sample_counts(Perimortem, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
PhylumLevelPerimortem = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%




df <- psmelt(PhylumLevelPerimortem) #FamilyLevel
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Sample_Type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata,"TrtdataWErrorBars.csv")
Trtdata2<-read.csv("TrtdataWErrorBarsCloacaErrorRemoved.csv",header=T)
Trtdata2$Sample_Type = factor(Trtdata2$Sample_Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))
Trtdata2<-subset(Trtdata2,Phylum=="Tenericutes"|Phylum=="Bacteroidetes")
PhylumFacetPlotByType=ggplot(Trtdata2, aes(x=Sample_Type,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+ylab("Relative Abundance (> 1%, SEM)") + 
  theme(axis.title.x = element_blank(),legend.position = "none",axis.text.x = element_text(angle = 0, hjust = 0.5,size=11),strip.text.x = element_text(size = 15))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+
  scale_x_discrete(labels=c("Ceca" = "Ceca", "Cloaca" = "Cloaca","Large Intestine"="Large\n Intestine","Small Intestine"="Small\n Intestine"))+
  scale_fill_manual(values=c("#56B4E9", "#0072B2", "#D55E00", "#000000","#CC79A7")
)+facet_wrap(~Phylum)+xlab("Time (hrs)")
PhylumFacetPlotByType




#####
# Combine figure
#####
PhylumAntePlot2A

PhylumFacetPlotByType

Figure3A # Across Time

# Perimortem vs Antimortem


Figure3B # NMDS
# 
# dev.off()
# tiff("UpdatedFigMay2022.tiff", width = 200, height = 200, units = 'mm', res = 600)
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

ggarrange(PhylumAntePlot2A, PhylumFacetPlotByType,Figure3A,Figure3B,
          labels = c("a","b","c","d"),
          ncol = 2, nrow = 2)

dev.off()
