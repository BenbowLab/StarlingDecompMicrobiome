
#Import######################

#Load packages
#most of the packages below can be installed with install.packages("package name here")
#For phyloseq and DeSeq2 you will have to use bioconductor, the easiest way is just to google the package name and install and follow the directions on the package's website
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
library(DESeq2)
library(tiff)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7", "#000000") #Manually defines a color blind friendly palette
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) #Sets plot options and font size ( to change font size change base_size=)

#To load files specify the direct paths
#running file.choose() will also let you pick the files you want without specifying the path to the file

biom=import_biom("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\table-with-taxonomy500.biom",parseFunction= parse_taxonomy_greengenes)
biom
set.seed(1117)
tax_table(biom)
biom=subset_taxa(biom, Kingdom== "Bacteria")
meta=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingMeta1.16.18.csv",header=TRUE)
head(meta)
#meta$Sample.Type = factor(meta$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))
meta
meta=sample_data(meta)
sample_names(meta)=meta$SampleID
tree=read.tree("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingTree7.31.19.nwk")
physeq=merge_phyloseq(biom,meta,tree)
physeq2=physeq
physeq=subset_taxa(physeq, Kingdom=="Bacteria")
physeq=subset_samples(physeq, Sample.Type != "Headgut Swab" & Sample.Type != "Small_Intestine_Contents") #Remove headgut swabs 

#physeqRar=rarefy_even_depth(physeq, sample.size = 500) #Files rarefied in QIIME see StarlingRarefactionCurves3.1.19
sample_data(physeq)
physeq



#Deseq

biom=import_biom("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\table-with-taxonomyNoR.biom",parseFunction= parse_taxonomy_greengenes)
biom
set.seed(1117)
tax_table(biom)
biom=subset_taxa(biom, Kingdom== "Bacteria")
meta=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingMeta1.16.18.csv",header=TRUE)
head(meta)
#meta$Sample.Type = factor(meta$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))
meta
meta=sample_data(meta)
sample_names(meta)=meta$SampleID
tree=read.tree("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingTree7.31.19.nwk")
physeq2=merge_phyloseq(biom,meta,tree)
physeq2=subset_taxa(physeq2, Kingdom=="Bacteria")
physeq2=subset_samples(physeq2, Sample.Type != "Headgut Swab" & Sample.Type != "Small_Intestine_Contents") #Remove headgut swabs 

#physeqRar=rarefy_even_depth(physeq, sample.size = 500) #Files rarefied in QIIME see StarlingRarefactionCurves3.1.19
sample_data(physeq2)
physeq2

########

#PiCRUSTImport#####
otufullPiCRUST=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingR\\StarlingPiCrustMatrixL2.txt",header=TRUE)

metadata=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingMeta1.16.18.csv",header=TRUE)
metadata$Sample.Type = factor(metadata$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))

groupPiCRUST=as.matrix(read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingR\\StarlingPiCrustGroupsL2.csv"))
head(groupPiCRUST)

OTU=otu_table(otufullPiCRUST, taxa_are_rows=TRUE)
#row.names(OTU)
TAX=tax_table(groupPiCRUST)

colnames(TAX)=("KEGG.Pathway")
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
row.names(OTU)=taxa_names(TAX)
PiCRUST=phyloseq(OTU,TAX,sampdat)
PiCRUST=subset_samples(PiCRUST, SampleID!="ST186"&SampleID!="ST189"&SampleID!="ST129") #remove samples that didn't pass taxonomy filtering
PiCRUST=subset_samples(PiCRUST, Sample.Type != "Headgut Swab" & Sample.Type != "Small_Intestine_Contents")
PiCRUST
head(sample_data(PiCRUST))


#Relative abundances
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
GPrGenus=tax_glom(GPr,"Species")
SpeciesLevel = filter_taxa(GPrGenus, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%

###########################################################
##Percentages for top phyla across all samples
##########################################################
#Produce a table of samples N by hours and type
metadataNumbers<-subset_samples(meta, SampleID!="ST186"&SampleID!="ST189"&SampleID!="ST129") #remove samples with < 500 reads from metadata
metadataNumbers<-subset_samples(metadataNumbers, Sample.Type != "Headgut Swab" & Sample.Type != "Small_Intestine_Contents") #Remove headgut swabs and small intestine contents from metadata

SampleNumbers <- ddply(metadataNumbers, c("Hours","Sample.Type"), summarise,
                         N    = length(LinkerSequence)
)
SampleNumbers
write.csv(SampleNumbers,"SampleNumbersStarling.csv")

#Phylum level
df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
PhylumRelAbu <- ddply(df, c("Phylum"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
PhylumRelAbu

#Family level
df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
FamilyRelAbuAll <- ddply(df, c("Family"), summarise,
                         N    = length(Abundance),
                         mean = mean(Abundance),
                         sd   = sd(Abundance),
                         se   = sd / sqrt(N)
)
FamilyRelAbuAll

#Relative Abus for antemortem samples

#Phylum
AntemortemPhylum=subset_samples(PhylumLevel, Anti_Post== "Antemortem")
sample_data(AntemortemPhylum)

df <- psmelt(AntemortemPhylum)
df$Abundance=df$Abundance*100
PhylumRelAbuAntemortem <- ddply(df, c("Phylum","Sample.Type"), summarise,
                         N    = length(Abundance),
                         mean = mean(Abundance),
                         sd   = sd(Abundance),
                         se   = sd / sqrt(N)
)
PhylumRelAbuAntemortem
write.csv(PhylumRelAbuAntemortem,"PhylumRelAbuAntemortem.csv")


#Phylum level All
df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
PhylumRelAbuAll <- ddply(df, c("Phylum","Sample.Type","Anti_Post"), summarise,
                         N    = length(Abundance),
                         mean = mean(Abundance),
                         sd   = sd(Abundance),
                         se   = sd / sqrt(N)
)
PhylumRelAbuAll
write.csv(FamilyRelAbuAll,"PhylumRelAbuAll.csv")

#Family level All
df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
FamilyRelAbuAll <- ddply(df, c("Family","Sample.Type","Anti_Post"), summarise,
                         N    = length(Abundance),
                         mean = mean(Abundance),
                         sd   = sd(Abundance),
                         se   = sd / sqrt(N)
)
head(FamilyRelAbuAll)
write.csv(FamilyRelAbuAll,"FamilyRelAbuAll.csv")
##################################################
# logFoldChange between Sample types antimortem  phylum
##################################################
#Large vs Small Intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Ceca")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, Phylum != "NA")
GPrIntestinePhylum=tax_glom(GPrIntestinePhylum, "Phylum")
GPrIntestinePhylum = filter_taxa(GPrIntestinePhylum, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestinePhylum)

diagdds = phyloseq_to_deseq2(GPrIntestinePhylum, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
res
alpha = 0.05
sigtab = res[which(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab


#Ceca vs Small Intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Large Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, Phylum != "NA")
GPrIntestinePhylum=tax_glom(GPrIntestinePhylum, "Phylum")
GPrIntestinePhylum = filter_taxa(GPrIntestinePhylum, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestinePhylum)

diagdds = phyloseq_to_deseq2(GPrIntestinePhylum, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

#Ceca vs large intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Small Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, Phylum != "NA")
GPrIntestinePhylum=tax_glom(GPrIntestinePhylum, "Phylum")
GPrIntestinePhylum = filter_taxa(GPrIntestinePhylum, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestinePhylum)

diagdds = phyloseq_to_deseq2(GPrIntestinePhylum, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab




##################################################
# logFoldChange between Sample types antimortem  phylum not filtered
##################################################
#Large vs Small Intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Ceca")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, Phylum != "NA")
GPrIntestinePhylum=tax_glom(GPrIntestinePhylum, "Phylum")
#GPrIntestinePhylum = filter_taxa(GPrIntestinePhylum, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestinePhylum)

diagdds = phyloseq_to_deseq2(GPrIntestinePhylum, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab


#Ceca vs Small Intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Large Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, Phylum != "NA")
GPrIntestinePhylum=tax_glom(GPrIntestinePhylum, "Phylum")
#GPrIntestinePhylum = filter_taxa(GPrIntestinePhylum, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestinePhylum)

diagdds = phyloseq_to_deseq2(GPrIntestinePhylum, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

#Ceca vs large intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Small Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, Phylum != "NA")
GPrIntestinePhylum=tax_glom(GPrIntestinePhylum, "Phylum")
#GPrIntestinePhylum = filter_taxa(GPrIntestinePhylum, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestinePhylum)

diagdds = phyloseq_to_deseq2(GPrIntestinePhylum, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab







##################################################
# logFoldChange between Sample types antimortem  family all
##################################################
#Large vs Small Intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Ceca")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestineFamily <- subset_taxa(GPrIntestineAnti, Family != "NA")
GPrIntestineFamily=tax_glom(GPrIntestineFamily, "Family")
#GPrIntestineFamily = filter_taxa(GPrIntestineFamily, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestineFamily)

diagdds = phyloseq_to_deseq2(GPrIntestineFamily, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineFamily)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab


#Ceca vs Small Intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Large Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestineFamily <- subset_taxa(GPrIntestineAnti, Family != "NA")
GPrIntestineFamily=tax_glom(GPrIntestineFamily, "Family")
#GPrIntestineFamily = filter_taxa(GPrIntestineFamily, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestineFamily)

diagdds = phyloseq_to_deseq2(GPrIntestineFamily, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineFamily)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

#Ceca vs large intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Small Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestineFamily <- subset_taxa(GPrIntestineAnti, Family != "NA")
GPrIntestineFamily=tax_glom(GPrIntestineFamily, "Family")
#GPrIntestineFamily = filter_taxa(GPrIntestineFamily, function(x) mean(x) > 1e-2, TRUE)

sample_data(GPrIntestineFamily)

diagdds = phyloseq_to_deseq2(GPrIntestineFamily, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineFamily)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab



##################################################
# logFoldChange between Sample types antimortem  genus all
##################################################
#Large vs Small Intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Ceca")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestineGenus <- subset_taxa(GPrIntestineAnti, Genus != "NA")
GPrIntestineGenus=tax_glom(GPrIntestineGenus, "Genus")


diagdds = phyloseq_to_deseq2(GPrIntestineGenus, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineGenus)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab


#Ceca vs Small Intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Large Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestineGenus <- subset_taxa(GPrIntestineAnti, Genus != "NA")
GPrIntestineGenus=tax_glom(GPrIntestineGenus, "Genus")


diagdds = phyloseq_to_deseq2(GPrIntestineGenus, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineGenus)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

#Ceca vs large intestine
GPrIntestine=subset_samples(physeq2, Sample.Type!="Cloaca"& Sample.Type != "Small Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestineGenus <- subset_taxa(GPrIntestineAnti, Genus != "NA")
GPrIntestineGenus=tax_glom(GPrIntestineGenus, "Genus")
#GPrIntestineFamily = filter_taxa(GPrIntestineFamily, function(x) mean(x) > 1e-2, TRUE)


diagdds = phyloseq_to_deseq2(GPrIntestineGenus, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineGenus)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab


############
#Taxonomic logfoldchange anti vs post all sample types combined
###########

AnteVSPostTax <- subset_taxa(physeq2, Phylum != "NA") #Change for other levels
AnteVSPostTax=tax_glom(AnteVSPostTax, "Phylum") #Change for other levels ex Family, Genus, Phylum


diagdds = phyloseq_to_deseq2(AnteVSPostTax, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(AnteVSPostTax)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab




################
#Taxonomic logfoldchange antemortem vs post by sample type
###############
#Large  Intestine

TaxGroup=subset_samples(physeq2, Sample.Type=="Large Intestine")
TaxGroup <- subset_taxa(TaxGroup, Genus != "NA")
TaxGroup=tax_glom(TaxGroup, "Genus")

sample_data(TaxGroup)

diagdds = phyloseq_to_deseq2(TaxGroup, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
#diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(TaxGroup)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab


#Small Intestine


TaxGroup=subset_samples(physeq2, Sample.Type=="Small Intestine")
TaxGroup <- subset_taxa(TaxGroup, Genus != "NA")
TaxGroup=tax_glom(TaxGroup, "Genus")

sample_data(TaxGroup)

diagdds = phyloseq_to_deseq2(TaxGroup, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(TaxGroup)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

#Cloaca

TaxGroup=subset_samples(physeq2, Sample.Type=="Cloaca")
TaxGroup <- subset_taxa(TaxGroup, Genus != "NA")
TaxGroup=tax_glom(TaxGroup, "Genus")

sample_data(TaxGroup)

diagdds = phyloseq_to_deseq2(TaxGroup, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(TaxGroup)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

##############
#####################
#Differences in alpha diversity antemortem
###################
#Shannon
Diversity=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingMeta4.18WDiversity.txt",header=TRUE) #File from QIIME OUTPUT
head(Diversity)
AntimortemDiv <- subset(Diversity, Anti_Post == 'Anti_Mortem'& Sample_Type!="Cloacal_Swab")
(AntimortemDiv)
Trtdata <- ddply(AntimortemDiv, c("Sample_Type"), summarise,
                 N    = length(shannon),
                 mean = mean(shannon),
                 sd   = sd(shannon),
                 se   = sd / sqrt(N))
Trtdata
labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
cdataplot=ggplot(Trtdata, aes(x=Sample_Type,y=mean))+geom_bar(aes(fill = Sample_Type),colour="black", stat="identity")+xlab("Sample Type")+ylab("Shannon Diversity") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1) + scale_fill_manual(values=cbPalette)#+facet_wrap(~Sample_Type)#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ 
cdataplot+guides(fill=FALSE)

kruskal.test(shannon~ Sample_Type, data=AntimortemDiv)

Means=compare_means(shannon ~ Sample_Type, data = AntimortemDiv, method= "wilcox",
                    , p.adjust.method = "fdr")
Means

#Faith's Diversity
head(Diversity)
Trtdata <- ddply(AntimortemDiv, c("Sample_Type"), summarise,
                 N    = length(faith_pd),
                 mean = mean(faith_pd),
                 sd   = sd(faith_pd),
                 se   = sd / sqrt(N))
Trtdata
labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
cdataplot=ggplot(Trtdata, aes(x=Sample_Type,y=mean))+geom_bar(aes(fill = Sample_Type),colour="black", stat="identity")+xlab("Sample Type")+ylab("Faith PD (SEM)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1) + scale_fill_manual(values=cbPalette)#+facet_wrap(~Sample_Type)#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ 
cdataplot+guides(fill=FALSE)

kruskal.test(faith_pd~ Sample_Type, data=AntimortemDiv)

Means=compare_means(shannon ~ Sample_Type, data = AntimortemDiv, method= "wilcox",
                    , p.adjust.method = "fdr")
Means

####################
#Differences in alpha diversity postmortem
###################
#Shannon By sample type
Diversity=read.table("StarlingMeta4.18WDiversity.txt",header=TRUE) #File from QIIME OUTPUT
head(Diversity)
PostmortemDiv <- subset(Diversity, Anti_Post == 'Post_Mortem')
head(PostmortemDiv)

kruskal.test(shannon~ Sample_Type, data=PostmortemDiv)
Means=compare_means(shannon ~ Sample_Type, data = AntimortemDiv, method= "wilcox",
                    , p.adjust.method = "fdr")
Means

#Faith's Diversity by sample type
head(Diversity)
kruskal.test(faith_pd~ Sample_Type, data=AntimortemDiv)
Means=compare_means(shannon ~ Sample_Type, data = AntimortemDiv, method= "wilcox",
                    , p.adjust.method = "fdr")
Means

#Shannon By hour Large Intestine
PostmortemDiv <- subset(Diversity, Anti_Post == 'Post_Mortem'&Sample_Type=="Large_Intestine_Tissue")
head(PostmortemDiv)

kruskal.test(shannon~ Hours, data=PostmortemDiv)

#Faith's Diversity by hour Large Intestine
kruskal.test(faith_pd~ Hours, data=PostmortemDiv)

#Shannon By hour Small Intestine
PostmortemDiv <- subset(Diversity, Anti_Post == 'Post_Mortem'&Sample_Type=="Small_Intestine_Tissue")
head(PostmortemDiv)

kruskal.test(shannon~ Hours, data=PostmortemDiv)

#Faith's Diversity by hour Small Intestine
kruskal.test(faith_pd~ Hours, data=PostmortemDiv)

#Shannon By hour Ceca
PostmortemDiv <- subset(Diversity, Anti_Post == 'Post_Mortem'&Sample_Type=="Ceca")
head(PostmortemDiv)

kruskal.test(shannon~ Hours, data=PostmortemDiv)

#Faith's Diversity by hour Ceca
kruskal.test(faith_pd~ Hours, data=PostmortemDiv)

#Shannon By hour Cloaca
PostmortemDiv <- subset(Diversity, Anti_Post == 'Post_Mortem'&Sample_Type=="Cloacal_Swab")
head(PostmortemDiv)
kruskal.test(shannon~ Hours, data=PostmortemDiv)

#Faith's Diversity by hour Cloaca
kruskal.test(faith_pd~ Hours, data=PostmortemDiv)

######################
#Differences in alpha diversity antemortem vs post split by sample type 
###################
#Shannon Ante-Post Large Intestine
PostmortemDiv <- subset(Diversity,Sample_Type=="Large_Intestine_Tissue")
head(PostmortemDiv)

kruskal.test(shannon~ Anti_Post, data=PostmortemDiv)
Trtdata <- ddply(PostmortemDiv, c("Anti_Post"), summarise,
                 N    = length(shannon),
                 mean = mean(shannon),
                 sd   = sd(shannon),
                 se   = sd / sqrt(N))
Trtdata



#Faith's Diversity Ante-Post Large Intestine
kruskal.test(faith_pd~ Anti_Post, data=PostmortemDiv)
Trtdata <- ddply(PostmortemDiv, c("Anti_Post"), summarise,
                 N    = length(faith_pd),
                 mean = mean(faith_pd),
                 sd   = sd(faith_pd),
                 se   = sd / sqrt(N))
Trtdata
#Shannon By Amte-Post Small Intestine
PostmortemDiv <- subset(Diversity, Sample_Type=="Small_Intestine_Tissue")
head(PostmortemDiv)

kruskal.test(shannon~ Anti_Post, data=PostmortemDiv)
Trtdata <- ddply(PostmortemDiv, c("Anti_Post"), summarise,
                 N    = length(shannon),
                 mean = mean(shannon),
                 sd   = sd(shannon),
                 se   = sd / sqrt(N))
Trtdata
#Faith's Diversity by Ante-Post Small Intestine
kruskal.test(faith_pd~ Anti_Post, data=PostmortemDiv)
Trtdata <- ddply(PostmortemDiv, c("Anti_Post"), summarise,
                 N    = length(faith_pd),
                 mean = mean(faith_pd),
                 sd   = sd(faith_pd),
                 se   = sd / sqrt(N))
Trtdata
#Shannon By Ante-Post Ceca
PostmortemDiv <- subset(Diversity, Sample_Type=="Ceca")
head(PostmortemDiv)

kruskal.test(shannon~ Anti_Post, data=PostmortemDiv)
Trtdata <- ddply(PostmortemDiv, c("Anti_Post"), summarise,
                 N    = length(shannon),
                 mean = mean(shannon),
                 sd   = sd(shannon),
                 se   = sd / sqrt(N))
Trtdata
#Faith's Diversity by Ante-Post Ceca
kruskal.test(faith_pd~ Anti_Post, data=PostmortemDiv)
Trtdata <- ddply(PostmortemDiv, c("Anti_Post"), summarise,
                 N    = length(faith_pd),
                 mean = mean(faith_pd),
                 sd   = sd(faith_pd),
                 se   = sd / sqrt(N))
Trtdata
#######################
#Beta dispersion split by sample type
#########################
#All grouped Antemortem vs post 
GPdist=phyloseq::distance(physeq, "jaccard")
beta=betadisper(GPdist, sample_data(physeq)$Anti_Post)
permutest(beta)
boxplot(beta)
beta

TukeyHSD(beta, which = "group", ordered = FALSE,
         conf.level = 0.95)
#Large Intestine
Subset<-subset_samples(physeq,Sample.Type=="Large Intestine")
Subset
GPdist=phyloseq::distance(Subset, "jaccard")
beta=betadisper(GPdist, sample_data(Subset)$Anti_Post)
permutest(beta)
boxplot(beta)
beta

TukeyHSD(beta, which = "group", ordered = FALSE,
         conf.level = 0.95)

#Small Intestine
Subset<-subset_samples(physeq,Sample.Type=="Small Intestine")
sample_data(Subset)
GPdist=phyloseq::distance(Subset, "jaccard")
beta=betadisper(GPdist, sample_data(Subset)$Anti_Post)
permutest(beta)
boxplot(beta)
beta

TukeyHSD(beta, which = "group", ordered = FALSE,
         conf.level = 0.95)

#Ceca
Subset<-subset_samples(physeq,Sample.Type=="Ceca")
sample_data(Subset)
GPdist=phyloseq::distance(Subset, "jaccard")
beta=betadisper(GPdist, sample_data(Subset)$Anti_Post)
permutest(beta)
boxplot(beta)
beta

TukeyHSD(beta, which = "group", ordered = FALSE,
         conf.level = 0.95)
########################
#Differences in beta diversity antemortem
#######################
Antimortem=subset_samples(physeq, Anti_Post== "Antemortem")
AntimortemJaccard=subset_samples(Antimortem, Sample.Type!= "Cloacal_Swab")

GPdistAntemortemJaccard=phyloseq::distance(AntimortemJaccard, "jaccard")
adonis(GPdistAntemortemJaccard ~ Sample.Type, as(sample_data(AntimortemJaccard), "data.frame"))
#F = 1.116, P = 0.125
#####################
#Differences in beta diversity antemortem vs postmortem by sample type
#####################
#Small Intestine
SmallIntestine=subset_samples(physeq,Sample.Type=="Small Intestine")
sample_data(SmallIntestine)
SmallIntestineOrd=phyloseq::distance(SmallIntestine, "jaccard") #Tree will not be a submultiple as some of the taxa are subset out
adonis(SmallIntestineOrd ~ Anti_Post, as(sample_data(SmallIntestine), "data.frame"))
#F= 1.87, P = 0.011
#Large Intestine
LargeIntestine=subset_samples(physeq,Sample.Type=="Large Intestine")
sample_data(LargeIntestine)
LargeIntestineOrd=phyloseq::distance(LargeIntestine, "jaccard")
adonis(LargeIntestineOrd ~ Anti_Post, as(sample_data(LargeIntestine), "data.frame"))
# F = 2.49, P = 0.001
#Cloaca Not run as Cloaca only has 2 antemortem samples
Cloaca=subset_samples(physeq,Sample.Type=="Cloaca")
sample_data(Cloaca)
CloacaOrd=phyloseq::distance(Cloaca, "jaccard")
adonis(CloacaOrd ~ Anti_Post, as(sample_data(Cloaca), "data.frame"))
#Ceca
Ceca=subset_samples(physeq,Sample.Type=="Ceca")
sample_data(Ceca)
CecaOrd=phyloseq::distance(Ceca, "jaccard")
adonis(CecaOrd ~ Anti_Post, as(sample_data(Ceca), "data.frame"))


#############################
#Picrust rel abundances antemortem
##############################
PicrustAntemortem<-subset_samples(PiCRUST, Anti_Post== "Antemortem")
PicrustAntemortem<- transform_sample_counts(PicrustAntemortem, function(x) x / sum(x) ) #transform samples based on relative abundance

head(sample_data(PicrustAntemortem))
df <- psmelt(PicrustAntemortem)
df$Abundance=df$Abundance*100
head(df)
Trtdata <- ddply(df, c("KEGG.Pathway","Sample.Type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N))
(Trtdata)
#write.csv(Trtdata, "AntemortemPiCrustByType.csv")

###########################
#PiCRUST logfold change sample types Antemortem
############################
#Large vs Small Intestine
PicrustAntemortem<-subset_samples(PiCRUST, Anti_Post== "Antemortem")

PicrustIntestinentemortem=subset_samples(PicrustAntemortem, Sample.Type!="Cloaca"& Sample.Type != "Ceca")
PicrustIntestinentemortem
PicrustIntestinentemortem <- subset_taxa(PicrustIntestinentemortem, KEGG.Pathway != "NA")

PicrustIntestinentemortem

diagdds = phyloseq_to_deseq2(PicrustIntestinentemortem, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PicrustIntestinentemortem)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab


#Ceca vs Small Intestine
GPrIntestine<-PicrustAntemortem #Change this instead of changing all variables below
GPrIntestine=subset_samples(GPrIntestine, Sample.Type!="Cloaca"& Sample.Type != "Large Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, KEGG.Pathway != "NA")

sample_data(GPrIntestinePhylum)

diagdds = phyloseq_to_deseq2(GPrIntestinePhylum, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

#Ceca vs large intestine
GPrIntestine<-PicrustAntemortem #Change this instead of changing all variables below
GPrIntestine=subset_samples(GPrIntestine, Sample.Type!="Cloaca"& Sample.Type != "Small Intestine")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, KEGG.Pathway != "NA")


sample_data(GPrIntestinePhylum)

diagdds = phyloseq_to_deseq2(GPrIntestinePhylum, ~ Type_Hours) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab



###########################
#PiCRUST logfold change by sample types Anti vs post
############################
#Large  Intestine

PiCRUStGroup=subset_samples(PiCRUST, Sample.Type=="Large Intestine")
PiCRUStGroup <- subset_taxa(PiCRUStGroup, KEGG.Pathway != "NA")

sample_data(PiCRUStGroup)

diagdds = phyloseq_to_deseq2(PiCRUStGroup, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
#diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PiCRUStGroup)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab


#Small Intestine


PiCRUStGroup=subset_samples(PiCRUST, Sample.Type=="Small Intestine")
PiCRUStGroup <- subset_taxa(PiCRUStGroup, KEGG.Pathway != "NA")

sample_data(PiCRUStGroup)

diagdds = phyloseq_to_deseq2(PiCRUStGroup, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PiCRUStGroup)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

#Cloaca

PiCRUStGroup=subset_samples(PiCRUST, Sample.Type=="Cloaca")
PiCRUStGroup <- subset_taxa(PiCRUStGroup, KEGG.Pathway != "NA")

sample_data(PiCRUStGroup)

diagdds = phyloseq_to_deseq2(PiCRUStGroup, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PiCRUStGroup)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab

##################
#Picrust logfold change all groups antemortem vs post
#################
#PiCRUST logfold change by sample types Anti vs post
############################
#All samples antemortem vs post

PiCRUStGroup=PiCRUST
PiCRUStGroup <- subset_taxa(PiCRUStGroup, KEGG.Pathway != "NA")

sample_data(PiCRUStGroup)

diagdds = phyloseq_to_deseq2(PiCRUStGroup, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
#diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
#res$padj
#res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PiCRUStGroup)[rownames(sigtab), ], "matrix"))
#head(sigtab)
sigtab



##############################
#Picrust relative abu anti vs post
##############################
PiCRUSTTransformed<- transform_sample_counts(PiCRUST, function(x) x / sum(x) ) #transform samples based on relative abundance

df <- psmelt(PiCRUSTTransformed)
df$Abundance=df$Abundance*100
head(df)
Trtdata <- ddply(df, c("KEGG.Pathway","Sample.Type","Anti_Post"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N))
(Trtdata)
write.csv(Trtdata, "PiCrustByTypeAntiVsPost.csv")

#####################
#Overall Picrust pathway abundance
####################
PiCRUSTTransformed<- transform_sample_counts(PiCRUST, function(x) x / sum(x) ) #transform samples based on relative abundance

df <- psmelt(PiCRUSTTransformed)
df$Abundance=df$Abundance*100
head(df)
Trtdata <- ddply(df, c("KEGG.Pathway"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N))
(Trtdata)
write.csv(Trtdata, "PiCrustRelAbuAll.csv")
############################
#Figures
############################
#Figure2A######################
Antimortem=subset_samples(physeq, Anti_Post== "Antemortem")

GPr  = transform_sample_counts(Antimortem, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%


df <- psmelt(PhylumLevel)
head(sample_data(PhylumLevel))
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Sample.Type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
write.csv(Trtdata, file = "AntemortemPhylumLevelWOOther.csv")
TrtdataWOther=read.csv( "C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingR\\AntemortemPhylumLevelIncludingOther.csv",header=TRUE)#Added in other category in excel to sum to 100%

levels(TrtdataWOther$Phylum)
TrtdataWOther$Phylum = factor(TrtdataWOther$Phylum, levels = c("Acidobacteria","Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria","Tenericutes","Other"))
TrtdataWOther$Sample.Type = factor(TrtdataWOther$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))

PhylumAntePlot2A=ggplot(TrtdataWOther, aes(x=Sample.Type,y=mean,color=Phylum))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+xlab("Sample Type")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 0,size=6))+ scale_fill_manual(values=cbPalette)+ scale_x_discrete(labels=c("Ceca" = "Ceca", "Cloaca" = "Cloaca","Large Intestine"="Large\n Intestine","Small Intestine"="Small\n Intestine"))+
  theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) #sets the plotting theme
PhylumAntePlot2A
dev.off()
tiff("Figures/Figure2A.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
PhylumAntePlot2A
dev.off()

#Figure2B#############
sample_data(physeq)$Sample.Type = factor(sample_data(physeq)$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))

Antimortem=subset_samples(physeq, Anti_Post== "Antemortem")
Antimortem
sample_data(Antimortem)
ord=ordinate(Antimortem,"PCoA", "jaccard") #other common methods besides for weighted unifrac include jaccard, bray, jsd
ordplot=plot_ordination(Antimortem, ord,"samples", color="Sample.Type")+geom_point(size=4)+scale_colour_manual(name="Sample Type",values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot2B<-ordplot+ theme(legend.position=c(0.7,0.75),legend.title = element_blank(),legend.background = element_rect(color="black"))
#Remove legend title+facet_wrap(~Anti_Post)+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Anti_Post))
ordplot2B

dev.off()
tiff("Figures/Figure2B.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
ordplot2B
dev.off()

#Figure2C######
PiCRUSTAntimortem=subset_samples(PiCRUST, Anti_Post== "Antemortem")

GPr  = transform_sample_counts(PiCRUSTAntimortem, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPiCRUST=tax_glom(GPr,"KEGG.Pathway")
KEGGLevel = filter_taxa(GPrPiCRUST, function(x) mean(x) > 5e-2, TRUE) #filter out any taxa lower tha 0.1%

df <- psmelt(KEGGLevel)
df$Abundance=df$Abundance*100

head(tax_table(KEGGLevel))
Trtdata <- ddply(df, c("KEGG.Pathway", "Sample.Type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
#write.csv(Trtdata, file = "AntemortemPiCRUSTLevelWOOther.csv")
TrtdataWOther=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingR\\AntemortemPiCRUSTLevelIncludingOther.csv",header=TRUE)

levels(TrtdataWOther$KEGG.Pathway)

TrtdataWOther$KEGG.Pathway = factor(TrtdataWOther$KEGG.Pathway, levels = c("Amino Acid Metabolism","Carbohydrate Metabolism","Cellular Processes and Signaling","Energy Metabolism","Lipid Metabolism","Membrane Transport","Metabolism of Cofactors and Vitamins","Nucleotide Metabolism","Poorly Characterized","Replication and Repair","Translation","Xenobiotics Biodegradation and Metabolism","Other"))

#PiCRUSTAntePlot=ggplot(TrtdataWOther, aes(x=Sample.Type,y=mean,color=KEGG.Pathway))+geom_bar(aes(fill = KEGG.Pathway),colour="black", stat="identity")+xlab("Sample Type")+ylab("Relative Abundance (> 3%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
#  theme(axis.text.x = element_text(angle = 0))+labs(fill="KEGG Pathway")#+scale_fill_manual(name="Sample Type")#+ scale_fill_manual(values=cbPalette)
#PhylumAntePlot

Trtdata$Sample.Type = factor(Trtdata$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))
Trtdata2C<-Trtdata
PiCrust2C=ggplot(Trtdata2C, aes(x=Sample.Type,y=mean))+ theme_bw(base_size = 7)+geom_bar(aes(fill = Sample.Type),colour="black", stat="identity")+xlab("Sample Type")+
  ylab("Relative Pathway Abundance (> 5%, SEM)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+facet_wrap(~KEGG.Pathway)+ 
  scale_fill_manual(values=cbPalette)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+
  guides(fill=FALSE)#+
   
PiCrust2C
dev.off()
tiff("Figures/Figure2C.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
PiCrust2C
dev.off()

#Figure2D####
sample_data(PiCRUSTAntimortem)$Sample.Type = factor(sample_data(PiCRUSTAntimortem)$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))

#theme_set(theme_bw(base_size = 8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) #sets the plotting theme

PiCRUSTAntimortem
ord=ordinate(PiCRUSTAntimortem,"PCoA", "jaccard") #other common methods besides for weighted unifrac include jaccard, bray, jsd
ordplot=plot_ordination(PiCRUSTAntimortem, ord,"samples", color="Sample.Type")+geom_point(size=4)+scale_colour_manual(name="Sample Type",values=cbPalette)+scale_fill_manual(values=cbPalette)
PiCrust2D<-ordplot+ theme_bw(base_size = 8)+theme(legend.position=c(0.8,0.2),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title = element_blank(),legend.background = element_rect(color="black"))+guides(fill=F)

dev.off()
tiff("Figures/Figure2D.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
PiCrust2D
dev.off()


#######
#Join Together Figure 2
#######

PhylumAntePlot2A
ordplot2B
PiCrust2C
PiCrust2D

dev.off()
tiff("Figures/Figure2.tiff", width = 6.85, height = 6.85, units = 'in', res = 300)
ggarrange(PhylumAntePlot2A,ordplot2B,PiCrust2C,PiCrust2D,
          labels = c("a", "b", "c","d"),
          ncol = 2, nrow = 2)
dev.off()


#Figure3A#####

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%

#PhylumSampleType=subset_samples(PhylumLevel, Sample_Type =="Cloacal_Swab")
df <- psmelt(PhylumLevel) #FamilyLevel
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum","Hours", "Sample.Type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
write.csv(Trtdata,file="PhylumLevelByHourAndSampleType.csv")
TrtdataOther<-read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingR\\PhylumLevelByHourAndSampleTypeIncludingOther.csv",header=T)#Added other category to make stackd bar graphs reach 100%
TrtdataOther$Sample.Type = factor(TrtdataOther$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))
TrtdataOther$Phylum = factor(TrtdataOther$Phylum, levels = c("Actinobacteria","Firmicutes","Proteobacteria","Tenericutes","Other"))
TrtdataOther$Hours <- as.character(TrtdataOther$Hours)
PhylumFacetPlotByType=ggplot(TrtdataOther, aes(x=Hours,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+ylab("Relative Abundance (> 1%)") + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
 scale_fill_manual(values=cbPalette)+facet_grid(~Sample.Type)+xlab("Time (hrs)")
Figure3A<-PhylumFacetPlotByType+  theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #sets the plotting theme
Figure3A
dev.off()
tiff("Figures/Figure3A.tiff", width = 6.6, height = 3.3, units = 'in', res = 300)
Figure3A
dev.off()

#Figure3B####
#Small Intestine
SmallIntestine=subset_samples(physeq,Sample.Type=="Small Intestine")
sample_data(SmallIntestine)
SmallIntestineOrd=phyloseq::distance(SmallIntestine, "jaccard") #Tree will not be a submultiple as some of the taxa are subset out
adonis(SmallIntestineOrd ~ Anti_Post, as(sample_data(SmallIntestine), "data.frame"))

#F1,12 = 1.87, P = 0.007

#Ceca
Ceca=subset_samples(physeq,Sample.Type=="Ceca")
sample_data(Ceca)
CecaOrd=phyloseq::distance(Ceca, "jaccard")
adonis(CecaOrd ~ Anti_Post, as(sample_data(Ceca), "data.frame"))
#F1,12 =1.79, P = 0.016
#Large Intestine
LargeIntestine=subset_samples(physeq,Sample.Type=="Large Intestine")
sample_data(LargeIntestine)
LargeIntestineOrd=phyloseq::distance(LargeIntestine, "jaccard")
adonis(LargeIntestineOrd ~ Anti_Post, as(sample_data(LargeIntestine), "data.frame"))
#F1,15 = 2.49, P < 0.001


#Cloaca Not run as Cloaca only has 2 antemortem samples
Cloaca=subset_samples(physeq,Sample.Type=="Cloaca")
sample_data(Cloaca)
CloacaOrd=phyloseq::distance(Cloaca, "jaccard")


PERMANOVAText<-c("F = 6.06,P = 0.033","F = 2.43, P = 0.007","F = 7.70, P < 0.001"," ")
PERMANOVAText
sample_data(physeq)$Sample.Type = factor(sample_data(physeq)$Sample.Type, levels = c("Small Intestine","Ceca", "Large Intestine","Cloaca"))

ord=ordinate(physeq,"PCoA", "jaccard") #other common methods besides for weighted unifrac include jaccard, bray, jsd
ordplot=plot_ordination(physeq, ord,"samples", color="Anti_Post",shape="Anti_Post")+geom_point(size=2)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)+
  facet_wrap(~Sample.Type,scales = "free")+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Anti_Post))+
  theme(legend.position=c(0.73,0.14),legend.title = element_blank(),legend.text = element_text(size = 6),legend.background = element_rect(color="black"))#+  annotate("text", label = PERMANOVAText, size = 2, x = 0, y = c(-0.15))
ordplot
Figure3B <-ordplot
Figure3B
dev.off()
tiff("Figures/Figure3B.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
Figure3B
dev.off()

#Figure3C####

PiCRUStGroup=PiCRUST
PiCRUStGroup <- subset_taxa(PiCRUStGroup, KEGG.Pathway != "NA")

head(sample_data(PiCRUStGroup))

diagdds = phyloseq_to_deseq2(PiCRUStGroup, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
res$padj
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PiCRUStGroup)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab$RelABu=sigtab$baseMean/sum(sigtab$baseMean)
sigtab
sigtabFivePercent<-subset.data.frame(sigtab,RelABu>0.05)
sigtab<-sigtabFivePercent
# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`KEGG.Pathway`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`KEGG.Pathway` = factor(as.character(sigtab$sigtab$`KEGG.Pathway`), levels=names(x))

sigtab$`KEGG.Pathway` <- reorder(sigtab$`KEGG.Pathway`, sigtab$log2FoldChange, sum)

DeseqPlot<-ggplot(sigtab, aes(x=KEGG.Pathway, y=log2FoldChange, color=`KEGG.Pathway` )) + geom_point(size=3)+
  ylab("log2FoldChange Anti vs Postmortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+xlab("KEGG Pathway")+
  scale_colour_manual(values=cbPalette)#+ theme(legend.position="none")
DeseqPlot
Figure3C<-DeseqPlot+ geom_hline(yintercept=0, linetype="dashed", color = "black",size=1)+geom_text(x=4,y=0.1, label="Higher Postmortem",color="black",size=5)+geom_text(x=4,y=-0.05,size=5,label="Higher Antemortem",color="black")+ 
  scale_x_discrete(labels=c("Amino Acid Metabolism" = "AA\n Metabolism", "Poorly Characterized" = "Poorly\n Characterized","Carbohydrate Metabolism" = "Carbohydrate\n Metabolism","Cellular Processes and Signaling"="Cell Processes \n & Signaling","Membrane Transport"="Membrane\n Transport"))+
  theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none",axis.title.x=element_blank(),axis.text.x = element_text(size = 5.5))+guides(Fill=F) #sets the plotting theme

dev.off()
tiff("Figures/Figure3C.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
Figure3C
dev.off()
#775x775

#Figure4A####

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrGenus=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GPrGenus, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%

df <- psmelt(GenusLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus", "Sample.Type","Anti_Post"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Subset<-subset.data.frame(Trtdata,Genus=="Serratia"|Genus=="Lactococcus"|Genus=="[Clostridium]")
Subset

Subset$Combined<-paste0(Subset$Sample.Type,": ",Subset$Anti_Post)
Subset
GenusFacetPlot=ggplot(Subset, aes(x=Sample.Type,y=mean))+geom_bar(aes(fill = Anti_Post),colour="black", stat="identity")+ facet_wrap(Anti_Post~Genus)+xlab("Sample Type")+ylab("Relative Abundance (%, SEM)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+ scale_fill_manual(values=cbPalette)+
  theme(legend.position=c(0.16,.87))+scale_colour_manual(values=cbPalette)+theme(legend.position = "none")
Figure4A<-GenusFacetPlot

dev.off()
tiff("Figures/Figure4A.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
Figure4A
dev.off()

#Figure4B####

GPr  = transform_sample_counts(PiCRUST, function(x) x / sum(x) ) #transform samples based on relative abundance
#FilteredPathway = filter_taxa(GPr, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%

df <- psmelt(GPr)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("KEGG.Pathway", "Sample.Type","Anti_Post"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
(Trtdata)
Subset<-subset.data.frame(Trtdata,KEGG.Pathway=="Cell Motility"|KEGG.Pathway=="Genetic Information Processing")
Subset
Subset$KEGG.Pathway


Subset$Combined<-paste0(Subset$Sample.Type,": ",Subset$Anti_Post)
Subset
GenusFacetPlot=ggplot(Subset, aes(x=Sample.Type,y=mean))+geom_bar(aes(fill = Anti_Post),colour="black", stat="identity")+ facet_wrap(Anti_Post~KEGG.Pathway)+xlab("Sample Type")+ylab("Relative Abundance (%,SEM)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+ scale_fill_manual(values=cbPalette)+
  theme(legend.position=c(0.85,.85))+scale_colour_manual(values=cbPalette)+theme(legend.position = "none")
Figure4B<-GenusFacetPlot
dev.off()
tiff("Figures/Figure4BNew.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
Figure4B
dev.off()


#Join together figures


dev.off()
tiff("Figures/Figure4.tiff", width = 6.6, height = 3.3, units = 'in', res = 300)
ggarrange(Figure4A, Figure4B,
          labels = c("a","b"),
          ncol = 2, nrow = 1)
dev.off()





#Supplemental FIgures

###########
#Alpha diversity ante vs post split by sample
##########
Diversity=read.csv("StarlingMeta4.18WDiversity.csv",header=TRUE) #File from QIIME OUTPUT
head(Diversity)

Trtdata <- ddply(Diversity, c("Sample_Type","Anti_Post"), summarise,
                 N    = length(shannon),
                 mean = mean(shannon),
                 sd   = sd(shannon),
                 se   = sd / sqrt(N))
Trtdata

ShannonAntePost=ggplot(Trtdata, aes(x=Anti_Post,y=mean))+geom_bar(aes(fill = Sample_Type),colour="black", stat="identity")+xlab("Sample Type")+ylab("Shannon Diversity") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1) + scale_fill_manual(values=cbPalette)+facet_wrap(~Sample_Type)+guides(fill=FALSE)
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) #Sets plot options and font size ( to change font size change base_size=)

dev.off()
tiff("Figures/FigureSXShannonAntePostByType.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
ShannonAntePost
dev.off()

#Faith's PD
Trtdata <- ddply(Diversity, c("Sample_Type","Anti_Post"), summarise,
                 N    = length(faith_pd),
                 mean = mean(faith_pd),
                 sd   = sd(faith_pd),
                 se   = sd / sqrt(N))
Trtdata

FaithAntePost=ggplot(Trtdata, aes(x=Anti_Post,y=mean))+geom_bar(aes(fill = Sample_Type),colour="black", stat="identity")+xlab("Sample Type")+ylab("Faith's PD") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1) + scale_fill_manual(values=cbPalette)+facet_wrap(~Sample_Type)+guides(fill=FALSE)
theme_set(theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) #Sets plot options and font size ( to change font size change base_size=)

dev.off()
tiff("Figures/FigureSXFaithAntePostByType.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
FaithAntePost
dev.off()


dev.off()
tiff("Figures/AlphaSupplemental.tiff", width = 6.85, height = 6.85, units = 'in', res = 300)
ggarrange(ShannonAntePost,FaithAntePost,
          labels = c("a", "b"),
          ncol = 2, nrow = 1)
dev.off()
