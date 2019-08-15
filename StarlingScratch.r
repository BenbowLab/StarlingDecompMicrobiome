#Start by import step in Starling2019.r
######################################################
#####Taxa bar plots, also see StarlingScript.rmd
####################################################
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
############Phylum level
sample_data(PhylumLevel)

###########################################################
##Percentages for top phyla across samples
##########################################################
df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
ggplot(Trtdata, aes(x=Phylum,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+ facet_wrap(Sample_Type~Anti_Post)+xlab("Sample Type")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_fill_manual(values=cbPalette)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)

df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata

############################
#Antimortem communities
###########################
Antimortem=subset_samples(physeq, Anti_Post== "Antemortem")
Antimortem
ord=ordinate(,"PCoA", "wunifrac") #other common methods besides for weighted unifrac include jaccard, bray, jsd
ordplot=plot_ordination(Antimortem, ord,"samples", color="Sample.Type")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   theme(legend.justification=c(1,0), legend.position=c(1,0))#+facet_wrap(~Anti_Post)+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Anti_Post))



GPr  = transform_sample_counts(Antimortem, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%
#PhylumLevel  = transform_sample_counts(PhylumLevel, function(x) x / sum(x) )
GPrFamily=tax_glom(GPr,"Family")
FamilyLevel = filter_taxa(GPrFamily, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
#FamilyLevel  = transform_sample_counts(FamilyLevel, function(x) x / sum(x) )
GPrGenus=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GPrGenus, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
GPrGenus=tax_glom(GPr,"Species")
SpeciesLevel = filter_taxa(GPrGenus, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
###############
#Alpha diversity
#################
#From Qiime outputs:
#Faith_PD_group_significance.qzv
#Shannon_group_significance.qzv



##Faceted graph of Antimortem only Phylum level
df <- psmelt(PhylumLevel)
sample_data(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Sample_Type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
PhylumAntiFacetPlot=ggplot(Trtdata, aes(x=Sample_Type,y=mean))+geom_bar(aes(fill = Sample_Type),colour="black", stat="identity")+ facet_wrap(~Phylum)+xlab("Sample Type")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
   theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+ scale_fill_manual(values=cbPalette)
PhylumAntiFacetPlot

Means=compare_means(Abundance ~ Hours, data = df, method = "kruskal.test",
                    group.by = "Phylum", p.adjust.method = "fdr")
head(Means)
keeps <- c("Phylum","group1","group2","p.format","p.adj","method","p.signif")
keeps=Means[keeps]
#keeps
test3 <- list('Phylum'= keeps$Phylum,'group1'=keeps$group1,'group2'= keeps$group2 ,'p'=keeps$p.format,'p.adj'=keeps$p.adj,p.signif=keeps$p.signif,'Method'=keeps$method)
test3= as.data.frame(test3)
test3
FilteredResults<-test3[!(test3$p.adj>0.5),]            
FilteredResults


labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
cdataplot=ggplot(Trtdata, aes(x=Anti_Post,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+ facet_grid(~Sample_Type, labeller=labeller(Sample_Type = labels))+xlab("Treatment")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)#+ scale_fill_manual(values=cbPalette)
cdataplot
#########################family level
FamilySampleType=subset_samples(FamilyLevel, Sample_Type =="Cloacal_Swab")
df <- psmelt(FamilySampleType) #FamilyLevel
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family","Hours"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
FamilyAntiFacetPlot=ggplot(Trtdata, aes(x=Hours,y=mean))+geom_bar(aes(fill = Hours),colour="black", stat="identity")+ facet_wrap(~Family)+xlab("Sample Type")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+ scale_fill_manual(values=cbPalette)
FamilyAntiFacetPlot

kruskal.test(Shannon~ Month, data=df)

Means=compare_means(Abundance ~ Hours, data = df, method= "kruskal.test",
                    group.by = "Family", p.adjust.method = "fdr")

head(Means)
keeps <- c("Family","group1","group2","p.format","p.adj","method","p.signif")
keeps=Means[keeps]
#keeps
test3 <- list('Family'= keeps$Family,'group1'=keeps$group1,'group2'= keeps$group2 ,'p'=keeps$p.format,'p.adj'=keeps$p.adj,p.signif=keeps$p.signif,'Method'=keeps$method)
test3= as.data.frame(test3)
test3
FilteredResults<-test3[!(test3$p.adj>0.5),]            
FilteredResults
################################





##################
##Compare sample types (Kruskal Wallis)
#################
#Phylum
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
df=psmelt(GPrPhylum)
Means=compare_means(Abundance ~ Sample.Type, data = df, 
                    group.by = "Phylum", p.adjust.method = "fdr")
#head(Means)
keeps <- c("Phylum","group1","group2","p.format","p.adj","method","p.signif")
keeps=Means[keeps]
#keeps
test3 <- list('Phylum'= keeps$Phylum,'group1'=keeps$group1,'group2'= keeps$group2 ,'p'=keeps$p.format,'p.adj'=keeps$p.adj,p.signif=keeps$p.signif,'Method'=keeps$method)
test3= as.data.frame(test3)
test3
FilteredResults<-test3[!(test3$p.adj>0.5),]            
FilteredResults
#Family
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrFamily=tax_glom(GPr,"Family")
df=psmelt(GPrFamily)
Means=compare_means(Abundance ~ Sample_Type, data = df, 
                    group.by = "Phylum", p.adjust.method = "fdr")
#head(Means)
keeps <- c("Phylum","group1","group2","p.format","p.adj","method","p.signif")
keeps=Means[keeps]
#keeps
test3 <- list('Phylum'= keeps$Phylum,'group1'=keeps$group1,'group2'= keeps$group2 ,'p'=keeps$p.format,'p.adj'=keeps$p.adj,p.signif=keeps$p.signif,'Method'=keeps$method)
test3= as.data.frame(test3)
test3
FilteredResults<-test3[!(test3$p.adj>0.5),]            
FilteredResults
##################################################
# logFoldChange between sample types Anti vs Postmortem
##################################################
diagdds = phyloseq_to_deseq2(GenusLevel, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GenusLevel)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`Family`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`Family` = factor(as.character(sigtab$sigtab$`Family`), levels=names(x))

sigtab$`Family` <- reorder(sigtab$`Family`, sigtab$log2FoldChange, sum)
sigtab$`FamilyGenus`= paste0(sigtab$Family,":",sigtab$Genus)
ggplot(sigtab, aes(x=`FamilyGenus`, y=log2FoldChange, color=`Family`)) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Anti vs Postmortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+ theme(legend.position="none")

##
##################################################
# logFoldChange between Sample_types antimortem  phylum
##################################################

sample_data(physeq)
#physeq2Intestine=subset_samples(physeq2, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
GPrIntestine=subset_samples(physeq, Sample.Type!="Ceca"& Sample.Type != "Cloacal_Swab")
GPrIntestine
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Antemortem")
GPrIntestineAnti
GPrIntestinePhylum <- subset_taxa(GPrIntestineAnti, Phylum != "NA")
GPrIntestinePhylum=tax_glom(GPrIntestinePhylum, "Phylum")
GPrIntestinePhylum = filter_taxa(GPrIntestinePhylum, function(x) mean(x) > 1e-2, TRUE)


physeqSampleType=subset_samples(physeq2, Sample_Type == "Small_Intestine_Tissue")
physeqSampleType <- subset_taxa(physeqSampleType, Phylum != "NA")
physeqSampleType=tax_glom(physeqSampleType, "Phylum")
sample_data(physeqSampleType)

diagdds = phyloseq_to_deseq2(physeqSampleType, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res$padj
res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeqSampleType)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab
# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`Phylum`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`Phylum` = factor(as.character(sigtab$sigtab$`Genus`), levels=names(x))

sigtab$`Genus` <- reorder(sigtab$`Genus`, sigtab$log2FoldChange, sum)
ggplot(sigtab, aes(x=`Genus`, y=log2FoldChange, color=`Phylum`)) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Large vs Small Intestine Antimortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+scale_colour_manual(values=cbPalette)#+ theme(legend.position="none")


sigtab

##################################################
# logFoldChange between  tissue type Antimortem glom Phylum
##################################################
#physeq2Intestine=subset_samples(physeq, Sample_Type!= "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
GPrIntestine=subset_samples(physeq, Sample.Type!="Cloaca"& Sample.Type != "Ceca")
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
res$padj
res
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestinePhylum)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab
# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`Genus`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`Genus` = factor(as.character(sigtab$sigtab$`Genus`), levels=names(x))

sigtab$`Genus` <- reorder(sigtab$`Genus`, sigtab$log2FoldChange, sum)
ggplot(sigtab, aes(x=`Genus`, y=log2FoldChange, color=`Phylum`)) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Large vs Small Intestine Antimortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+scale_colour_manual(values=cbPalette)#+ theme(legend.position="none")


sigtab

##################################################
# logFoldChange between  tissue type Antimortem glom Family
##################################################
physeq2Intestine=subset_samples(physeq2, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
GPrIntestine=subset_samples(physeq2Intestine, Sample_Type!="Small_Intestine_Tissue"& Sample_Type != "Large_Intestine_Tissue")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Post_Mortem")
GPrIntestineAnti
GPrIntestineFamily <- subset_taxa(GPrIntestineAnti, Family != "NA")
GPrIntestineFamily=tax_glom(GPrIntestineFamily, "Family")
GPrIntestineFamily = filter_taxa(GPrIntestineFamily, function(x) mean(x) > 1e-2, TRUE)
sample_data(GPrIntestineFamily)
diagdds = phyloseq_to_deseq2(GPrIntestineFamily, ~ Sample_Type) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res$padj
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineFamily)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab
# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`Family`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`Family` = factor(as.character(sigtab$sigtab$`Family`), levels=names(x))

sigtab$`Family` <- reorder(sigtab$`Family`, sigtab$log2FoldChange, sum)
ggplot(sigtab, aes(x=`Family`, y=log2FoldChange, color=`Phylum`)) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Large vs Small Intestine Antimortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+scale_colour_manual(values=cbPalette)#+ theme(legend.position="none")

res
sigtab



##################################################
# logFoldChange Anti vs Postmortem Phylum
##################################################
physeq2Genus=subset_samples(physeq, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
#sample_data(physeq2Genus)
GPrGenus <- subset_taxa(physeq2Genus, Phylum != "NA")
GPrIntestineFamily=tax_glom(GPrGenus, "Phylum")
GPrIntestineFamily = filter_taxa(GPrIntestineFamily, function(x) mean(x) > 1e-2, TRUE)
#sample_data(GPrIntestineFamily)
diagdds = phyloseq_to_deseq2(GPrIntestineFamily, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res$padj
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineFamily)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab
# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`Phylum`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`Family` = factor(as.character(sigtab$sigtab$`Phylum`), levels=names(x))

sigtab$`Family` <- reorder(sigtab$`Family`, sigtab$log2FoldChange, sum)
ggplot(sigtab, aes(x=`Phylum`, y=log2FoldChange, color=`Phylum`)) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Large vs Small Intestine Antimortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+scale_colour_manual(values=cbPalette)#+ theme(legend.position="none")

res
sigtab


##################################################
# logFoldChange Anti vs Postmortem Family
##################################################
physeq2Genus=subset_samples(physeq, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
#sample_data(physeq2Genus)
GPrGenus <- subset_taxa(physeq2Genus, Family != "NA")
GPrIntestineFamily=tax_glom(GPrGenus, "Family")
GPrIntestineFamily = filter_taxa(GPrIntestineFamily, function(x) mean(x) > 1e-2, TRUE)
#sample_data(GPrIntestineFamily)
diagdds = phyloseq_to_deseq2(GPrIntestineFamily, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res$padj
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineFamily)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab
# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`Family`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`Family` = factor(as.character(sigtab$sigtab$`Family`), levels=names(x))

sigtab$`Family` <- reorder(sigtab$`Family`, sigtab$log2FoldChange, sum)
ggplot(sigtab, aes(x=`Family`, y=log2FoldChange, color=`Phylum`)) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Large vs Small Intestine Antimortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+scale_colour_manual(values=cbPalette)#+ theme(legend.position="none")

res
sigtab




##################################################
# logFoldChange Anti vs Postmortem Genus
##################################################
physeq2Genus=subset_samples(physeq, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
#sample_data(physeq2Genus)
GPrGenus <- subset_taxa(physeq2Genus, Genus != "NA")
GPrIntestineFamily=tax_glom(GPrGenus, "Genus")
GPrIntestineFamily = filter_taxa(GPrIntestineFamily, function(x) mean(x) > 1e-2, TRUE)
#sample_data(GPrIntestineFamily)
diagdds = phyloseq_to_deseq2(GPrIntestineFamily, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res$padj
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineFamily)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab
# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`Family`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`Family` = factor(as.character(sigtab$sigtab$`Family`), levels=names(x))

sigtab$`Family` <- reorder(sigtab$`Family`, sigtab$log2FoldChange, sum)
ggplot(sigtab, aes(x=`Family`, y=log2FoldChange, color=`Phylum`)) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Large vs Small Intestine Antimortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+scale_colour_manual(values=cbPalette)#+ theme(legend.position="none")

res
sigtab



##################################################
# logFoldChange between Intestine tissue type Antimortem glom Genus
##################################################
physeq2Intestine=subset_samples(physeq2, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
GPrIntestine=subset_samples(physeq2Intestine, Sample_Type!="Ceca"& Sample_Type != "Cloacal_Swab")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Anti_Mortem")
GPrIntestineAnti
GPrIntestineGenus <- subset_taxa(GPrIntestineAnti, Genus != "NA")
GPrIntestineGenus=tax_glom(GPrIntestineGenus, "Genus")

sample_data(GPrIntestineGenus)
diagdds = phyloseq_to_deseq2(GPrIntestineGenus, ~ Sample_Type) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors
diagdds$TypeHours <- factor(paste0(diagdds$Sample_Type, diagdds$Hours))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)

diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
res$padj
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(GPrIntestineGenus)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab
# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`Genus`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`Genus` = factor(as.character(sigtab$sigtab$`Genus`), levels=names(x))

sigtab$`Genus` <- reorder(sigtab$`Genus`, sigtab$log2FoldChange, sum)
ggplot(sigtab, aes(x=`Genus`, y=log2FoldChange, color=`Phylum`)) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Large vs Small Intestine Antimortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+scale_colour_manual(values=cbPalette)#+ theme(legend.position="none")


sigtab
###
sample_data(physeq)
#########################################################
####Beta diversity Intestine
########################################################
physeq=merge_phyloseq(biom,meta,tree)
physeq=subset_samples(physeq, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
GPrIntestine=subset_samples(physeq, Sample_Type!="Ceca"& Sample_Type != "Cloacal_Swab")
GPrIntestineAnti=subset_samples(GPrIntestine, Anti_Post=="Anti_Mortem")
GPrIntestineAnti=rarefy_even_depth(GPrIntestineAnti, sample.size = 100) #Sets the rarefaction depth

GPdist=phyloseq::distance(GPrIntestineAnti, "wunifrac")
adonis(GPdist ~ Sample_Type, as(sample_data(GPrIntestineAnti), "data.frame"))

#####
GPrIntestine=rarefy_even_depth(GPrIntestine, sample.size = 100) #Sets the rarefaction depth
ord=ordinate(GPrIntestine,"PCoA", "wunifrac") #other common methods besides for weighted unifrac include jaccard, bray, jsd
ordplot=plot_ordination(GPrIntestine, ord,"samples", color="Sample_Type")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Sample_Type))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~Anti_Post)+ theme(legend.justification=c(1,0), legend.position=c(1,0))

###############
#Taxa plots large vs small intestine
#################
physeq4=subset_samples(physeq, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents") #Remove headgut swabs 
GPrIntestine=subset_samples(physeq4, Sample_Type!="Ceca"& Sample_Type != "Cloacal_Swab")
GPr  = transform_sample_counts(GPrIntestine, function(x) x / sum(x) ) #transform samples based on relative abundance
GPrPhylum=tax_glom(GPr,"Phylum")
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%
GPrFamily=tax_glom(GPr,"Family")
FamilyLevel = filter_taxa(GPrFamily, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
#FamilyLevel  = transform_sample_counts(FamilyLevel, function(x) x / sum(x) )
GPrGenus=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GPrGenus, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
GPrGenus=tax_glom(GPr,"Species")
SpeciesLevel = filter_taxa(GPrGenus, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 1%
PhylumLevel=subset_samples(PhylumLevel, Anti_Post =="Post_Mortem")
df <- psmelt(PhylumLevel)
head(sample_data(PhylumLevel))
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Sample_Type","Hours"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N),
)
Trtdata

cdataplot=ggplot(Trtdata, aes(x=Hours,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+xlab("")+ylab("Relative Abundance Antemortem (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+facet_wrap(~Sample_Type)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))#+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)#+ scale_fill_manual(values=cbPalette)+facet_grid(~Phylum)
cdataplot

cdataplot=ggplot(Trtdata, aes(x=Sample_Type,y=mean))+geom_bar(aes(fill = Sample_Type),colour="black", stat="identity")+xlab("")+ylab("Relative Abundance Postmortem (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+ scale_fill_manual(values=cbPalette)+facet_wrap(~Phylum)
cdataplot

#Genus level
df <- psmelt(GenusLevel)
head(sample_data(GenusLevel))
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus", "Sample_Type","Hours"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N),
)
Trtdata

cdataplot=ggplot(Trtdata, aes(x=Hours,y=mean))+geom_bar(aes(fill = Genus),colour="black", stat="identity")+xlab("")+ylab("Relative Abundance Antemortem (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+facet_wrap(~Sample_Type)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))#+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)#+ scale_fill_manual(values=cbPalette)+facet_grid(~Phylum)
cdataplot

cdataplot=ggplot(Trtdata, aes(x=Hours,y=mean))+geom_bar(aes(fill = Sample_Type),colour="black", stat="identity")+xlab("")+ylab("Relative Abundance Postmortem (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)+ scale_fill_manual(values=cbPalette)+facet_wrap(SampleType~Genus)
cdataplot




###############################
###Random Forest Intestine
###############################
ForestData=physeq#Change this one so you dont have to rewrite all variables
GenusLevel
#GenusLevel=subset_samples(GenusLevel, Anti_Post == "Post_Mortem")
#GenusLevel=subset_samples(GenusLevel, Sample_Type != "Cloacal_Swab")
taxa_names(ForestData)=c("Morganella","Serratia","Escherichia","Enterococcus","Streptococcus","Lactococcus","Clostridium")
row.names(otu_table(ForestData))=c("Morganella","Serratia","Escherichia","Enterococcus","Streptococcus","Lactococcus","Clostridium")

predictors=t(otu_table(ForestData))
row.names
(otu_table(ForestData))
row.names(otu_table(ForestData))
(otu_table(ForestData))
#taxa_names(ForestData)=c("Morganella","Serratia","Escherichia","Enterococcus","Streptococcus","Lactococcus","Clostridium")
#row.names(otu_table(ForestData))=c("Morganella","Serratia","Escherichia","Enterococcus","Streptococcus","Lactococcus","Clostridium")

dim(predictors)
response <- as.factor(sample_data(ForestData)$Hours) #This is where you change the response variable ex Hours
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
imp.20
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() #+
 # ggtitle("Most important functional groups for classifying  samples\n by hours")#\n in a string tells it to start a new line
imp.20$MeanDecreaseGini
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification








#To subset samples and only look at certain sample types 

PhylumLevelCeca=subset_samples(PhylumLevel, Sample_Type== "Ceca") #Alternatively you can use Sample_Type != "Ceca" to include everything besides Ceca
PhylumLevelCeca
df <- psmelt(PhylumLevelCeca)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Hours"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
#(Trtdata)
cdataplot=ggplot(Trtdata, aes(x=Hours,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+ facet_grid(~Phylum)+xlab("")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1)#+ scale_fill_manual(values=cbPalette)
cdataplot

#Stacked barplot of Phylum level abundance by sample location and time
df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Hours","Sample_Type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
#(Trtdata)
labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
cdataplot=ggplot(Trtdata, aes(x=Hours,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+ facet_grid(~Sample_Type, labeller=labeller(Sample_Type = labels))+xlab("Treatment")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
   scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
cdataplot
#########################################

GPdist=phyloseq::distance(physeq, "wunifrac")
adonis(GPdist ~ Sample_Type*Anti_Post, as(sample_data(physeq), "data.frame"))

Post=subset_samples(physeq, Anti_Post =="Post_Mortem")
Anti=subset_samples(physeq, Anti_Post =="Anti_Mortem")
GPdist=phyloseq::distance(Post, "wunifrac")
adonis(GPdist ~ Sample_Type, as(sample_data(Post), "data.frame"))

ord=ordinate(physeq,"PCoA", "wunifrac") 
ordplot=plot_ordination(physeq, ord,"samples", color="Anti_Post")+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Anti_Post))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   theme(legend.justification=c(1,0), legend.position=c(1,0))+facet_wrap(~Sample_Type)

#############################################################
###Beta dispersal
#################################
GPdist=phyloseq::distance(physeq, "wunifrac")
beta=betadisper(GPdist, sample_data(physeq)$Anti_Post)
permutest(beta)
boxplot(beta)

##########################
##PiCRUST 
###############################
otufullPiCRUST=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingR\\StarlingPiCrustMatrixL2.txt",header=TRUE)

metadata=read.table("C:\\Users\\Joe Receveur\\Documents\\Virtual Box\\Starling\\StarlingMeta4.18.txt",header=TRUE)
groupPiCRUST=as.matrix(read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingR\\StarlingPiCrustGroupsL2.csv"))
head(groupPiCRUST)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 20)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())) #sets the plotting theme

OTU=otu_table(otufullPiCRUST, taxa_are_rows=TRUE)
row.names(OTU)
TAX=tax_table(groupPiCRUST)

colnames(TAX)=("KEGG Pathway")
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
row.names(OTU)=taxa_names(TAX)
PiCRUST=phyloseq(OTU,TAX,sampdat)
PiCRUST=subset_samples(PiCRUST, SampleID!="ST163" &SampleID!="ST186") #remove samples with less than 2,500 total reads
PiCRUST=subset_samples(PiCRUST, Sample_Type != "Headgut_Swab" & Sample_Type != "Small_Intestine_Contents")

sample_data(PiCRUST)

PiCRUST

##########################
##PiCRUST Deseq2
###############################
PiCRUSTIntestineAnti=subset_samples(PiCRUST, Anti_Post== "Anti_Mortem")
PiCRUSTIntestineAnti=subset_samples(PiCRUSTIntestineAnti, Sample_Type!="Small_Intestine_Tissue"& Sample_Type!="Cloacal_Swab"&Sample_Type!="Headgut_Swab")
sample_data(PiCRUSTIntestineAnti)
PiCRUST = filter_taxa(PiCRUST, function(x) mean(x) > 1e-2, TRUE) 

diagdds = phyloseq_to_deseq2(PiCRUST, ~ Anti_Post) #Here is where you choose what variable you want to look at ex Anti_Post, log fold change can only deal with two groups so you will have to subset your samples
# calculate geometric means prior to estimate size factors

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
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PiCRUST)[rownames(sigtab), ], "matrix"))
head(sigtab)
sigtab
sigtab$`KEGG Pathway`
theme_set(theme_bw())



# Group log fold change plot
x = tapply(sigtab$log2FoldChange, sigtab$`KEGG Pathway`, function(x) max(x))
x = sort(x, TRUE)
sigtab$sigtab$`KEGG Pathway` = factor(as.character(sigtab$sigtab$`KEGG Pathway`), levels=names(x))

sigtab$`KEGG Pathway` <- reorder(sigtab$`KEGG Pathway`, sigtab$log2FoldChange, sum)
ggplot(sigtab, aes(x=`KEGG Pathway`, y=log2FoldChange, color=`KEGG Pathway`, )) + geom_point(size=6)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ ylab("log2FoldChange Anti vs Postmortem (SEM)")+
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1)+ theme(legend.position="none")
PiCRUSTsum  = transform_sample_counts(PiCRUST, function(x) x / sum(x) ) #transform samples based on relative abundance

PiCRUSTKEGG = filter_taxa(PiCRUSTsum, function(x) mean(x) > 1e-2, TRUE) 
PiCRUSTKEGGAnti=subset_samples(PiCRUSTKEGG, Anti_Post != "Post_Mortem")
PiCRUSTKEGGAnti
#df <- psmelt(PiCRUSTKEGGAnti)
sample_data(PiCRUSTKEGGAnti)
df <- psmelt(PiCRUSTKEGG)
df$Abundance=df$Abundance*100
head(df)
Trtdata <- ddply(df, c("KEGG.Pathway","Anti_Post"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N))
#(Trtdata)
labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
cdataplot=ggplot(Trtdata, aes(x=Anti_Post,y=mean))+geom_bar(aes(fill = KEGG.Pathway),colour="black", stat="identity")+xlab("Treatment")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette) + facet_wrap(~Sample_Type, labeller=labeller(Sample_Type = labels))
cdataplot+guides(fill=FALSE)


Trtdata <- ddply(df, c("KEGG.Pathway", "Sample_Type"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N))
#head(Trtdata)
labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
cdataplot=ggplot(Trtdata, aes(x=Anti_Post,y=mean))+geom_bar(aes(fill = Anti_Post),colour="black", stat="identity")+ facet_wrap(~KEGG.Pathway,labeller = label_wrap_gen(width=10),nrow = 2)+xlab("Treatment")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1) + scale_fill_manual(values=cbPalette)#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ 
cdataplot+guides(fill=FALSE)

###############Richness tests
Diversity=read.table("C:\\Users\\Joe Receveur\\Documents\\MSU data\\Starling\\StarlingMeta4.18WDiversity.txt",header=TRUE)
head(Diversity)
AntimortemDiv <- subset(Diversity, Sample_Type == 'Ceca')
(AntimortemDiv)
Trtdata <- ddply(Diversity, c("Anti_Post", "Sample_Type"), summarise,
                 N    = length(shannon),
                 mean = mean(shannon),
                 sd   = sd(shannon),
                 se   = sd / sqrt(N))
Trtdata
head(Trtdata)
labels <- c(Ceca = "Ceca", Cloacal_Swab = "Cloacal Swab",Large_Intestine_Tissue="Large Intestine Tissue", Small_Intestine_Contents="Small Intestine Contents",Headgut_Swab="Headgut Swab",Small_Intestine_Tissue="Small Intestine Tissue")
cdataplot=ggplot(Trtdata, aes(x=Anti_Post,y=mean))+geom_bar(aes(fill = Sample_Type),colour="black", stat="identity")+xlab("Sample Type")+ylab("Shannon Diversity") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1) + scale_fill_manual(values=cbPalette)+facet_wrap(~Sample_Type)#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ 
cdataplot+guides(fill=FALSE)

kruskal.test(shannon~ Anti_Post, data=AntimortemDiv)

Means=compare_means(shannon ~ Hours, data = AntimortemDiv, method= "wilcox",
                   , p.adjust.method = "fdr")
Means



PhylumFacetPlotByType=ggplot(dat[Abundance > 0], aes(x=Hours,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+xlab("Sample Type")+ylab("Relative Abundance (> 1%)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(angle = 0))+ scale_fill_manual(values=cbPalette)+facet_grid(~Sample.Type)
PhylumFacetPlotByType+ scale_x_continuous(breaks=c(0,24,48,72))
