library(dplyr)
library(vegan)
library(permute)
library(car)
library(MASS)
library(tidyverse)
library (ggplot2)
library(RColorBrewer)
library(ape)
library(MuMIn)
library(ggsignif)
library(asbio)
library(doBy)
library(emmeans)
library(multcomp)
library(gridExtra)
library(rcompanion)

# set working directory

setwd('/Users/will1809/OneDrive - purdue.edu/Repositories.GitHub/WilliamsGinzel.JnGmFusarium.2021')

source("Code/OTU.functions.R")

# read in clusters and assignments from mothur

OTUS.mothur <- read.csv("Data/mothur.otus.ForestPlant.Nov20.2020.csv", header=F)
names(OTUS.mothur) <- c('otu','Isolate')
OTUS <- convert.otu.table_(OTUS.mothur)
OTUS$Isolate <- OTUS$Isolate %>% gsub('(?<=Rh)_', '', ., perl=T)

# read in experiment metadata

design <- read.csv('Data/Isolate.metadata.ForestPlant.4.9.19.SC.csv')

inoc.design <- read.csv('Data/forest plant smbiome treat design.csv') %>%
	mutate(Soil=str_replace_all(string=Soil, pattern="Sterile", replacement="Control"))

inoc.design$plant <- paste(inoc.design$Soil, inoc.design$Rep)
names(inoc.design)[names(inoc.design)=="Treatment"] <- "Inoc"

inoc.design<-inoc.design[-c(24,67, 103,107),]

# the following code blocks manipulates the data to sort it as desired and rename codes

pattern <- "(?<=^Rh)[0-9]{1,3}(?=([ab])|(_[12]))?"

design$Parent.isolate <- regmatches(design$Isolate, regexpr(pattern, design$Isolate, perl=T))

OTUS$Parent.isolate <- regmatches(OTUS$Isolate, regexpr(pattern, OTUS$Isolate, perl=T))

fulld <- na.omit(plyr::join(OTUS, design[,3:6], type="right",by="Parent.isolate", match = "all"))
dim(design)
dim(fulld)

fulld<-fulld[-10,]
partd <- fulld %>%
mutate(Treat=str_replace_all(string = Treat, pattern='S', 'Control')) %>%
mutate(Treat=str_replace_all(string = Treat, pattern='P', 'Plantation')) %>%
mutate(Treat=str_replace_all(string = Treat, pattern='F', 'Forest'))

partd$Treat <- factor(partd$Treat, levels=c('Control','Plantation','Forest'))

partd$plant <- paste(partd$Treat, partd$Rep)

as.factor(design[,"Isolate"]) %>% setdiff(as.factor(partd$Isolate))

# look at taxon summary by treatment

sum(partd$Treat=="Forest")
sum(partd$Treat=="Plantation")
sum(partd$Treat=="Control")

# rename taxonomy and sort it into things

taxids <- read_tsv("Data/ForPlant.ITS.MesquiteAssembly.Final.mothur.unique.opti_mcc.0.05.cons.taxonomy") %>% dplyr::rename(otu=OTU) %>% dplyr::rename(taxonomy=Taxonomy) %>% 
	mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
	mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
	mutate(taxonomy=str_replace_all(string=taxonomy, pattern="[kpocfgs]\\__", replacement="")) %>%
	mutate(taxonomy=str_replace_all(string=taxonomy, pattern="_unclassified", replacement="_sp")) %>%
	separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep=";")

taxids$otu <- as.factor(taxids$otu)

# quality check the data

taxids %>% anti_join(., fulld) 

taxids %>% inner_join(., fulld) %>% write.csv('OTU.table.11.30.20.csv')

# look at the relative abundances of some of the more common Nectriaceae pathogen genera

data.frame(isolate=c((taxids %>% inner_join(., fulld))$Isolate[grep('Fusarium', (taxids %>% inner_join(., fulld))$genus)],
(taxids %>% inner_join(., fulld))$Isolate[grep('Neocosmospora', (taxids %>% inner_join(., fulld))$genus)],
(taxids %>% inner_join(., fulld))$Isolate[grep('Cylindrocarpon', (taxids %>% inner_join(., fulld))$genus)],
(taxids %>% inner_join(., fulld))$Isolate[grep('Nectriaceae_sp', (taxids %>% inner_join(., fulld))$genus)]))

# look at rel abunds of Trichoderma spp.

data.frame(isolate=(taxids %>% inner_join(., fulld))$Isolate[grep('Trichoderma', (taxids %>% inner_join(., fulld))$genus)])

# select out taxa we want to work with (more than one instance, not a contaminant)
# we assume Malasseziomycetes is a contaminant
non.singletons <- taxids %>% inner_join(., partd) %>% filter (class != "Malasseziomycetes") %>% group_by (otu) %>% filter (n() > 1) %>% ungroup()

# organize by genus
otu.by.genus <- non.singletons %>% filter (!is.na(genus)) %>%
	group_by (Treat, genus) %>% dplyr::summarize (count = n()) %>%
	arrange(Treat, (desc(count))) #%>% dplyr::top_n(8, count)

# look at unique OTUs and genera
unique(taxids$otu)
unique(otu.by.genus$genus)

# create data frames at different taxonomic levels
# order, family
# then summarize families, genera, and OTUs

otu.by.order <- taxids %>% inner_join(., partd) %>% filter (!is.na(order)) %>%
	filter (class != "Malasseziomycetes") %>%
	group_by (Treat, order) %>% dplyr::summarize (count = n()) %>%
	arrange(Treat, (desc(count)))

otu.by.family <- taxids %>% inner_join(., partd) %>% filter (!is.na(family)) %>%
	filter (class != "Malasseziomycetes") %>%
	group_by (Treat, family) %>% dplyr::summarize (count = n()) %>%
	arrange(Treat, (desc(count))) %>% dplyr::top_n(10, count)

taxids %>% inner_join(., partd) %>% filter (!is.na(family)) %>%
	filter (class != "Malasseziomycetes") %>%
	group_by (Treat, family) %>% dplyr::summarize (count = n()) %>% summary()

taxids %>% inner_join(., partd) %>% filter (!is.na(genus)) %>%
	filter (class != "Malasseziomycetes") %>%
	group_by (Treat, genus) %>% dplyr::summarize (count = n()) %>% summary()

taxids %>%inner_join(., partd) %>% filter (!is.na(otu)) %>%
	filter (class != "Malasseziomycetes") %>%
	group_by (Treat, otu) %>% dplyr::summarize (count = n())%>% summary()

# how many different orders, families, and genera

(taxids$order %>% unique() %>% length())-(taxids$order %>% unique() %>% grep("_sp|unclassified",.) %>% length())
(taxids$family %>% unique()%>% length())- (taxids$family %>% unique()%>% grep("_sp|unclassified",.) %>% length())
(taxids$genus %>% unique() %>% length())-(taxids$genus %>% unique()%>% grep("_sp|unclassified",.) %>% length())

# OTUs found in each soil amendment treatment

control.otus <- ((taxids %>%inner_join(., partd) )[(taxids %>%inner_join(., partd) )$Treat=='Control', 'otu'] %>% unique())$otu
plantation.otus <- ((taxids %>%inner_join(., partd) )[(taxids %>%inner_join(., partd) )$Treat=='Plantation', 'otu'] %>% unique())$otu
forest.otus <- ((taxids %>%inner_join(., partd) )[(taxids %>%inner_join(., partd) )$Treat=='Forest', 'otu'] %>% unique())$otu

# look at which OTUs occur in different treatments, are shared between them, or are unique

intersect(control.otus,plantation.otus)
intersect(forest.otus,plantation.otus)
intersect(control.otus,forest.otus)
intersect(forest.otus,plantation.otus) %>% intersect (control.otus)

library(VennDiagram)
venn.diagram(list(control.otus,plantation.otus,forest.otus),"VennNov30.2020.png",imagetype="png",category.names=c("Control","Plantation","Forest"), cat.pos=c(0,0,180))

print(otu.by.genus, n=80)

# more renaming and variable stuff

otu.by.family$family[grep("_sp|_fam_Incertae_sedis", otu.by.family$family, perl=T)] <- 'Other'

otu.by.family$family <- as.factor(otu.by.family$family)

otu.by.family$family <- factor(otu.by.family$family, levels = c(levels(otu.by.family$family)[which(levels(otu.by.family$family) !='Other')],'Other'))

# make a family stacked bar plot

colors.n.fam <- length(levels(otu.by.family$family))

p<-otu.by.family %>% left_join(.,
otu.by.family %>% group_by (Treat) %>% dplyr::summarize (t = sum(count))) %>%
mutate(relabund= count/t) %>% dplyr::rename(Family=family) %>%
	ggplot(aes(x=Treat, y=relabund, fill=Family)) +
		geom_bar(aes(), stat="identity", position="fill")+
		guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  		scale_fill_manual(values=colorRampPalette(brewer.pal(12, 'Set1'))(colors.n.fam)) +
  		xlab("") +
  		ylab("") +
		theme(axis.text.x=element_text(size=14, angle = 45, hjust = 1),
			axis.text.y=element_text(size = 14),
			axis.line.y=element_line(),
			axis.title=element_text(size=24),
			axis.line.x = element_blank(),
			axis.ticks.x = element_blank(),
			panel.border = element_blank(),
			panel.background = element_rect(fill = "transparent"),
    			plot.background = element_rect(fill = "transparent", color = NA),
    			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
    		legend.background = element_rect(fill = "transparent"),
    		legend.box.background = element_blank(),
    		legend.text = element_text(size=12),
    		legend.spacing.x = unit(0.5, 'cm'),
    		legend.spacing.y = unit(0.5, 'cm'),
    		legend.title = element_text(size=14))

#ggsave(p, width=10, height=7, filename="Family.bars.3.31.2020.png", bg = "transparent")

# genus stacked bar plot with colors from Colombia tanagers

otu.by.genus$genus <- as.factor(otu.by.genus$genus)

library(tanagR)

colors.n.genus <- length(levels(otu.by.genus$genus))

otu.by.genus$genus %>% unique()
levels(otu.by.genus$genus)[which(levels(otu.by.genus$genus)=="Neocosmospora")]<-"Fusarium Solani Species Complex"
otu.by.genus$genus %>% unique()

gplot.genus <- otu.by.genus %>% left_join(.,
otu.by.genus %>% group_by (Treat) %>% dplyr::summarize (t = sum(count))) %>%
mutate(relabund= count/t) %>% dplyr::rename(Genus=genus) %>%
	ggplot(aes(x=Treat, y=relabund, fill=Genus)) +
		geom_bar(aes(), stat="identity", position="fill")+
		guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  		xlab("") +
  		ylab("Relative abundance\n") +
  		scale_x_discrete(labels = c("Steam-treated","Plantation","Forest")) +
		theme(axis.text.x=element_text(size=14, angle = 45, hjust = 1),
			axis.text.y=element_text(size = 14),
			axis.line.y=element_line(),
			axis.title=element_blank(),#text(size=14),
			axis.line.x = element_blank(),
			axis.ticks.x = element_blank(),
			panel.border = element_blank(),
			panel.background = element_rect(fill = "transparent"),
    			plot.background = element_rect(fill = "transparent", color = NA),
    			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
    		legend.background = element_rect(fill = "transparent"),
    		legend.box.background = element_blank(),
    		legend.text = element_text(size=12),
    		legend.spacing.x = unit(0.5, 'cm'),
    		legend.spacing.y = unit(0.5, 'cm'),
    		legend.title = element_text(size=14))

# order stacked bar plot

gplot <- otu.by.order %>% left_join(.,
otu.by.order %>% group_by (Treat) %>% dplyr::summarize (t = sum(count)))%>%
mutate(relabund= count/t) %>% dplyr::rename(Order=order) %>%
	ggplot(aes(x=Treat, y=relabund, fill=Order)) +
		geom_bar(aes(), stat="identity", position="fill")+
		guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  		xlab("") +
  		ylab("Relative abundance\n") +
  		scale_x_discrete(labels = c("Steam-treated","Plantation","Forest")) +
		theme(axis.text.x=element_text(size=14, angle = 45, hjust = 1),
			axis.text.y=element_text(size = 14),
			axis.line.y=element_line(),
			axis.title=element_blank(),#text(size=14),
			axis.line.x = element_blank(),
			axis.ticks.x = element_blank(),
			panel.border = element_blank(),
			panel.background = element_rect(fill = "transparent"),
    			plot.background = element_rect(fill = "transparent", color = NA),
    			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
    		legend.background = element_rect(fill = "transparent"),
    		legend.box.background = element_blank(),
    		legend.text = element_text(size=12),
    		legend.spacing.x = unit(0.5, 'cm'),
    		legend.spacing.y = unit(0.5, 'cm'),
    		legend.title = element_text(size=14))

#gplot + scale_fill_manual(values=colorRampPalette(c(tanagr_palette("tangara_parzudakii"),tanagr_palette("chlorochrysa_nitidissima"),tanagr_palette("tangara_chilensis")))(colors.n.genus))

# genus stacked bar plot with nice colors

p1 <- gplot.genus + scale_fill_manual(values=colorRampPalette(c(tanagr_palette("tangara_seledon"),tanagr_palette("ramphocelus_sanguinolentus"),tanagr_palette("chlorochrysa_nitidissima"),tanagr_palette("tangara_chilensis")))(colors.n.genus))
#ggsave(p1, width=12, height=7, filename="Figure_S1.pdf", bg = "transparent", device="pdf")

Figure_3B <- p1

# order stacked bar plot with nice colors

custom_tanager <- c(tanagr_palette("tangara_seledon"),tanagr_palette("ramphocelus_sanguinolentus"),tanagr_palette("chlorochrysa_nitidissima"),tanagr_palette("tangara_chilensis"))[-c(13,7)]

length(custom_tanager)

p2 <- gplot + scale_fill_manual(values=custom_tanager)
#ggsave(p2, width=10, height=7, filename="Order.bars.11.30.2020.png", bg = "transparent")

### OTUs found in each treat

s.otu <- with(partd, otu[Treat == 'Control'])
p.otu <- with(partd, otu[Treat == 'Plantation'])
f.otu <- with(partd, otu[Treat == 'Forest'])

# plantation but not control
p.only <- setdiff(p.otu,s.otu)

# forest but not control
f.only <- setdiff(f.otu,s.otu)

p.only
f.only
setdiff(f.only,p.only)
setdiff(p.only,f.only)
intersect(f.only,p.only)

## data and plots with Gm Inoculation

## inoculation-soil amendment treatment interactive stacked barplot by order

joined_d <- taxids %>% inner_join(., partd) %>% plyr::join(., inoc.design[,c(3,6)], by='plant', type='left', match="first")

otu.by.order.inoc.treat <-joined_d %>% filter (!is.na(order)) %>% filter (!is.na(Inoc)) %>% filter (class != "Malasseziomycetes") %>% group_by (Treat, Inoc, order) %>% dplyr::summarize (count = n())

soil.labs <- c("Steam-treated","Plantation","Forest")
names(soil.labs) <- c("Control","Plantation","Forest")

interaction.stackedgplot <- otu.by.order.inoc.treat %>% left_join(., otu.by.order.inoc.treat %>% group_by(Treat,Inoc) %>% dplyr::summarize (t = sum(count))) %>% mutate(relabund= count/t)%>% dplyr::rename(Order=order) %>%
	ggplot(aes(x=Inoc, y=relabund, fill=Order)) +
		geom_bar(aes(), stat="identity", position="fill")+ facet_grid(~Treat, labeller=labeller(Treat = soil.labs))+
		guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  		xlab("") +
  		ylab("Relative abundance\n") +
  		scale_x_discrete(labels = c("Agar-only control", expression(paste(italic("G. morbida"))))) +
		theme(
			strip.text.x = element_text(size = 14, color="black", face="bold"),
			strip.background = element_blank(),
			axis.text.x=element_text(size=14, angle = 45, hjust = 1),
			axis.text.y=element_text(size = 14),
			axis.line.y=element_line(),
			axis.title=element_blank(),#text(size=20, face="bold"),
			axis.line.x = element_blank(),
			axis.ticks.x = element_blank(),
			panel.border = element_blank(),
			panel.background = element_rect(fill = "transparent"),
    			plot.background = element_rect(fill = "transparent", color = NA),
    			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
    		legend.background = element_rect(fill = "transparent"),
    		legend.box.background = element_blank(),
    		legend.text = element_text(size=12),
    		legend.spacing.x = unit(0.5, 'cm'),
    		legend.spacing.y = unit(0.5, 'cm'),
    		legend.title = element_text(size=14))+scale_fill_manual(values=custom_tanager)
#ggsave(interaction.stackedgplot , width=12, height=7, filename="Figure_3.pdf", bg = "transparent", device="pdf")

Figure_3A <- interaction.stackedgplot

t <- text_grob("Relative abundance\n", rot=90, size=16, family="Arial")
lay <- rbind(c(1,rep(2,11)), c(1,rep(3,11)))
pdf("Fig_3_composite.pdf", width=12, height=12)
grid.arrange(t, Figure_3A + labs(title="A")+theme(plot.title=element_text(size=14, face="bold", family="Arial")), Figure_3B + labs(title="B")+theme(plot.title=element_text(size=14, face="bold", family="Arial")), layout_matrix=lay)
dev.off()

##### ADONIS


s.otu <- with(partd, otu[Treat == 'Control'])
p.otu <- with(partd, otu[Treat == 'Plantation'])
f.otu <- with(partd, otu[Treat == 'Forest'])

# plantation but not control
p.only <- setdiff(p.otu,s.otu)

# forest but not control
f.only <- setdiff(f.otu,s.otu)

p.only
f.only
setdiff(f.only,p.only)
setdiff(p.only,f.only)
intersect(f.only,p.only)

##### ADONIS

# define singletons as ocurring more than 5 times

non.singletons2 <- taxids %>% inner_join(., partd) %>% filter (class != "Malasseziomycetes") %>% group_by (otu) %>% filter (n() > 5) %>% ungroup()

# need to add which ones had no isolates

xtra<-read.csv("Data/No.isolates.ForPlant.Nov30.csv")

joined_d <- taxids %>% inner_join(., partd) %>% plyr::join(., inoc.design[,c(3,6)], by='plant', type='left', match="first")

fullmatrix<-with(joined_d, table(plant,otu))

dim(fullmatrix)
dim(xtra)

# just singletons

nonsinglematrix <- fullmatrix[,colnames(fullmatrix) %in% non.singletons2$otu]

rownames(inoc.design) <- inoc.design$plant

communitymatrix<-nonsinglematrix [-which(rowSums(nonsinglematrix)==0),]

taxids[taxids$otu %in% colnames(communitymatrix),c('otu','family','genus','species')]

#check matrix
rowSums(communitymatrix)
dim(communitymatrix)

# prepare data for adonis
gm2<-inoc.design[rownames(communitymatrix),"Inoc"] 
soil2<-as.factor(inoc.design[rownames(communitymatrix),"Soil"])

# check matrix dimensions
length(gm2);length(soil2);dim(communitymatrix)
summary(gm2)
summary(soil2)

# distance matrix jaccard
d <- vegdist(as.matrix(communitymatrix), method='jaccard')

# conduct adonis with different orders of effects
adonis2(d ~ soil2:gm2 + soil2 + gm2, permutations=9999)
adonis2(d ~ soil2:gm2 + gm2 + soil2, permutations=9999)
adonis2(d ~ soil2 + soil2:gm2 + gm2, permutations=9999)
adonis2(d ~ gm2 + soil2:gm2 + soil2, permutations=9999)
adonis2(d ~ soil2 + soil2:gm2 + gm2, permutations=9999)
adonis2(d ~ gm2 + soil2:gm2 +  soil2, permutations=9999)

adonis2(d ~ gm2 +  soil2+ soil2:gm2 , permutations=9999)
adonis2(d ~ soil2+gm2 +  soil2:gm2 , permutations=9999)

# without interaction effect which was not significant
adonis2(d ~ gm2 +  soil2 , permutations=9999)


#### Poisson regression on Fusarium + Rhizoctonia spp.

# commented out for now, allows looking at different taxonomic levels/groupings

#f.j_d <- joined_d$genus %in% c("Neocosmospora","Fusarium","Gibberella","Nectria","Nectriaceae_sp")
#f.j_d <- joined_d$otu == "Otu001"
f.j_d <- joined_d$genus %in% c("Neocosmospora","Fusarium")
#f.j_d <- joined_d$family == "Nectriaceae"
#r.j_d <- joined_d$otu == "Otu002"
#r.j_d <- joined_d$family == "Ceratobasidiaceae"
r.j_d <- joined_d$genus == "Rhizoctonia"

# look at relabunds of diff Rhizoctonia and Fusarium
xtabs(~f.j_d + Inoc + Treat, joined_d)
xtabs(~r.j_d + Inoc + Treat, joined_d)

#create data for use in analysis
joined_d$Nectriaceae <- f.j_d
joined_d$Ceratobasidiaceae <- r.j_d

occurrence.data <- summaryBy(Ceratobasidiaceae + Nectriaceae ~ Treat + Inoc + plant, data=joined_d, FUN=sum , keep.names=T)
levels(occurrence.data$Treat)[1]<-"Sterile"

occurrence.data <- rbind(occurrence.data, xtra)
xtabs(~Inoc+Treat, occurrence.data)
dim(occurrence.data)

## RHIZOCTONIA

# Fusarium poisson model
pois.fus <- glm(Nectriaceae ~ Treat * Inoc, "poisson", occurrence.data)
Anova(pois.fus)

# Fusarium calculate dispersion
phihat.fus<-sum(residuals(pois.fus, type="pearson")^2) / pois.fus$df.residual
phihat.fus

# Fusarium quasipoisson
summ.qpf<-summary(pois.fus, dispersion=phihat.fus)
summ.qpf
qpois.fus <- glm(Nectriaceae ~ Treat * Inoc, data = occurrence.data, family = quasipoisson(link = "log"))
Anova(qpois.fus, Type=3)

# Fusarium group comparisons
pairs(emmeans(qpois.fus, ~ Treat * Inoc), adjust='none', type="response")
pairs(emmeans(qpois.fus, ~ Treat * Inoc), adjust='none', type="response", reverse="T")
fus.groups <- lsmeans(qpois.fus, ~ Treat * Inoc, type="response") %>% cld(Letters=letters, adjust='none')

# Fusarium look at R2
with(summ.qpf, 1-deviance/null.deviance)

### RHIZOCTONIA

# Rhizoctonia poisson model
pois.rhz <- glm(Ceratobasidiaceae ~ Treat * Inoc, "poisson", occurrence.data)
Anova(pois.rhz)

# Rhizoctonia calculate dispersion
phihat.rhz<-sum(residuals(pois.rhz, type="pearson")^2) / pois.rhz$df.residual
phihat.rhz

# Rhizoctonia quasipoisson
summ.qpr<-summary(pois.rhz, dispersion=phihat.rhz)
summ.qpr
qpois.rhz <- glm(Ceratobasidiaceae ~ Treat * Inoc, data = occurrence.data, family = quasipoisson(link = "log"))
Anova(qpois.rhz, Type=3)

# Rhizoctonia group comparisons
pairs(emmeans(qpois.rhz, ~ Treat * Inoc), adjust='none', type="response")
pairs(emmeans(qpois.rhz, ~ Treat * Inoc), adjust='none', type="response", reverse="T")
rhz.groups <- lsmeans(qpois.rhz, ~ Treat * Inoc, type="response") %>% cld(Letters=letters, adjust='none')
with(summ.qpr, 1-deviance/null.deviance)

# plots

# Fusarium plot

yl.fus<-expression(paste("Isolation rate of \n", italic("Fusarium\n"), " spp.\n"))
fusplot<-ggplot(data = fus.groups, aes(x = Treat, y = rate, fill=Inoc)) +
  geom_bar( stat="identity", color="black", position=position_dodge(), size=0.75) +
  geom_errorbar(aes(ymin=rate, ymax=rate-SE), position=position_dodge(0.9), width=0.2, color=c("black","black","black","black","black","white"))+
  geom_errorbar(aes(ymin=rate, ymax=rate+SE), position=position_dodge(0.9), width=0.2)+
theme_bw() +
  scale_x_discrete(labels = c("Steam-treated","Plantation","Forest")) +
  scale_fill_manual(values = c("black","grey"), labels=c("Agar-only control",expression(italic("G. morbida")))) + 
  theme(
  	plot.margin = margin(1,1,1,1,"cm"),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none") +
  ylab(yl.fus)+
  xlab("\nSoil Amendment")

# Rhizoctonia plot

yl.rhz<-expression(paste("Isolation rate of \n", italic("Rhizoctonia\n"), " spp.\n"))
rhzplot<-ggplot(data = rhz.groups, aes(x = Treat, y = rate, fill=Inoc)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), size=.75) +
  geom_errorbar(aes(ymin=rate, ymax=rate-SE), position=position_dodge(0.9), width=0.2, color=c("white","black","white","black","white","black"))+
  geom_errorbar(aes(ymin=rate, ymax=rate+SE), position=position_dodge(0.9), width=0.2)+
  theme_bw() +
  scale_x_discrete(labels = c("Steam-treated","Plantation","Forest")) +
  scale_fill_manual(values = c("black","grey"), labels=c("Agar-only control",expression(italic("G. morbida")))) + 
  theme(
  	plot.margin = margin(1,1,1,1,"cm"),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = c(.8,.9),
    legend.text.align = 0) +
  ylab(yl.rhz) +
  xlab("\nSoil Amendment") +
  labs(fill = "Inoculation")
#quartz("Figure_4.pdf", width=12, height=6)
#pdf("Figure_4.pdf", width=12, height=6)
#grid.arrange(fusplot + labs(title="(a)"), rhzplot + labs(title="(b)"), nrow=1)
#dev.off()

# Composite Fusarium-Rhizoctonia plot

pdf("/Users/will1809/OneDrive - purdue.edu/Dissertation/GMW.Dissertation.Analyses/Figures/CH2PhytobiomesJ/Figure_4.pdf", width=12, height=6)
grid.arrange(fusplot + labs(title="A")+theme(plot.title=element_text(size=14, face="bold", family="Arial")), rhzplot + labs(title="B")+theme(plot.title=element_text(size=14, face="bold", family="Arial")), nrow=1)
dev.off()