setwd("..")

fusarium<-read.csv("Data/Fus.Clust1.All.Tissues.3.31.2020.csv" )
cankers<-read.csv("Data/FT.Soil.Gm.2019.cankers.csv" )

library(doBy)
library(emmeans)
library(multcomp)

source("Code/Canker.functions.R")

fusarium
fusarium$Treatment <- factor(fusarium$Treatment, levels=c("NS","C","T","F","FT"))

max.soil <- fusarium[fusarium$Source == "Soil","Colonies"] %>% max()
fusarium.soiltrans <- fusarium
fusarium.soiltrans[fusarium$Source == "Soil","Colonies"] <- fusarium[fusarium$Source == "Soil","Colonies"] / max.soil

fusarium.summary <- summaryBy(Colonies ~ Treatment + Source, data=fusarium.soiltrans, FUN=c(sd,mean,length))

fusarium.summary$HSD<-NA

reiso.roots <- read.csv("Data/FT.Root.Reisolation.csv" )
reiso.roots$Treatment <- factor(reiso.roots$Treatment, levels=c("NS","C","T","F","FT"))
reiso.roots
pois.root<- glm(formula = Isolates ~ Treatment + Media, data = reiso.roots, family = poisson(link = "log"), weights=Root.pieces)
phihat.roots <- sum(residuals(pois.root, type="pearson")^2) / pois.root$df.residual
phihat.roots
root.summary<-summary(pois.root, dispersion=phihat.roots)
qpois.root <- glm(formula = Isolates ~ Treatment + Media, data = reiso.roots, family = quasipoisson(link = "log"), weights=Root.pieces)
Anova(qpois.root, Type=2)
root.groups <- lsmeans(qpois.root, ~ Treatment) %>% cld(Letters=letters, adjust='none')
rownames(root.groups)<-root.groups$Treatment
1-root.summary$deviance/root.summary$null.deviance

reiso.canks <- read.csv("Data/Fus.Canker.Reisolation.csv" )
reiso.canks$Treatment <- factor(reiso.canks$Treatment, levels=c("NS","C","T","F","FT"))
reiso.canks

pois.cank<-glm(formula = Fus ~ Treatment, data = reiso.canks, family = poisson(link = "log"), weights=Pieces)
phihat.canks <- sum(residuals(pois.cank, type="pearson")^2) / pois.cank$df.residual
phihat.canks
cank.summary<-summary(pois.cank, dispersion=phihat.canks)
qpois.cank <- glm(formula = Fus ~ Treatment, data = reiso.canks, family = quasipoisson(link = "log"), weights=Pieces)
Anova(qpois.cank, Type=2)
cank.groups <- lsmeans(qpois.cank, ~ Treatment) %>% cld(Letters=letters, adjust='none')
rownames(cank.groups)<-cank.groups$Treatment
1-cank.summary$deviance/cank.summary$null.deviance

pois.gm<-glm(formula = Gm ~ Treatment, data = reiso.canks, family = poisson(link = "log"), weights=Pieces)
phihat.gm.canks <- sum(residuals(pois.gm, type="pearson")^2) / pois.gm$df.residual
phihat.gm.canks
gm.summary <- summary(pois.gm, dispersion=phihat.gm.canks)
gm.summary
qpois.gm <- glm(formula = Gm ~ Treatment, data = reiso.canks, family = quasipoisson(link = "log"), weights=Pieces)
Anova(qpois.gm, Type=2)
gm.groups <- lsmeans(qpois.gm, ~ Treatment) %>% cld(Letters=letters, adjust='none')
rownames(gm.groups)<-gm.groups$Treatment
1-gm.summary$deviance/gm.summary$null.deviance


pois.soil<- glm(formula = Colonies ~ Treatment , data = fusarium[fusarium$Source=="Soil",], family = poisson(link = "log"))
phihat.soil <- sum(residuals(pois.soil, type="pearson")^2) / pois.soil$df.residual
phihat.soil
soil.summary<-summary(pois.soil, dispersion=phihat.soil)
qpois.soil <- glm(formula = Colonies ~ Treatment , data = fusarium[fusarium$Source=="Soil",], family = quasipoisson(link = "log"))
Anova(qpois.soil, Type=2)
soil.groups <- lsmeans(qpois.soil, ~ Treatment) %>% cld(Letters=letters, adjust='none')
rownames(soil.groups)<-soil.groups$Treatment
1-soil.summary$deviance/soil.summary$null.deviance

pairs(emmeans(qpois.root, ~ Treatment, type='response'), adjust='none') 
pairs(emmeans(qpois.root, ~ Treatment, type='response'), adjust='none', reverse="T") 
pairs(emmeans(qpois.cank, ~ Treatment, type='response'), adjust='none') 
pairs(emmeans(qpois.cank, ~ Treatment, type='response'), adjust='none', reverse="T") 
pairs(emmeans(qpois.soil, ~ Treatment, type='response'), adjust='none') 
pairs(emmeans(qpois.soil, ~ Treatment, type='response'), adjust='none', reverse="T") 

pairs(emmeans(qpois.gm, ~ Treatment, type='response'), adjust='none') 
pairs(emmeans(qpois.gm, ~ Treatment, type='response'), adjust='none', reverse="T") 



fusarium.summary[fusarium.summary$Source=="Roots",'HSD'] <- root.groups[c("NS","C","T","F","FT"),".group"] %>% trimws()
fusarium.summary[fusarium.summary$Source=="Canker",'HSD'] <- cank.groups[c("NS","C","T","F","FT"),".group"] %>% trimws()
fusarium.summary[fusarium.summary$Source=="Soil",'HSD'] <- soil.groups[c("NS","C","T","F","FT"),".group"] %>% trimws()

fusarium.summary

pdf("Figures/FTG.reiso.Final.bw.pdf", width=9, height=5)
ggplot(data = fusarium.summary, aes(x = Treatment, y = Colonies.mean, fill=Source)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=Colonies.mean, ymax=Colonies.mean+Colonies.sd/sqrt(Colonies.length)), position=position_dodge(0.9), width=0.2) +
  scale_y_continuous(sec.axis = sec_axis(~.*max.soil, name="Number of colonies from soil\n")) +
  geom_text(mapping = aes(x=Treatment, y=Colonies.mean+Colonies.sd/sqrt(Colonies.length)+0.05, label=HSD), position=position_dodge(0.9)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = c(0.5,.9))  +
    scale_fill_manual(values=c("black","grey","white"))+
  ylab("Reisolation rate from plant\n") +
  xlab("\nTreatment")
dev.off()

####
cank <- melt(id.vars=c("Plant","Treatment"), cankers[-which(cankers$Treatment=='-'),c(1:2,10,13,16)])
cank$Treatment <- factor(cank$Treatment, levels=c("NS","C","T","F","FT"))
names(cank)[4]<-"mm2"
names(cank)[3]<-"Canker"

## reisolation data

## join reisolation data to cankers

cank.fus.join.sfa<-left_join(cank[,c("Plant","Canker","mm2")], reiso[reiso$Media=="SFA",], by="Plant")
cank.fus.join.tsm<-left_join(cank[,c("Plant","Canker","mm2")], reiso[reiso$Media=="TSM",], by="Plant")

# plot reisolation data by treatment
boxplot(Percent ~ Treatment, data=reiso[reiso$Media=="SFA",])
boxplot(Percent ~ Treatment, data=reiso[reiso$Media=="TSM",])

# plot canker size by reisolation frequency
plot(cank.fus.join.sfa$Percent, cank.fus.join.sfa$mm2, col=cank.fus.join.sfa$Treatment)

# by treatment
with(cank.fus.join.sfa[cank.fus.join.sfa$Treatment == "F",], plot(Percent, mm2, col=Treatment))
with(cank.fus.join.sfa[cank.fus.join.sfa$Treatment == "FT",], plot(Percent, mm2, col=Treatment))
with(cank.fus.join.sfa[cank.fus.join.sfa$Treatment == "NS",], plot(Percent, mm2, col=Treatment))
with(cank.fus.join.sfa[cank.fus.join.sfa$Treatment == "T",], plot(Percent, mm2, col=Treatment))
with(cank.fus.join.sfa[cank.fus.join.sfa$Treatment == "C",], plot(Percent, mm2, col=Treatment))

