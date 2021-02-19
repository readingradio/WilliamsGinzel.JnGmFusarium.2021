setwd("..")

# read in the data
cankers<-read.csv("Data/FT.Soil.Gm.2019.Final.cankers.csv" )

library(reshape2)
library(car)
library(lme4)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(agricolae)
library(gridExtra)
library(emmeans)
library(nlme)
library(asbio)
library(MuMIn)

source("Code/Canker.functions.randomcorrected2.R")

## cankers

cank <- melt(id.vars=c("Plant","Treatment"), cankers[,c(1:2,10,13,16)])
cank$Treatment <- factor(cank$Treatment, levels=c("NS","C","T","F","FT"))
names(cank)[4]<-"mm2"
names(cank)[3]<-"Canker"

cank$Fus <- cank$Treat %in% c("F","FT")
cank$Trich <- cank$Treat  %in% c("T","FT")
cank$Sterile <- cank$Treat == "NS"

# try an interactive model with Fusarium and trichoderma
lm(mm2 ~ Fus * Trich, data=cank[cank$Sterile == F,]) %>% Anova(type=3)

# look at untransformed model
untrans<-analyze.necrosis(cank$Treatment, cank$mm2)
untrans$Regression %>% Anova()
untrans$Regression %>% summary()

# effect sizes in untransformed model
untrans.effects <- function(w,x,y,z)
	print(w$lc[untrans$lc[,'contrast']==paste(x,'-', y), 2:3] / untrans$Tukey$means[z,1])

# look at group mean differences in untransformed data
untrans.effects(untrans, "C", "F", "C")
untrans.effects(untrans, "T", "F", "T")
untrans.effects(untrans, "NS", "C", "C")
untrans.effects(untrans, "NS", "T", "T")
untrans.effects(untrans, "T", "FT", "FT")

# set 0s to very small number
cank$mm2[cank$mm2 == 0] <- 0.00000001

# box cox transformed analysis without random effects
out<-analyze.necrosis.bc.outliers(cank$Treatment, cank$mm2)
out$lambda
out$outliers

# with random effects
out.random<-analyze.necrosis.bc.outliers(cank$Treatment, cank$mm2, cank$Plant)
out.random$lambda
out.random$outliers
cank[out.random$outliers,]

# calculate R squared different methods
D<-REMLcrit(out.random$Regression)
D.0<- REMLcrit(lmer(mm2^out$lambda ~ 1 + (1|Plant), data=cank[-out$outliers,]))

# deviance based
1-D/D.0
1-D/deviance(lm(mm2^out$lambda ~ 1, data=cank[-out.random$outliers,]))

# another calculation
r.squaredGLMM(out.random$Regression)

# Nagelkerke, McFadden, and one other R square and chi square calculations for model
nagelkerke(out.random$Regression, null=lm(mm2^out$lambda ~ 1, data=cank[-out.random$outliers,]))
nagelkerke(out.random$Regression, null=lmer(mm2^out$lambda ~ 1 + (1|Plant), data=cank[-out.random$outliers,]))

# analysis of deviance by comparing a nested model
anova(lmer(mm2^out$lambda ~ Treatment+(1|Plant), data=cank[-out$outliers,]), lmer(mm2^out$lambda ~ (1|Plant), data=cank[-out$outliers,]))

# linear mixed effect model
lmer(mm2^out$lambda ~ Treatment+(1|Plant), data=cank[-out$outliers,]) %>% Anova()

# use R squared from lmer function
summary(out$Regression)$r.squared

# what were the outliers
cank[out$outliers,]
out.random

# analysis of deviance (type II)
out.random$Regression %>% Anova(Type=2)

d<-attr(out.random$Regression, "frame")

# more ANOVAs
anova(lmer(n^out$lambda ~ treat + (1|random), data=d), lm(n^out$lambda ~ 1, data=d))
anova(lmer(n^out$lambda ~ treat + (1|random), data=d), lmer(n^out$lambda ~ 1+ (1|random), data=d))

# look at regression output for mixed model
out.random$Regression %>% summary()
out.random$Regression %>% emmeans("treat") %>% pairs(adjust="none")

# put special characters and italics in the axis labels
Axis.Labels.Italic<-list(	"Non-sterile  ",
						"Control  ",
						expression(paste(italic("Trichoderma"), " Rh-366  ")),
						"FSSC Rh-217  ",
						expression(paste(italic("Fus")," + ",italic("Trich"),"  ")))

Axis.Labels.Long<-c("Non-sterile  ","Control  ","Trichoderma (Rh-366)  ","FSSC (Rh-217)  ","FSSC + Trich  ")
Axis.Labels.Short<-c("NS  ","C  ","T  ","F  ","FT  ")

# plots necrosis

necrplot.bw.long.tuk<- plot.necrosis(cank$Treatment, cank$mm2, xlab="Seed Inoculation Treatment", adjust="tukey", bc=TRUE, tick.labels=Axis.Labels.Long, random=cank$Plant)

plot.necrosis(cank$Treatment, cank$mm2, xlab="Seed Inoculation Treatment", adjust="tukey", bc=TRUE, tick.labels=Axis.Labels.Italic, random=cank$Plant)

## callus

callus <- melt(id.vars=c("Plant","Treatment"), cankers[,c(1:2,11,14,17)])
callus$Treatment <- factor(cank$Treatment, levels=c("NS","C","T","F","FT"))
names(callus)[4]<-"Rating"
names(callus)[3]<-"Canker"

# fixed effects only proportional odds (ordinal) generalized linear model

wilcox <- with(callus, analyze.callus(treat=Treatment, ratings=Rating))
polr   <- with(callus, analyze.callus.polr(treat=Treatment, ratings=as.factor(Rating)), adj='none')

polr

library(ordinal)

callus$Plant<-as.factor(callus$Plant)
callus$Rating<-as.factor(callus$Rating)

# mixed model
polr.rand <-clmm(Rating ~ Treatment + (1|Plant), data=callus)
summary(polr.rand)
str(polr)
polr$clm$model
polr.rand
str(polr.rand)

# check deviance and likelihood
deviance(polr.rand)
logLik(polr.rand)

# create null and partial models for deviance based model analytics
polr.rand.null <- clmm(Rating ~ 1 + (1|Plant), data=callus)
polr.null <- clm(Rating ~ 1, data=callus)
polr.fixed<- clm(Rating ~ Treatment, data=callus)

# saturated model
satvar<-factor(1:dim(callus)[1])
polr.satur<- clm(Rating ~ satvar, data=callus)

# analyses of deviance and R square calculations
anova(polr.rand, polr.rand.null)
nagelkerke(polr.rand, null = polr.null)
nagelkerke(polr.rand, null = polr.rand.null)
nagelkerke(polr.fixed, null = polr.null)

# from scratch analyses of deviance and R square calculations
dev.rand<- -2*logLik(polr.rand)+2*logLik(polr.satur)
dev.fixed<--2*logLik(polr.fixed)+2*logLik(polr.satur)
dev.null<- -2*logLik(polr.null)+2*logLik(polr.satur)
dev.rand.null<--2*logLik(polr.rand.null)+2*logLik(polr.satur)
1-(dev.rand[1]/dev.null[1]) # fixed+randomeffects
1-(dev.fixed[1]/dev.null[1])# fixed effects only
dev.rand.null[1] - dev.rand[1]
pchisq( dev.rand.null[1] - dev.rand[1], 4, lower.tail=F)

str(polr.rand)
fitted.var<-var(polr.rand$fitted.values)
rand.var<-polr.rand$ST$Plant[1]^2
resid.var<-callus$Rating

# group comparisons

polr.rand.emm <- emmeans(polr.rand, "Treatment", type='response')
contrasts.callus.rand <- pairs(polr.rand.emm, adjust='none') %>% as.data.frame()
contrasts.callus.rand.rev <- pairs(polr.rand.emm, adjust='none', reverse="T") %>% as.data.frame()
contrasts.nonresp <- pairs(emmeans(polr.rand, "Treatment"), adjust='none')
comp <- cldList(p.value ~ contrast, data=contrasts.callus.rand, threshold=0.05)

with(contrasts.callus.rand, data.frame(contrast, estimate=exp(estimate), lower=exp(estimate+1.96*SE),upper=exp(estimate-1.96*SE)))
with(contrasts.callus.rand.rev, data.frame(contrast, estimate=exp(estimate), lower=exp(estimate+1.96*SE),upper=exp(estimate-1.96*SE)))

with(polr$contrasts, data.frame(contrast, estimate=exp(estimate), lower=exp(estimate+1.96*SE),upper=exp(estimate-1.96*SE)))

with(polr$contrasts.rev, data.frame(contrast, estimate=exp(estimate), lower=exp(estimate+1.96*SE),upper=exp(estimate-1.96*SE)))

# PLOTS OF CALLUS DATA

calplot<-with(callus, plot.callus(treat=Treatment, ratings=as.factor(Rating), y.axis=F, bw='F', letters=comp))

calplot.bw<-with(callus, plot.callus(treat=Treatment, ratings=as.factor(Rating), y.axis=F, letters=comp))

## multipanel

#pdf("Figures/FTGm.multi.Final.noTukey.pdf", width=12, height=6)
#grid.arrange(necrplot, calplot, nrow=1)
#dev.off()

## multipanel bw
pdf("Figures/CH2PhytobiomesJ/Figure_5_final.pdf", width=10, height=4)
grid.arrange(necrplot.bw.long.tuk + labs(title="A")+theme(plot.title=element_text(size=14, face="bold", family="Arial")), calplot.bw + labs(title="B")+theme(plot.title=element_text(size=14, face="bold", family="Arial")), nrow=1)
dev.off()

random.effects.cank <- lme(mm2 ~ Treatment, random = ~1 | Canker, data=cank, method='ML', contrasts=)
summary(random.effects.cank)
Anova(random.effects.cank, type=2)

## is canker position significant?

fixed.effects.cank <- lm(mm2 ~ Treatment + Canker, data=cank)
summary(fixed.effects.cank)
Anova(fixed.effects.cank, type=3)

fixed.effects.cank.interactive <- lm(mm2 ~ Treatment * Canker, data=cank)
summary(fixed.effects.cank.interactive)
Anova(fixed.effects.cank.interactive, type=3)