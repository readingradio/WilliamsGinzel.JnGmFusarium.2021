setwd("..")
the.data<- read.csv('Data/forest plant cankers raw out 10 28 18.csv')

library(plyr)
library(data.table)library(car)library(ggplot2)library(doBy)library(plotrix)library(tidyverse)
library(gridExtra)
library(asbio)

source("Code/Canker.functions.R")

with(the.data, t.test(A.necr[Inoc=="Control"]))
with(the.data, t.test(A.necr[Inoc=="Control"],A.necr[Inoc=="Gm"]))

# need to set negative canker obtained through subtracting area of cork borer to zero

to.zero <- the.data$A.necr <= 0

# divide into gm and control and test for significant effect of Gm inoculation
# need to set 0 canker area to some small value > 0 for box-cox
the.data$A.necr[to.zero] <- 0.000001model <- lm(A.necr ~ Treat * Inoc, the.data)boxCox(model)lambda<-with(boxCox(model, plotit=FALSE), x[which.max(y)])the.data$A.bc <- the.data$A.necr^lambdamodel.bc <- lm(A.bc ~ Treat * Inoc, the.data)summary(model.bc)Anova(model.bc, type=3)reduced.bc <- lm(A.bc ~ Treat + Inoc, the.data)summary(reduced.bc)Anova(reduced.bc, type=3)##### within treatmentsF.bc <- lm(A.necr^lambda ~ Inoc, the.data[the.data$Treat == "F",])summary(F.bc)
summary(lm(A.necr ~ Inoc, the.data[the.data$Treat == "F",]))Anova(F.bc, type=3)
pairs(lsmeans(F.bc, "Inoc", type="response"))
P.bc <- lm(A.necr^lambda ~ Inoc, the.data[the.data$Treat == "P",])summary(P.bc)
summary(lm(A.necr ~ Inoc, the.data[the.data$Treat == "P",]))Anova(P.bc, type=3)
pairs(lsmeans(P.bc, "Inoc", type="response"))S.bc <- lm(A.necr^lambda ~ Inoc, the.data[the.data$Treat == "S",])summary(S.bc)
summary(lm(A.necr ~ Inoc, the.data[the.data$Treat == "S",]))Anova(S.bc, type=3)
pairs(lsmeans(S.bc, "Inoc", type="response"))###the.data$A.necr[to.zero] <- 0
gm <- the.data[the.data$Inoc == 'Gm',]con <-the.data[the.data$Inoc == 'Control',]

gm$Treat <- factor(gm$Treat, levels=c("S","P","F"))

# same results for all the transformations
gm.model <- analyze.necrosis.bc.outliers(treat=gm$Treat, necr=gm$A.necr)
gm.model
gm.model$Regression %>% summary()
ss <- gm.model$Regression %>% Anova(Type=2)
ss
ss$'Sum Sq'[1]/(ss$'Sum Sq'[1]+ss$'Sum Sq'[2])

gm[gm.model$outliers,]

# untransformed group means and se after removal of outliers
untrans <- analyze.necrosis(treat=gm$Treat[-gm.model$outliers], necr=gm$A.necr[-gm.model$outliers])
untrans
untrans$HSD.p$treat['S-F', 1:3] / untrans$Tukey$means['S',1]
untrans$lc[untrans$lc[,'contrast']=='S - F', 2:3] / untrans$Tukey$means['S',1]


levels(gm$Treat)<-c("Control","Plantation","Forest")

p1 <- plot.necrosis(gm$Treat, gm$A.necr, units="cm\u00B2", xlab="Soil Amendment Treatment", bw="F")
p1.bw <- plot.necrosis(gm$Treat, gm$A.necr, units="cm\u00B2", xlab="Soil Amendment", tick.labels=c("Steam-treated ","Plantation ","Forest "))

## callus analysis

design <- read.csv('Data/forest plant smbiome treat design.csv')
design$Treat <- regmatches(design$Soil,regexpr('(?<=^)[FPS]', design$Soil, perl= TRUE))

design <- design[c('Treat', 'Rep', 'Treatment')]
names(design)[3] <- 'Inoc'

callus <- read.csv('Data/forest plant callus data.csv')
callus.full <- join(callus, design, by=c('Treat','Rep'), type = 'left', match = 'first')
callus.gm <- callus.full[callus.full$Inoc == 'Gm',]

## reorder callus score

callus.gm$Callus <- as.factor(callus.gm$Callus)
levels(callus.gm$Callus)<-c(rev(levels(callus.gm$Callus)))
callus.gm$Treat <- factor(callus.gm$Treat, levels=c("S","P","F"))
ratings<<-callus.gm$Callus
treat <<- callus.gm$Treat

polr<-analyze.callus.polr(callus.gm$Treat, callus.gm$Callus)

with(polr$contrasts, data.frame(contrast, estimate=exp(estimate), lower=exp(estimate+1.96*SE),upper=exp(estimate-1.96*SE)))

with(polr$contrasts.rev, data.frame(contrast, estimate=exp(estimate), lower=exp(estimate+1.96*SE),upper=exp(estimate-1.96*SE)))

p2<-plot.callus(callus.gm$Treat, callus.gm$Callus, y.axis=F, bw="F")
p2.bw<-plot.callus(callus.gm$Treat, callus.gm$Callus, y.axis=F)

pdf("Figures/ForstPlantGm.smbiome.multi.bw.Final.pdf", width=10, height=4)
grid.arrange(p1.bw + labs(title="A"), p2.bw + labs(title="B"), nrow=1)
dev.off()


pdf("Figures/ForstPlantGm.smbiome.multi.Final.pdf", width=12, height=6)
grid.arrange(p1, p2, nrow=1)
dev.off()
