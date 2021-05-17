setwd("..")

# read in the data
the.data<- read.csv('Data/forest plant cankers raw out 10 28 18.csv')

library(plyr)
library(data.table)library(car)library(ggplot2)library(doBy)library(plotrix)library(tidyverse)
library(gridExtra)
library(asbio)

# helper code for analysis of cankers and healing
source("Code/Canker.functions.randomcorrected2.R")

the.data %>% xtabs(~ Treat + Inoc,.)

# t test on necrosis between inoculated and non inoculated cankers
with(the.data, t.test(A.necr[Inoc=="Control"]))
with(the.data, t.test(A.necr[Inoc=="Control"],A.necr[Inoc=="Gm"]))

# need to set negative canker obtained through subtracting area of cork borer to zero

# set some negative canker areas to zero (due to subtraction of the area of the cork borer)
to.zero <- the.data$A.necr <= 0

# divide into gm and control and test for significant effect of Gm inoculation
# need to set 0 canker area to some small value > 0 for box-coxthe.data$A.necr[to.zero] <- 0.000001

# initial box cox lambda calcultatemodel <- lm(A.necr ~ Treat * Inoc, the.data)boxCox(model)lambda<-with(boxCox(model, plotit=FALSE), x[which.max(y)])

# lambda transformed data columnthe.data$A.bc <- the.data$A.necr^lambda

# treat by inoc model with interactionmodel.bc <- lm(A.bc ~ Treat * Inoc, the.data)summary(model.bc)Anova(model.bc, type=3)

# no interactionreduced.bc <- lm(A.bc ~ Treat + Inoc, the.data)summary(reduced.bc)Anova(reduced.bc, type=3)##### within treatmentsF.bc <- lm(A.necr^lambda ~ Inoc, the.data[the.data$Treat == "F",])summary(F.bc)
summary(lm(A.necr ~ Inoc, the.data[the.data$Treat == "F",]))Anova(F.bc, type=3)
pairs(lsmeans(F.bc, "Inoc", type="response"))
P.bc <- lm(A.necr^lambda ~ Inoc, the.data[the.data$Treat == "P",])summary(P.bc)
summary(lm(A.necr ~ Inoc, the.data[the.data$Treat == "P",]))Anova(P.bc, type=3)
pairs(lsmeans(P.bc, "Inoc", type="response"))S.bc <- lm(A.necr^lambda ~ Inoc, the.data[the.data$Treat == "S",])summary(S.bc)
summary(lm(A.necr ~ Inoc, the.data[the.data$Treat == "S",]))Anova(S.bc, type=3)
pairs(lsmeans(S.bc, "Inoc", type="response"))## separate gm from control data
gm <- the.data[the.data$Inoc == 'Gm',]con <-the.data[the.data$Inoc == 'Control',]

# rename treatments
gm$Treat <- factor(gm$Treat, levels=c("S","P","F"))

# same results for all the transformations
gm.model <- analyze.necrosis.bc.outliers(treat=gm$Treat, necr=gm$A.necr)
gm.model

gm$Treat %>% table
gm$Treat %>% length

gm.model$Regression %>% summary()

$ R squared
ss <- gm.model$Regression %>% Anova(Type=2)
ss
ss$'Sum Sq'[1]/(ss$'Sum Sq'[1]+ss$'Sum Sq'[2])

# what are the outliers
gm[gm.model$outliers,]

# what are the outliers
gm[-gm.model$outliers,] %>% xtabs(~Treat+Inoc,.)

# untransformed group means and se after removal of outliers
untrans <- analyze.necrosis(treat=gm$Treat[-gm.model$outliers], necr=gm$A.necr[-gm.model$outliers])
untrans
untrans$HSD.p$treat['S-F', 1:3] / untrans$Tukey$means['S',1]
untrans$lc[untrans$lc[,'contrast']=='S - F', 2:3] / untrans$Tukey$means['S',1]

# rename treatment codes
levels(gm$Treat)<-c("Control","Plantation","Forest")

# plot
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

data.frame(treat=callus.gm$Treat, ratings=callus.gm$Callus) %>% xtabs(~ ratings+ treat,.)

## reorder callus score - error corrected May 14 2021 (see below)

callus.gm$Callus <- as.factor(callus.gm$Callus)
levels(callus.gm$Callus)<-c(rev(levels(callus.gm$Callus)))
callus.gm$Treat <- factor(callus.gm$Treat, levels=c("F","P","S")) ## was "levels=c("S","P","F")"
ratings<<-callus.gm$Callus
treat <<- callus.gm$Treat

data.frame(treat=callus.gm$Treat,ratings) %>% xtabs(~ ratings+ treat,.)

polr<-analyze.callus.polr(callus.gm$Treat, callus.gm$Callus)

# look at contrasts 95% confidence intervals
with(polr$contrasts, data.frame(contrast, estimate=exp(estimate), lower=exp(estimate+1.96*SE),upper=exp(estimate-1.96*SE)))
with(polr$contrasts.rev, data.frame(contrast, estimate=exp(estimate), lower=exp(estimate+1.96*SE),upper=exp(estimate-1.96*SE)))

# callus plots

p2<-plot.callus(callus.gm$Treat, callus.gm$Callus, y.axis=F, bw="F")
p2.bw<-plot.callus(callus.gm$Treat, callus.gm$Callus, y.axis=F, scale_labs=c("0: No healing", "1: Reactive margin", "2: Partly healed", "3: Fully healed"))

# composite plots
pdf("Figure_2_final_correct.pdf", width=11, height=4)
#pdf("Figures/CH2PhytobiomesJ/Figure_2_final_correct.pdf", width=10, height=4)
grid.arrange(p1.bw + labs(title="A")+theme(plot.title=element_text(size=14, face="bold", family="Arial")), p2.bw + labs(title="B")+theme(plot.title=element_text(size=14, face="bold", family="Arial")), nrow=1)
dev.off()

pdf("Figures/ForstPlantGm.smbiome.multi.Final_correct.pdf", width=12, height=6)
grid.arrange(p1, p2, nrow=1)
dev.off()
