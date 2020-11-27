library(tidyverse)
library(RColorBrewer)
library(vegan)
library(rcompanion)
library(car)
library(agricolae)
library(extrafont)

library(ordinal)
library(effects)
library(nnet)
library(MASS)
library(MuMIn)
library(emmeans)
loadfonts()

plot.callus <- function(treat, ratings, y.axis=TRUE, test=c("polr","wilcox")) {
	
	if (y.axis) {
		axis.line.y = element_line(colour = "black")
		axis.title.y = element_text(size = 16, face = "bold")
		axis.text.y = element_text(size = 14)
		axis.ticks.y = element_line()
	} else {
		axis.line.y = element_blank()
		axis.title.y = element_blank()
		axis.text.y = element_blank()
		axis.ticks.y = element_blank()
	}
	
	if (test == "wilcox") le <- wilcox.groups(analyze.callus(treat, ratings)$Pairwise.Wilcox)
	else le <- analyze.callus.polr(treat, ratings)$groups

	caltab <- xtabs(~ Rating + TREAT, data.frame(TREAT=treat, Rating=ratings))
	reltab <- decostand(caltab, 'total', 2)
	toplot <- melt(reltab)
	toplot
	toplot$TREAT <- as.factor(toplot$TREAT)
	toplot$Rating <- as.factor(toplot$Rating)
	toplot %>% ggplot() +
		geom_bar(aes(x=TREAT, y=value, fill=Rating), stat="identity", position="stack")+
		geom_text(data = le, aes(x=Group, y=1.1, label=Letter))+
		scale_fill_brewer("Canker healing", palette='YlOrBr', direction=-1)+
		labs(title=NULL,
			x="Seed Inoculation Treatment\n",
			y="\nProportion of replicates") +
		theme(
			#text = element_text(family="Montserrat"),
			panel.border = element_blank(),
			panel.background = element_rect(fill = "transparent"), # bg of the panel
    			plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
			panel.grid = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.line.x = element_line(colour = "black"),
			axis.title.x = element_text(size = 16),
			axis.text.x = element_text(size = 14),
			axis.line.y = axis.line.y,
			axis.title.y = axis.title.y,
			axis.text.y = axis.text.y,
			axis.ticks.y = axis.ticks.y,
   	 		legend.background = element_rect(fill = "transparent"),
   	 		legend.box.background = element_blank(),
   	 		legend.text = element_text(size=14),
   	 		legend.title = element_text(size=18))+
		coord_flip()
}

plot.necrosis <- function(treat, necr, trans.test = function (x) x, units="mm\u00B2", y.axis=TRUE, xlab="Treatment", remove.0s = FALSE, remove.0s.test = FALSE, bc=TRUE, adjust=c("Tukey", "none")) {

	if (bc) {
		analysis <- analyze.necrosis.bc.outliers(treat=treat, necr=necr)
		treat <- treat[-analysis$outliers]
		necr <- necr[-analysis$outliers]
	}
	else {
		analysis <- analyze.necrosis(treat=treat, necr=necr, trans=trans.test, remove.0s=(remove.0s.test)) 
	}

	if (remove.0s) {
		treat <- treat[necr > 0]
		necr <- necr[necr > 0]
	}
	
	spacer <- vector('numeric', length=2)
	names(spacer) <- c("mm\u00B2", "cm\u00B2")
	spacer[]<-c(2, 0.02)
	#print(spacer)
	
	bp <- boxplot(necr ~ treat, plot=F)

	n.bp <- with(bp,
		data.frame(Treatment=names, pos=stats[5,], hsd.group = analysis$Tukey$groups[names,'groups']))
	n.bp <- left_join(n.bp, data.frame(Treatment=analysis$lcgroups$Group, LCGroup = analysis$lcgroups$Letter))
	if (adjust == "none") n.bp$hsd.group <- n.bp$LCGroup

	data.frame(n=necr, treat=treat) %>% ggplot(., mapping=aes(y=necr, x=treat)) +
		geom_boxplot(aes(fill = treat), outlier.shape=NA) +
		geom_text(data=n.bp, aes(x=Treatment, y=pos + spacer[units], label=hsd.group), vjust=0) +
		coord_flip() +
		geom_jitter(width=0.2) +
		theme_bw() +
		theme(
			#text = element_text(family="Montserrat"),
			panel.border = element_blank(),
			panel.background = element_rect(fill = "transparent"), # bg of the panel
    			plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      	  	legend.position = "none",
			panel.grid = element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title = element_text(size = 16),
			axis.text = element_text(size = 14)) +
			ylab(paste("\nNecrotic Area (", units, ")", sep="")) +
			xlab(paste(xlab, "\n", sep=""))

}

analyze.callus <- function(treat, ratings) {
	krusk <- kruskal.test(ratings ~ as.factor(treat))
	wilxpair <- pairwise.wilcox.test(ratings, treat, paired = FALSE, p.adjust.method = 'BH')
	list(Kruskal.Wallis = krusk, Pairwise.Wilcox = wilxpair)
}

analyze.callus.polr <- function(treat, ratings, adj="none") {
	
	d <<- data.frame(rg=ratings,tr=treat)
	
	mod1 <- polr(rg ~ tr, data=d, method="logistic")
	null1<- polr(rg ~ 1, data=d, method="logistic")
	
	mod2 <- clm(rg ~ tr, data=d, link="logit")
	
	R2 <- 1-(mod1$deviance / null1$deviance)
	mod1.emm <- emmeans(mod1, "tr")
	contrasts.callus <- pairs(mod1.emm, adjust=adj) %>% as.data.frame()
	comp <- cldList(p.value ~ contrast, data=contrasts.callus, threshold=0.05)

	ret <- list(polr=summary(mod1), clm=summary(mod2), anova=Anova(mod1), R2=R2, contrasts=contrasts.callus, groups=comp)
	
	# remove d from global environment
	remove(d,envir=.GlobalEnv)
	
	ret
}

wilcox.groups <- function (pw) {
	comp <- na.omit(melt(pw$p.value))
	comp$c <- with(comp, paste (Var1, Var2,  sep =" - "))
	cldList(value ~ c, data=comp, threshold=0.05)
}

analyze.necrosis <- function(treat, necr, trans = function (x) x, back = function(y) y, remove.0s=FALSE) {
	if (remove.0s) {
		treat <- treat[necr > 0]
		necr <- necr[necr > 0]
	}
	n <- trans(necr)
	model <- lm(n ~ treat)
	anv <- aov(n ~ treat)
	t3 <- Anova(model, type=3)
	tuk <- TukeyHSD(anv, ordered=TRUE)
	hsd <- HSD.test(anv, "treat")
	tuk$means[,1:2] <- back(tuk$means[,1:2]); tuk$groups[,1] <- back(tuk$groups[,1])#; hsd$treat[,1:3] <- back(tuk$treat[,1:3])
	lcpairs<-as.data.frame(pairs(emmeans(model,'treat'),adjust='none'))
	lcpairs[,2] <- back(lcpairs[,2])
	lcpairs[,3] <- back(lcpairs[,2]) - back(lcpairs[,2]-lcpairs[,3])
	lcpairs[,7] <- back(lcpairs[,2]+lcpairs[,3]) - back(lcpairs[,2])
	colnames(lcpairs)[3] <-"SE.low"
	colnames(lcpairs)[7] <- "SE.high"
	lcpairs<-lcpairs[,c(1:3,7,4:6)]
	comp <- cldList(p.value ~ contrast, data=lcpairs, threshold=0.05)
	list(Regression=model, ANOVA=t3, Tukey=hsd, HSD.p=tuk, lc=lcpairs, lcgroups=comp)
}

analyze.necrosis.bc.outliers <- function(treat, necr, cook.denom=4) {
	
	# store data in dataframe to global environment because 
	# boxCox is poorly written and relies heavily on global
	# environment, causing an errror; see the github page
	# https://stackoverflow.com/questions/31921425/apply-and-boxcox-function-error-in-r
	d <<- data.frame(ne=necr,tr=treat)

	# start out by fitting model of gm-inoculated
	model.genv <- lm(ne ~ tr, data=d)

	# refit with box cox transformation
	lambda<- with(boxCox(model.genv, plotit=FALSE), x[which.max(y)])
	model2<- lm(ne^lambda ~ tr, data=d)

	# identify outliers
	outliers <- cooks.distance(model2) > cook.denom / length(necr)

	# refit model and refit with boxcox transformation
	# and perform Tukey test with analyze.necrosis
	
	# remove d from global environment
	remove(d,envir=.GlobalEnv)
	
	c(analyze.necrosis(treat[-outliers], necr[-outliers], function(x) x^lambda, function(y) y^(1/lambda)), list(outliers=which(as.vector(outliers))), lambda=lambda)
}

#analyze.sporulation <- function() {}