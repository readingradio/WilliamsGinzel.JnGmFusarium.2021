library(tidyverse)
library(RColorBrewer)
library(vegan)
library(rcompanion)
library(car)
library(agricolae)
library(extrafont)
library(lme4)
library(ordinal)
library(effects)
library(nnet)
library(MASS)
library(MuMIn)
library(emmeans)
loadfonts()

plot.callus <- function(treat, ratings, y.axis=TRUE, test=c("polr","wilcox"), adj='none', bw='T', letters=NULL) {
	
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
	else if (is.null(letters)) le <- analyze.callus.polr(treat, ratings, adj)$groups
	else le <- letters

	caltab <- xtabs(~ Rating + TREAT, data.frame(TREAT=treat, Rating=ratings))
	reltab <- decostand(caltab, 'total', 2)
	toplot <- melt(reltab)
	toplot
	toplot$TREAT <- as.factor(toplot$TREAT)
	toplot$Rating <- as.factor(toplot$Rating)
	
	if (bw) base_plot <- toplot %>% ggplot() +
		geom_bar(aes(x=TREAT, y=value, fill=Rating), stat="identity", position="stack")+
		geom_text(data = le, aes(x=Group, y=1.1, label=Letter))+
		scale_fill_grey("Canker healing")#, direction=-1)
	else base_plot <- toplot %>% ggplot() +
		geom_bar(aes(x=TREAT, y=value, fill=Rating), stat="identity", position="stack")+
		geom_text(data = le, aes(x=Group, y=1.1, label=Letter))+
		scale_fill_brewer("Canker healing", palette='YlOrBr', direction=-1)

	 base_plot +
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

plot.necrosis <- function(treat, necr, trans.test = function (x) x, units="mm\u00B2", y.axis=TRUE, xlab="Treatment", remove.0s = FALSE, remove.0s.test = FALSE, bc=TRUE, adjust=c("Tukey", "none"), bw='T', tick.labels = levels(treat), random=NULL, cov=NULL, cook.denom=4) {


	#bc TRUE is if there was a box cox test
	#cook.denom is the denominator for outlier removal

	if (bc) {
		if (is.null(cov)) analysis <- analyze.necrosis.bc.outliers(treat=treat, necr=necr, random=random, cook.denom)
		else analysis <- analyze.necrosis.bcout.cov(treat=treat, necr=necr, cov=cov, random=random, cook.denom)
		#if (!(is.null(cook.denom) | length(analysis$outliers == 0))) {
			treat <- treat[-analysis$outliers]
			necr <- necr[-analysis$outliers]
		#}
	} else {
		if (is.null(cov)) analysis <- analyze.necrosis(treat=treat, necr=necr, trans=trans.test, remove.0s=(remove.0s.test), random=random)
		else analysis <- analyze.necrosis.cov(treat=treat, necr=necr, cov=cov, trans=trans.test, remove.0s=(remove.0s.test), random=random) 
	}

#pointstyle <- vector(length=length(necr))
#pointstyle <- rep(0,length(necr))
#if (length(analysis$outliers > 0)) pointstyle[analysis$outliers] <- 1
#pointstyle <- as.factor(pointstyle)

	if (remove.0s) {
		treat <- treat[necr > 0]
		necr <- necr[necr > 0]
	}
	
	print(analysis)
	
	spacer <- vector('numeric', length=2)
	names(spacer) <- c("mm\u00B2", "cm\u00B2")
	spacer[]<-c(2, 0.02)
	#print(spacer)
	
	bp <- boxplot(necr ~ treat, plot=F)
	#print(bp)
	#print(analysis$lcgroups$Letter)

	if (adjust == "none") {
		if(!is.null(analysis$lcgroups)) n.bp <- with(bp, data.frame(Treatment=names, pos=stats[5,], letter.group = analysis$lcgroups$Letter))
		else n.bp<-NULL
	} else {
		if (is.null(random)) n.bp <- with(bp, data.frame(Treatment=names, pos=stats[5,], letter.group = analysis$Tukey$groups[names,'groups']))
		else {
			n.bp <- with(bp, data.frame(Treatment=names, pos=stats[5,], letter.group = analysis$lctukgroups$Letter))
	}}

	#print(n.bp)

	if (bw) base_plot <- data.frame(n=necr, treat=treat) %>% ggplot(., mapping=aes(y=necr, x=treat)) +
		geom_boxplot(outlier.shape=NA) 
	else base_plot <- data.frame(n=necr, treat=treat) %>% ggplot(., mapping=aes(y=necr, x=treat)) +
		geom_boxplot(aes(fill = treat), outlier.shape=NA) 

	style<-theme_bw() +
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
			axis.text = element_text(size = 14))

	#base_plot+style

	if (is.null(n.bp)) base_plot +
		scale_x_discrete(labels = tick.labels) +
		coord_flip() +
		geom_jitter(width=0.2)+#, aes(shape=pointstyle)) +
		#scale_shape_manual(values=c(16,1)) +
		style+
			ylab(paste("\nNecrotic Area (", units, ")", sep="")) +
			xlab(paste(xlab, "\n", sep=""))
		
	else base_plot +
		geom_text(data=n.bp, aes(x=Treatment, y=pos + spacer[units], label=letter.group), vjust=0) +
		scale_x_discrete(labels = tick.labels) +
		coord_flip() +
		geom_jitter(width=0.2)+#, aes(shape=pointstyle)) +
		#scale_shape_manual(values=c(16,1)) +
		style+
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
	mod1.emm <- emmeans(mod1, "tr", type='response')
	contrasts.callus <- pairs(mod1.emm, adjust=adj) %>% as.data.frame()
	rev <- pairs(mod1.emm, adjust=adj, reverse="T") %>% as.data.frame()
	contrasts.nonresp <- pairs(emmeans(mod1, "tr"), adjust=adj)
	comp <- cldList(p.value ~ contrast, data=contrasts.callus, threshold=0.05)

	ret <- list(polr=summary(mod1), clm=summary(mod2), anova=Anova(mod1), R2=R2, contrasts=contrasts.callus, contrasts.rev=rev, logodds = contrasts.nonresp,groups=comp)
	
	# remove d from global environment
	remove(d,envir=.GlobalEnv)
	
	ret
}

wilcox.groups <- function (pw) {
	comp <- na.omit(melt(pw$p.value))
	comp$c <- with(comp, paste (Var1, Var2,  sep =" - "))
	cldList(value ~ c, data=comp, threshold=0.05)
}

analyze.necrosis <- function(treat, necr, trans = function (x) x, remove.0s=FALSE, random=NULL) {
	if (remove.0s) {
		treat <- treat[necr > 0]
		necr <- necr[necr > 0]
	}

	n <- trans(necr)

	if (is.null(random)) model <- lm(n ~ treat)

	else model <- lmer(n ~ treat + (1|treat:random))

	t3 <- Anova(model, type=3)

	lcpairs<-as.data.frame(pairs(emmeans(model,'treat'),adjust='none'))
	lctuk <- as.data.frame(pairs(emmeans(model,'treat')))
	#print(lcpairs)
	#print(data.frame(p=as.numeric(lcpairs$p.value),c=as.character(lcpairs$contrast)))
	comp <- tryCatch({with(lcpairs,cldList(p ~ c, data=data.frame(p=as.numeric(p.value),c=as.character(contrast)), threshold=0.05))}, error = function (err) {return(NULL)})
	
	comptuk<- tryCatch({cldList(p.value ~ contrast, data=lctuk, threshold=0.05)}, error = function (err) {return(NULL)})
	if (is.null(random)) {
		anv <- aov(n ~ treat)
		tuk <- TukeyHSD(anv, ordered=TRUE)
		hsd <- HSD.test(anv, "treat")
		return(list(Regression=model, ANOVA=t3, Tukey=hsd, HSD.p=tuk, lc=lcpairs, lcgroups=comp, lctk=lctuk, lctukgroups=comptuk))
	} else {
		return(list(Regression=model, ANOVA=t3, lc=lcpairs, lcgroups=comp, lctk=lctuk, lctukgroups=comptuk))
	}
}

analyze.necrosis.cov <- function(treat, necr, cov, trans = function (x) x, remove.0s=FALSE, random=NULL) {
	if (remove.0s) {
		treat <- treat[necr > 0]
		necr <- necr[necr > 0]
	}

	n <- trans(necr)

	if (is.null(random)) model <- lm(n ~ treat + cov)

	else model <- lmer(n ~ treat + cov + (1|treat:random))

	t3 <- Anova(model, type=3)

	lcpairs<-as.data.frame(pairs(emmeans(model,'treat'),adjust='none'))
	lctuk <- as.data.frame(pairs(emmeans(model,'treat')))
	#print(lcpairs)
	#print(data.frame(p=as.numeric(lcpairs$p.value),c=as.character(lcpairs$contrast)))
	
	comp <- tryCatch({with(lcpairs,cldList(p ~ c, data=data.frame(p=as.numeric(p.value),c=as.character(contrast)), threshold=0.05))}, error = function (err) {return(NULL)})
	
	#print(with(lcpairs,cldList(p ~ c, data=data.frame(p=as.numeric(p.value),c=as.character(contrast)), threshold=0.05)))
	
	comptuk<- tryCatch({cldList(p.value ~ contrast, data=lctuk, threshold=0.05)}, error = function (err) {return(NULL)})
	if (is.null(random)) {
		anv <- aov(n ~ treat)
		tuk <- TukeyHSD(anv, ordered=TRUE)
		hsd <- HSD.test(anv, "treat")
		return(list(Regression=model, ANOVA=t3, Tukey=hsd, HSD.p=tuk, lc=lcpairs, lcgroups=comp, lctk=lctuk, lctukgroups=comptuk))
	} else {
		return(list(Regression=model, ANOVA=t3, lc=lcpairs, lcgroups=comp, lctk=lctuk, lctukgroups=comptuk))
	}
}

analyze.necrosis.bc.outliers <- function(treat, necr, random=NULL, cook.denom=4) {
	
	# store data in dataframe to global environment because 
	# boxCox is poorly written and relies heavily on global
	# environment, causing an errror; see the github page
	# https://stackoverflow.com/questions/31921425/apply-and-boxcox-function-error-in-r
	d <<- data.frame(ne=necr, tr=treat)

	# start out by fitting model of gm-inoculated
	model.genv <- lm(ne ~ tr, data=d)

	# refit with box cox transformation
	lambda<- with(boxCox(model.genv, plotit=FALSE), x[which.max(y)])
	model2<- lm(ne^lambda ~ tr, data=d)

	# identify outliers
	if (is.null(cook.denom)) outliers <- rep(FALSE, length=length(necr))
	else outliers <- cooks.distance(model2) > cook.denom / length(necr)

	# remove d from global environment
	remove(d,envir=.GlobalEnv)

	# refit model and refit with boxcox transformation
	# and perform Tukey test with analyze.necrosis

	c(analyze.necrosis(treat[!outliers], necr[!outliers], function(x) x^lambda, random=random[!outliers]), list(outliers=which(as.vector(outliers))), lambda=lambda)
}


analyze.necrosis.bcout.cov <- function(treat, necr, cov, random=NULL, cook.denom=4) {
	
	# store data in dataframe to global environment because 
	# boxCox is poorly written and relies heavily on global
	# environment, causing an errror; see the github page
	# https://stackoverflow.com/questions/31921425/apply-and-boxcox-function-error-in-r
	d <<- data.frame(ne=necr, tr=treat, c=cov)

	# start out by fitting model of gm-inoculated
	model.genv <- lm(ne ~ tr + c, data=d)

	# refit with box cox transformation
	lambda<- with(boxCox(model.genv, plotit=FALSE), x[which.max(y)])
	model2<- lm(ne^lambda ~ tr +c, data=d)

	# identify outliers
	if (is.null(cook.denom)) outliers <- rep(FALSE, length=length(necr))
	else outliers <- cooks.distance(model2) > cook.denom / length(necr)

	# refit model and refit with boxcox transformation
	# and perform Tukey test with analyze.necrosis

	# remove d from global environment
	remove(d,envir=.GlobalEnv)

	c(analyze.necrosis.cov(treat[!outliers], necr[!outliers], cov[!outliers], function(x) x^lambda, random=random[!outliers]), list(outliers=which(as.vector(outliers))), lambda=lambda)
}


#analyze.sporulation <- function() {}