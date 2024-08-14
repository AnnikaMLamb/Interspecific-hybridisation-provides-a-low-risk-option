##R script for analyses of hybrid and purebred coral performance in the ocean
#Data obtained using methods described in manuscript -  Interspecific hybridisation provides a low-risk option for increasing genetic diversity of reef-building corals - methods
#Written by Annika Lamb

#Load packages and functions
library(lme4)
library(ggplot2)
library(multcomp)
library(ggrepel)
library(lsmeans)
library(nlme)
library(tidyr)
library(survival)
library(survminer)
library(brms)
library(rstanarm)
library(coxme)
library(loo)
library(ggforce)
library(tidyverse)


# Summary function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  require(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

###Fertilisation data analyses###
##Import data
Fert <- read.table('HybridFieldProject_FertilisationData.csv', header=T, sep=',')
Fert$Cross<-as.factor(Fert$Cross)
Fert$Mother<-as.factor(Fert$Mother)
Fert$SampleID<-as.factor(Fert$SampleID)
Fert <- na.omit(Fert)
##Summarise data
Fert_ByDam<-summarySE(Fert, measurevar="Percentage", groupvars=c("Mother","Cross"),conf.interval=0.95)
Fert_ByDam
Fert_ByCross<-summarySE(Fert, measurevar="Percentage", groupvars=c("Cross"),conf.interval=0.95)
Fert_ByCross
##Graphic
ggplot(Fert,aes(x = Cross, y = Percentage, col = Cross))+geom_boxplot()+scale_x_discrete(labels=expression(italic('A. florida'),'FS hybrid', 'SF hybrid', italic('A. sarmentosa')))+
  ylim(0, 100) + ylab("Percentage fertilised") + xlab("Offspring group")+
  scale_color_manual(name="Offspring group",labels = c(expression(italic("A. florida")), "FS hybrid", "SF hybrid",expression(italic("A. sarmentosa"))),  values = c("darkorchid3", "darkolivegreen3","tomato","steelblue2"))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(data=Fert_ByCross, (aes(x=factor(Cross), y=min(Percentage)*0.7, label=paste("n =",N))))
##Generalised linear mixed effects modelling
fert_glmm_mod<-glmer(Fertilised~Cross + (1|Mother), data=Fert, family=poisson())
summary(fert_glmm_mod)
summary(glht(fert_glmm_mod, mcp(Cross="Tukey")))

###Survival analyses###
##Kaplan-Meier estimation
Survivorship <- read.table('HybridFieldExperiment_Survival_KMDataset.csv', header=T, sep=',')
Survivorship<- Survivorship %>% drop_na(Survival)
Survivorship$Frame<-as.factor(Survivorship$Frame)
Survivorship$UniqueID<-as.factor(Survivorship$UniqueID)
Survivorship$TileID<-as.factor(Survivorship$TileID)
Survivorship$Date<-as.Date(Survivorship$Date, "%Y-%m-%d")
Survivorship$Days <- as.numeric(Survivorship$Days)
Survivorship$Time <- as.numeric(Survivorship$Time)
#Build model
surv_object <- Surv(time = Survivorship$Days, event = Survivorship$Survival)
surv_object 
fit1 <- survfit(surv_object ~ Cross, data = Survivorship)
summary(fit1)
#Graphic
labels=c("*A. florida*", "FS hybrid", "SF hybrid", "*A. sarmentosa*")
KmGraph<-ggsurvplot(fit1, data = Survivorship, pval = TRUE, conf.int = TRUE, conf.int.style = "ribbon",  conf.int.alpha = 0.3, ncensor.plot = F, risk.table = T,  xlab = "Days post-deployment",break.time.by = 50, legend.labs = labels, palette = c("darkorchid3", "darkolivegreen3","tomato","steelblue2"))
KmGraph
#Compare offspring groups
res <- pairwise_survdiff(Surv(Days, Survival) ~ Cross,
                         data = Survivorship)
res
##Bayesian Generalised Linear Mixed Effects Modelling
#Build dataset without time point 0
Data <- read.table('HybridFieldExperiment_Survival_BGLMMDataset.csv', header=T, sep=',')
selected<-c("2","4","6","10")
Data_filtered_T24610<-Data[Data$TimePoint %in% selected,]
Data_filtered_T24610$TimePoint <- as.factor(Data_filtered_T24610$TimePoint)
#Build and run models
survival.form1 <- bf(Survival ~ Cross * TimePoint + (1|UniqueID),
                     family=bernoulli)
survival.brms1 <- brm(survival.form1, data=Data_filtered_T24610,prior=set_prior('normal(0,3)'),
                      iter=5000, warmup=2000, thin=5,
                      refresh=0, chains=4)
survival.form2 <- bf(Survival ~ Cross * TimePoint + (1|Frame/TileID/UniqueID),
                     family=bernoulli)
survival.brms2 <- brm(survival.form2, data=Data_filtered_T24610,prior=set_prior('normal(0,2.5)'),
                      iter=5000, warmup=2000, thin=5,
                      refresh=0, chains=4, control = list(adapt_delta = 0.95))
survival.form3 <- bf(Survival ~ Cross * TimePoint + (1|TileID/UniqueID),
                     family=bernoulli)
survival.brms3 <- brm(survival.form3, data=Data_filtered_T24610,prior=set_prior('normal(0,2.5)'),
                      iter=5000, warmup=2000, thin=5,
                      refresh=0, chains=4, control = list(adapt_delta = 0.95))
survival.form4 <- bf(Survival ~ Cross * TimePoint + (1|Frame/UniqueID),
                     family=bernoulli)
survival.brms4 <- brm(survival.form4, data=Data_filtered_T24610,prior=set_prior('normal(0,2.5)'),
                      iter=5000, warmup=2000, thin=5,
                      refresh=0, chains=4, control = list(adapt_delta = 0.95))
#Compare models
(loo.1 = loo(survival.brms1))
plot(loo.1)
(loo.2 = loo(survival.brms2))
plot(loo.2)
(loo.3 = loo(survival.brms3))
plot(loo.3)
(loo.4 = loo(survival.brms4))
plot(loo.4)
Comp1v2 <- loo_compare(loo.1,loo.2)
print(Comp1v2, simplify = FALSE)
Comp1v3 <- loo_compare(loo.1,loo.3)
print(Comp1v3, simplify = FALSE)
Comp1v4 <- loo_compare(loo.1,loo.4)
print(Comp1v4, simplify = FALSE)
Comp2v3 <- loo_compare(loo.2,loo.3)
print(Comp2v3, simplify = FALSE)
Comp2v4 <- loo_compare(loo.2,loo.4)
print(Comp2v4, simplify = FALSE)
Comp3v4 <- loo_compare(loo.3,loo.4)
print(Comp3v4, simplify = FALSE)
#Summarise best-performing model
summary(survival.brms1)
plot(survival.brms1)
emm_s.t <- emmeans(survival.brms1, pairwise ~ Cross|TimePoint)
emm_s.t 

####Colour analyses####
##Import data
Data <- read.table('HybridFieldExperiment_ColourDataset.csv', header=T, sep=',')
Data_filtered<- Data %>% drop_na(RelativeGrey)
Data_filtered$RelativeGrey<- as.numeric(Data_filtered$RelativeGrey)
Data_filtered$TimePoint<- as.numeric(Data_filtered$TimePoint)
Data_filtered$Date<-as.Date(Data_filtered$Date, "%Y-%m-%d")
##Graphic
summaryRelativeGrey<-summarySE(Data_filtered, measurevar="RelativeGrey", groupvars=c("Date","Cross"),conf.interval=0.95)
RelativeGreyColourGraph=ggplot(summaryRelativeGrey, aes(y=RelativeGrey, x=Date, color=Cross))+
  geom_ribbon(aes(ymin=RelativeGrey-se, ymax=RelativeGrey+se),position=position_dodge(0),fill = "grey90", alpha=0.5)+
  geom_point(size=2.5)+theme_bw()+labs(x="Month", y="Mean relative grey value", color = "Cross")+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b-%y") +
  scale_color_manual(name="Offspring group",labels = c(expression(italic("A. florida")), "FS hybrid", "SF hybrid",expression(italic("A. sarmentosa"))),  values = c("darkorchid3", "darkolivegreen3","tomato","steelblue2"))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))
RelativeGreyColourGraph
RelativeGreyColourGraph + geom_label_repel(label=summaryRelativeGrey$N, direction=c("both"), box.padding = 0.1, label.padding = 0.15, point.padding = 3, label.size = 0.1, arrow = NULL,show.legend = FALSE)
###LME
Data_filtered$TimePoint <- as.factor(Data_filtered$TimePoint)
Data_filtered$Frame <- as.factor(Data_filtered$Frame)
Data_filtered$TileID <- as.factor(Data_filtered$TileID)
Data_filtered$UniqueID <- as.factor(Data_filtered$UniqueID)
##Build and run models
Colour_mod1<-lme(RelativeGrey~Cross*TimePoint, random= ~1|UniqueID, data=Data_filtered)
Colour_mod2<-lme(RelativeGrey~Cross*TimePoint, random= ~1|TileID/UniqueID, data=Data_filtered)
Colour_mod3<-lme(RelativeGrey~Cross*TimePoint, random= ~1|Frame/TileID/UniqueID, data=Data_filtered)
Colour_mod4<-lme(RelativeGrey~Cross*TimePoint, random= ~1|Frame/UniqueID, data=Data_filtered)
##Compare models
anova(Colour_mod1,Colour_mod2)
anova(Colour_mod1,Colour_mod3)
anova(Colour_mod1,Colour_mod4)
anova(Colour_mod2,Colour_mod3)
anova(Colour_mod2,Colour_mod4)
anova(Colour_mod3,Colour_mod4)
##Analyse changes within groups using best-performing model
emm_s.ct <- emmeans(Colour_mod2, consec ~ TimePoint|Cross,  adjust = "bonferroni")
emm_s.ct


####Growth####
##Import data
Data_filtered <- read.table('HybridFieldExperiment_GrowthDataset.csv', header=T, sep=',')
Data_filtered<- Data_filtered %>% drop_na(Size)
Data_filtered$Cross<-as.factor(Data_filtered$Cross)
Data_filtered$Date<-as.Date(Data_filtered$Date, "%Y-%m-%d")
##Graphic 
summarySize<-summarySE(Data_filtered, measurevar="Size", groupvars=c("Date","TimePoint","Cross"),conf.interval=0.95)
Data_filtered$zoom <- ifelse(Data_filtered$TimePoint < 4, TRUE, FALSE)
SizeGraph=ggplot(summarySize, aes(y=Size, x=Date, color=Cross))+
  geom_ribbon(aes(ymin=Size-se, ymax=Size+se),position=position_dodge(0),fill = "grey90", alpha=0.5)+
  geom_point(size=2.5)+theme_bw()+labs(x="Month", y=expression(Area ~ mm^2), color = "Cross")+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b-%y") +
  scale_color_manual(name="Offspring group",labels = c(expression(italic("A. florida")), "FS hybrid", "SF hybrid",expression(italic("A. sarmentosa"))),  values = c("darkorchid3", "darkolivegreen3","tomato","steelblue2"))+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))
SizeGraph
SizeGraph + geom_label_repel(label=summarySize$N, direction=c("both"), box.padding = 0.1, label.padding = 0.15, point.padding = 2.5, label.size = 0.1, arrow = NULL,show.legend = FALSE)
SizeGraph + facet_zoom(xlim = c(as.Date("2020-09-01"), as.Date("2021-01-29")), ylim = c(0,200), horizontal = FALSE, zoom.size=0.25, shrink = F) 
SizeGraph + facet_zoom(xlim = c(as.Date("2020-09-01"), as.Date("2021-01-29")), ylim = c(0,200), horizontal = FALSE, zoom.size=0.25, shrink = F) +
  geom_label_repel(label=summarySize$N, direction=c("both"), box.padding = 0.1, label.padding = 0.15, point.padding = 2.5, label.size = 0.1, arrow = NULL,show.legend = FALSE)
##LME
#Build models
Data_filtered$TimePoint<-as.factor(Data_filtered$TimePoint)
ctrl <- lmeControl(opt='optim');
Growth_mod1<-lme(Size~Cross*TimePoint, random= ~1|Frame/TileID/UniqueID, data=Data_filtered,control=ctrl, weights=varIdent(form=~TimePoint))
summary(Growth_mod1)
Growth_mod2<-lme(Size~Cross*TimePoint, random= ~1|UniqueID, data=Data_filtered,control=ctrl, weights=varIdent(form=~TimePoint))
summary(Growth_mod2)
Growth_mod3<-lme(Size~Cross*TimePoint, random= ~1|Frame/UniqueID, data=Data_filtered,control=ctrl, weights=varIdent(form=~TimePoint))
summary(Growth_mod3)
Growth_mod4<-lme(Size~Cross*TimePoint, random= ~1|TileID/UniqueID, data=Data_filtered,control=ctrl, weights=varIdent(form=~TimePoint))
summary(Growth_mod4)
#Compare models
anova(Growth_mod1,Growth_mod2)
anova(Growth_mod1,Growth_mod3)
anova(Growth_mod1,Growth_mod4)
anova(Growth_mod2,Growth_mod3)
anova(Growth_mod2,Growth_mod4)
anova(Growth_mod3,Growth_mod4)
#Summarise results
#Amongst offspring groups at each time point
emm_s.ct <- emmeans(Growth_mod4, pairwise ~ Cross |TimePoint,adjust = "none")
emm_s.ct
#Between consecutive time points for each offspring group
emm_s.tc <- emmeans(Growth_mod4, consec ~  TimePoint|Cross, adjust = "none")
emm_s.tc 
