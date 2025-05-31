
##################################### meta-analysis
######################## Map
library(dplyr)
library(ggplot2)
library(maps)
library(tidyverse)
library(readxl)
world.dat<-map_data("world")
ggplot() +
  geom_polygon(data=world.dat,aes(x=long,y=lat,group=group),
               fill="#dedede")+
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(x=NULL,y=NULL)-> world.map
world.map
### Map
df<-read.csv(file.choose())
world.map + 
  geom_point(data = df, 
             aes(x = Longitude, 
                 y = Latitude, 
                 color = Type, 
                 shape = Type, 
                 size = 10), 
             stroke = 0.5,
             fill = NA) + 
  scale_color_manual(values = c("AMF"="#000000",
                                "Bacteria"="#C52A20",
                                "Fungi"="#448DCD",
                                "Mix"="#E9DDD3",
                                "Multiple"="#F7AF34"),
                     name = "EWEs") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25),
                     name = "EWEs") + 
  theme(legend.position = c(0.1, 0.3))
# output  8*4


######################## Linear mixed effects model and plots
library(metafor)
library(boot)
library(parallel)

MWD<- read.csv('MWD.csv')
# Check data
head(MWD)

# 1. Check the number of observations
total_number <- nrow(MWD)
cat("Total number of observations in the dataset:", total_number, "\n")
# Total number of observations in the dataset: 464

# 2. Check the number of studies
unique_studyid_number <- length(unique(MWD$StudyID))
cat("Number of unique StudyID:", unique_studyid_number, "\n")
# Number of unique StudyID: 67



library(lme4) #package for mixed effect model
library(MuMIn)
library(lmerTest)
###Background Character
##########
##ExperimentalDuration
mdExperimentalDuration<-lmer(lnRR~ scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdExperimentalDuration)
# Number of obs: 458, groups:  StudyID, 65
anova(mdExperimentalDuration) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
# scale(ExperimentalDuration)  63.143  63.143     1 112.997  2.4256 0.12217  
# Inoculanttype               214.494  71.498     3  73.232  2.7465 0.04898 *


##Sand
mdtime<-lmer(lnRR~ scale(Sand)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 222, groups:  StudyID, 18
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Sand)                   8.212   8.212     1  15.599  0.2844 0.6013
# scale(ExperimentalDuration)  10.923  10.923     1 144.441  0.3783 0.5395
# Inoculanttype               178.822  59.607     3   9.617  2.0645 0.1711


##Silt
mdtime<-lmer(lnRR~ scale(Silt)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 222, groups:  StudyID, 18
anova(mdtime) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(Silt)                   0.007   0.007     1  31.155  0.0002 0.9880
# scale(ExperimentalDuration)   6.547   6.547     1 139.370  0.2266 0.6348
# Inoculanttype               167.812  55.937     3   9.934  1.9365 0.1881  


##Clay
mdtime<-lmer(lnRR~ scale(Clay)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 222, groups:  StudyID, 18
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(Clay)                  51.099  51.099     1  9.867  1.7568 0.2149
# scale(ExperimentalDuration)  15.468  15.468     1 72.311  0.5318 0.4682
# Inoculanttype               262.554  87.518     3  6.524  3.0089 0.1095


##ColonizationRate
mdtime<-lmer(lnRR~ scale(ColonizationRate)+scale(ExperimentalDuration)+(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 163, groups:  StudyID, 44
anova(mdtime) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)   
# scale(ColonizationRate)     225.924 225.924     1 152.852  7.8394 0.005772 **
# scale(ExperimentalDuration)  10.217  10.217     1  45.472  0.3545 0.554514


##TGRSPCK
mdtime<-lmer(lnRR~ scale(TGRSPCK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 144, groups:  StudyID, 30
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(TGRSPCK)               7.9644  7.9644     1  32.209  0.5119 0.4795
# scale(ExperimentalDuration)  1.3766  1.3766     1  47.006  0.0885 0.7674
# Inoculanttype               30.3002 30.3002     1 135.574  1.9473 0.1652


##EEGRSPCK
mdtime<-lmer(lnRR~ scale(EEGRSPCK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 132, groups:  StudyID, 29
anova(mdtime) 
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
# scale(EEGRSPCK)             41.360  41.360     1  19.196  2.4887 0.1310
# scale(ExperimentalDuration)  1.363   1.363     1  20.671  0.0820 0.7774
# Inoculanttype               14.247  14.247     1 120.590  0.8572 0.3564


##pHCK
mdtime<-lmer(lnRR~ scale(pHCK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 32, groups:  StudyID, 5
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# scale(pHCK)                 6.7505  6.7505     1    27  0.8277 0.3710
# scale(ExperimentalDuration) 4.9661  4.9661     1    27  0.6089 0.4420
# Inoculanttype               0.4103  0.2052     2    27  0.0252 0.9752


##SOCCK
mdtime<-lmer(lnRR~ scale(SOCCK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 137, groups:  StudyID, 25
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# scale(SOCCK)                83.893  83.893     1 33.357  2.5155 0.1222
# scale(ExperimentalDuration)  9.820   9.820     1 13.381  0.2944 0.5963
# Inoculanttype                5.812   2.906     2 86.503  0.0871 0.9166


##TNCK
mdtime<-lmer(lnRR~ scale(TNCK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 32, groups:  StudyID, 4
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
# scale(TNCK)                 707.58  707.58     1 26.8098 20.7626 0.0001017 ***
# scale(ExperimentalDuration) 184.65  184.65     1  2.8257  5.4182 0.1077532    
# Inoculanttype               136.79   68.39     2 25.2959  2.0069 0.1552433 


##NO3CK
mdtime<-lmer(lnRR~ scale(NO3CK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 12, groups:  StudyID, 3
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
# scale(NO3CK)                21.438  21.438     1 0.1418  2.7805 0.7381  
# scale(ExperimentalDuration)  1.741   1.741     1 0.1115  0.2258 0.8832  
# Inoculanttype               41.561  41.561     1 7.0909  5.3903 0.0528 .


##NH4CK
mdtime<-lmer(lnRR~ scale(NH4CK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 12, groups:  StudyID, 3
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
# scale(NH4CK)                35.052  35.052     1     8  4.7449 0.06103 .
# scale(ExperimentalDuration) 16.065  16.065     1     8  2.1747 0.17853  
# Inoculanttype               42.377  42.377     1     8  5.7366 0.04351 *


##APCK
mdtime<-lmer(lnRR~ scale(APCK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 50, groups:  StudyID, 10
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# scale(APCK)                 579.30  579.30     1    45 15.0863 0.0003339 ***
# scale(ExperimentalDuration)   0.04    0.04     1    45  0.0011 0.9733092    
# Inoculanttype               283.37  141.69     2    45  3.6898 0.0328184 *  


##AKCK
mdtime<-lmer(lnRR~ scale(AKCK)+scale(ExperimentalDuration)+Inoculanttype +(1|StudyID), weights=Wr, data=MWD)
summary(mdtime)
# Number of obs: 48, groups:  StudyID, 9
anova(mdtime)
# Type III Analysis of Variance Table with Satterthwaite's method
#                              Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
# scale(AKCK)                 1286.00 1286.00     1 38.057 123.6206 1.622e-13 ***
# scale(ExperimentalDuration)   13.95   13.95     1  6.481   1.3406    0.2878    
# Inoculanttype                 23.49   11.75     2 36.920   1.1292    0.3342 



### Regression
library(ggplot2)    # For creating plots
library(ggpmisc)    # For stat_poly_eq and displaying equations
library(ggpubr)     # For stat_cor to calculate Pearson correlation
library(dplyr)      # Optional, useful for data manipulation

##Sand
sum(!is.na(MWD$Sand)) ## n = 222
p1 <- ggplot(MWD, aes(y=lnRR, x=Sand)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="Sand")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,  
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="Sand n = 222" , y="lnMWD")
p1
pdf("Sand.pdf",width=8,height=8)
p1
dev.off()

##Silt
sum(!is.na(MWD$Silt)) ## n = 222
p2 <- ggplot(MWD, aes(y=lnRR, x=Silt)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="Silt")+
  theme(panel.grid=element_blank())+ 
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="Silt n = 222" , y="lnMWD")
p2
pdf("Silt.pdf",width=8,height=8)
p2
dev.off()


##Clay
sum(!is.na(MWD$Clay)) ## n = 222
p3 <- ggplot(MWD, aes(y=lnRR, x=Clay)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="Clay")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="Clay n = 222" , y="lnMWD")
p3
pdf("Clay.pdf",width=8,height=8)
p3
dev.off()  

##Experimental Duration
sum(!is.na(MWD$ExperimentalDuration)) ## n = 458
p4 <- ggplot(MWD, aes(y=lnRR, x=ExperimentalDuration)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="ExperimentalDuration")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="ExperimentalDuration n = 458" , y="lnMWD")
p4
pdf("ExperimentalDuration.pdf",width=8,height=8)
p4
dev.off()    

##ColonizationRate
sum(!is.na(MWD$ColonizationRate)) ## n = 167
p5 <- ggplot(MWD, aes(y=lnRR, x=ColonizationRate)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="ColonizationRate")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="ColonizationRate n = 167" , y="lnMWD")
p5
pdf("ColonizationRate.pdf",width=8,height=8)
p5
dev.off()    


##TGRSPCK
sum(!is.na(MWD$TGRSPCK)) ## n = 144
p6 <- ggplot(MWD, aes(y=lnRR, x=TGRSPCK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="TGRSPCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="TGRSPCK n = 144" , y="lnMWD")
p6
pdf("TGRSPCK.pdf",width=8,height=8)
p6
dev.off()    


##EEGRSPCK
sum(!is.na(MWD$EEGRSPCK)) ## n = 132
p7 <- ggplot(MWD, aes(y=lnRR, x=EEGRSPCK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="EEGRSPCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="EEGRSPCK n = 132" , y="lnMWD")
p7
pdf("EEGRSPCK.pdf",width=8,height=8)
p7
dev.off()    

##pHCK
sum(!is.na(MWD$pHCK)) ## n = 36
p8 <- ggplot(MWD, aes(y=lnRR, x=pHCK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="pHCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="pHCK n = 36" , y="lnMWD")
p8
pdf("pHCK.pdf",width=8,height=8)
p8
dev.off()    



##SOCCK
sum(!is.na(MWD$SOCCK)) ## n = 139
p9 <- ggplot(MWD, aes(y=lnRR, x=SOCCK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="SOCCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="SOCCK n = 139" , y="lnMWD")
p9
pdf("SOCCK.pdf",width=8,height=8)
p9
dev.off()    


##TNCK
sum(!is.na(MWD$TNCK)) ## n = 38
p10 <- ggplot(MWD, aes(y=lnRR, x=TNCK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="TNCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="TNCK n = 38" , y="lnMWD")
p10
pdf("TNCK.pdf",width=8,height=8)
p10
dev.off()    


##NO3CK
sum(!is.na(MWD$NO3CK)) ## n = 14
p11 <- ggplot(MWD, aes(y=lnRR, x=NO3CK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="NO3CK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="NO3CK n = 14" , y="lnMWD")
p11
pdf("NO3CK.pdf",width=8,height=8)
p11
dev.off()    


##NH4CK
sum(!is.na(MWD$NH4CK)) ## n = 14
p12 <- ggplot(MWD, aes(y=lnRR, x=NH4CK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="NH4CK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="NH4CK n = 14" , y="lnMWD")
p12
pdf("NH4CK.pdf",width=8,height=8)
p12
dev.off()    


##APCK
sum(!is.na(MWD$APCK)) ## n = 56
p13 <- ggplot(MWD, aes(y=lnRR, x=APCK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="APCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="APCK n = 56" , y="lnMWD")
p13
pdf("APCK.pdf",width=8,height=8)
p13
dev.off()    


##AKCK
sum(!is.na(MWD$AKCK)) ## n = 54
p14 <- ggplot(MWD, aes(y=lnRR, x=AKCK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="AKCK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="AKCK n = 54" , y="lnMWD")
p14
pdf("AKCK.pdf",width=8,height=8)
p14
dev.off()    


##lnRRTGRSP
sum(!is.na(MWD$lnRRTGRSP)) ## n = 144
p15 <- ggplot(MWD, aes(y=lnRR, x=lnRRTGRSP)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRTGRSP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRTGRSP n = 144" , y="lnMWD")
p15
pdf("lnRRTGRSP.pdf",width=8,height=8)
p15
dev.off()    


##lnRREEGRSP
sum(!is.na(MWD$lnRREEGRSP)) ## n = 132
p16 <- ggplot(MWD, aes(y=lnRR, x=lnRREEGRSP)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRREEGRSP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRREEGRSP n = 132" , y="lnMWD")
p16
pdf("lnRREEGRSP.pdf",width=8,height=8)
p16
dev.off()    


##lnRRpH
sum(!is.na(MWD$lnRRpH)) ## n = 36
p17 <- ggplot(MWD, aes(y=lnRR, x=lnRRpH)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRpH")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRpH n = 36" , y="lnMWD")
p17
pdf("lnRRpH.pdf",width=8,height=8)
p17
dev.off()    


##lnRRSOC
sum(!is.na(MWD$lnRRSOC)) ## n = 139
p18 <- ggplot(MWD, aes(y=lnRR, x=lnRRSOC)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRSOC")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRSOC n = 139" , y="lnMWD")
p18
pdf("lnRRSOC.pdf",width=8,height=8)
p18
dev.off()    


##lnRRTN
sum(!is.na(MWD$lnRRTN)) ## n = 38
p19 <- ggplot(MWD, aes(y=lnRR, x=lnRRTN)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRTN")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRTN n = 38" , y="lnMWD")
p19
pdf("lnRRTN.pdf",width=8,height=8)
p19
dev.off()    


##lnRRNO3
sum(!is.na(MWD$lnRRNO3)) ## n = 14
p20 <- ggplot(MWD, aes(y=lnRR, x=lnRRNO3)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRNO3")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRNO3 n = 14" , y="lnMWD")
p20
pdf("lnRRNO3.pdf",width=8,height=8)
p20
dev.off()    


##lnRRNH4
sum(!is.na(MWD$lnRRNH4)) ## n = 14
p21 <- ggplot(MWD, aes(y=lnRR, x=lnRRNH4)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRNH4")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRNH4 n = 14" , y="lnMWD")
p21
pdf("lnRRNH4.pdf",width=8,height=8)
p21
dev.off()    


##lnRRAP
sum(!is.na(MWD$lnRRAP)) ## n = 56
p22 <- ggplot(MWD, aes(y=lnRR, x=lnRRAP)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRAP")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRAP n = 56" , y="lnMWD")
p22
pdf("lnRRAP.pdf",width=8,height=8)
p22
dev.off()    


##lnRRAK
sum(!is.na(MWD$lnRRAK)) ## n = 54
p23 <- ggplot(MWD, aes(y=lnRR, x=lnRRAK)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRAK")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRAK n = 54" , y="lnMWD")
p23
pdf("lnRRAK.pdf",width=8,height=8)
p23
dev.off()    


##lnRRUrease
sum(!is.na(MWD$lnRRUrease)) ## n = 11
p24 <- ggplot(MWD, aes(y=lnRR, x=lnRRUrease)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRUrease")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRUrease n = 11" , y="lnMWD")
p24
pdf("lnRRUrease.pdf",width=8,height=8)
p24
dev.off()    


##lnRRInvertase
sum(!is.na(MWD$lnRRInvertase)) ## n = 9
p25 <- ggplot(MWD, aes(y=lnRR, x=lnRRInvertase)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRInvertase")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRInvertase n = 9" , y="lnMWD")
p25
pdf("lnRRInvertase.pdf",width=8,height=8)
p25
dev.off()    


##lnRRCatalase
sum(!is.na(MWD$lnRRCatalase)) ## n = 10
p26 <- ggplot(MWD, aes(y=lnRR, x=lnRRCatalase)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRCatalase")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRCatalase n = 10" , y="lnMWD")
p26
pdf("lnRRCatalase.pdf",width=8,height=8)
p26
dev.off()    


##lnRRAlkp
sum(!is.na(MWD$lnRRAlkp)) ## n = 12
p27 <- ggplot(MWD, aes(y=lnRR, x=lnRRAlkp)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRAlkp")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRAlkp n = 12" , y="lnMWD")
p27
pdf("lnRRAlkp.pdf",width=8,height=8)
p27
dev.off()    


##lnRRAcip
sum(!is.na(MWD$lnRRAcip)) ## n = 12
p28 <- ggplot(MWD, aes(y=lnRR, x=lnRRAcip)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRAcip")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRAcip n = 12" , y="lnMWD")
p28
pdf("lnRRAcip.pdf",width=8,height=8)
p28
dev.off()    


##lnRRDH
sum(!is.na(MWD$lnRRDH)) ## n = 0
p29 <- ggplot(MWD, aes(y=lnRR, x=lnRRDH)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRDH")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRDH n = 220" , y="lnMWD")
p29
pdf("lnRRDH.pdf",width=8,height=8)
p29
dev.off()    


##lnRRBG
sum(!is.na(MWD$lnRRBG)) ## n = 6
p30 <- ggplot(MWD, aes(y=lnRR, x=lnRRBG)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRBG")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRBG n = 6" , y="lnMWD")
p30
pdf("lnRRBG.pdf",width=8,height=8)
p30
dev.off()    


##lnRRRootWeight
sum(!is.na(MWD$lnRRRootWeight)) ## n = 170
p31 <- ggplot(MWD, aes(y=lnRR, x=lnRRRootWeight)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRRootWeight")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRRootWeight n = 170" , y="lnMWD")
p31
pdf("lnRRRootWeight.pdf",width=8,height=8)
p31
dev.off()    


##lnRRShootWeight
sum(!is.na(MWD$lnRRShootWeight)) ## n = 154
p32 <- ggplot(MWD, aes(y=lnRR, x=lnRRShootWeight)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRShootWeight")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRShootWeight n = 154" , y="lnMWD")
p32
pdf("lnRRShootWeight.pdf",width=8,height=8)
p32
dev.off()    


##lnRRPlantWight
sum(!is.na(MWD$lnRRPlantWight)) ## n = 165
p33 <- ggplot(MWD, aes(y=lnRR, x=lnRRPlantWight)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRPlantWight")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5,
    label.x = 0.05,
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRPlantWight n = 165" , y="lnMWD")
p33
pdf("lnRRPlantWight.pdf",width=8,height=8)
p33
dev.off()    


##lnRRPlantHeight
sum(!is.na(MWD$lnRRPlantHeight)) ## n = 88
p34 <- ggplot(MWD, aes(y=lnRR, x=lnRRPlantHeight)) +
  geom_point(color="gray", size=15, shape=21) +geom_smooth(method=lm , color="black", linewidth=2.0, se=TRUE, level=0.95, linetype="dashed") + theme_bw()+theme(text = element_text(family = "serif",
                                                                                                                                                                                    size=20))+
  geom_hline(aes(yintercept=0), colour="black", linewidth=0.5, linetype="dashed")+
  labs(y="lnMWD", x="lnRRPlantHeight")+
  theme(panel.grid=element_blank())+
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~')),
    formula = y ~ x,  parse = TRUE,color="black",
    size = 5, 
    label.x = 0.05, 
    label.y = 0.85) + stat_cor(method = "pearson", size = 5) +
  labs(x="lnRRPlantHeight n = 88" , y="lnMWD")
p34
pdf("lnRRPlantHeight.pdf",width=8,height=8)
p34
dev.off()    



########################################################### Greenhouse experiment
####################### CPCoA
cpcoa <- function (otutab, metadata, dis = "bray", groupID = "Group", 
                   ellipse = T, label = F) 
{
  p_list = c("ggplot2", "vegan", "ggrepel")
  for (p in p_list) {
    if (!requireNamespace(p)) {
      install.packages(p)
    }
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
  idx = rownames(metadata) %in% colnames(otutab)
  metadata = metadata[idx, , drop = F]
  otutab = otutab[, rownames(metadata)]
  sampFile = as.data.frame(metadata[, groupID], row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"
  if (length(unique(sampFile$group)) > 2) {
    variability_table = function(cca) {
      chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
      variability_table = cbind(chi, chi/chi[1])
      colnames(variability_table) = c("inertia", "proportion")
      rownames(variability_table) = c("total", "constrained", 
                                      "unconstrained")
      return(variability_table)
    }
    capscale.gen = capscale(t(otutab) ~ group, data = sampFile, 
                            add = F, sqrt.dist = T, distance = dis)
    set.seed(123)
    perm_anova.gen = anova.cca(capscale.gen, permutations = 1000, 
                               parallel = 4)
    var_tbl.gen = variability_table(capscale.gen)
    eig = capscale.gen$CCA$eig
    variance = var_tbl.gen["constrained", "proportion"]
    p.val = perm_anova.gen[1, 4]
    points = as.data.frame(capscale.gen$CCA$wa)
    points = cbind(sampFile, points[rownames(points), ])
    
  } 
  return(list(eig = eig, points = points, var = variance, p = p.val))
}
library(ggplot2)
library(vegan)
library(ape)
library(plyr)
library(openxlsx)
otu <- read.csv('genes.csv', row.names = 1)
metadata <- read.csv('design.csv')
colnames(metadata) <- c("samples", "Treatment", "Group","Inter", "Duration" )
otu <- otu[, metadata$samples]
metadata <- metadata[metadata$samples %in% colnames(otu), ]
sampFile <- as.data.frame(metadata$Treatment, row.names = metadata$samples)
colnames(sampFile)[1] <- "treatment"
cpcoa_result <- capscale(t(otu) ~ treatment, data = sampFile, add = FALSE, sqrt.dist = TRUE, distance = "bray")
cppoints <- as.data.frame(cpcoa_result$CCA$wa)
cppoints$samples <- row.names(cppoints)
colnames(cppoints) <- c("CPCoA1", "CPCoA2", "CPCoA3", "CPCoA4", "CPCoA5", "samples")
cppoints <- cppoints[, c("CPCoA1", "CPCoA2", "samples")]
eigenvalues <- cpcoa_result$CCA$eig
pca1 <- round(100 * eigenvalues[1] / sum(eigenvalues), 2)
pca2 <- round(100 * eigenvalues[2] / sum(eigenvalues), 2)
variability_table <- function(cca) {
  chi <- c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi/chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained")
  return(variability_table)
}

var_tbl <- variability_table(cpcoa_result)
variance <- var_tbl["constrained", "proportion"]
p_val <- anova.cca(cpcoa_result, permutations = 1000)[1, 4]
df <- merge(cppoints, metadata, by = "samples")
color <- c('black','black','black',"#F7AF34", "#448DCD", "#ffcd85", "#a1b8d5", "#FFE5BE",'#CDD9E8')
shapes <- c(17, 15, 16, 18)
p1 <- ggplot(df, aes(CPCoA1, CPCoA2, shape = Group)) +
  geom_point(aes(color = Treatment), size = 8) +  #
  xlab(paste("CPCoA1 (", pca1, "%)", sep = "")) + 
  ylab(paste("CPCoA2 (", pca2, "%)", sep = "")) +
  stat_ellipse(aes(group = Inter, color = Inter), level = 0.95, linetype = "dashed", linewidth = 1) +  # 
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = color) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black', size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Variance explained: ", round(variance * 100, 2), "%; P = ", format(p_val, digits = 2), sep = ""))

p1
################# PERMANOVA
library(vegan)
otutab <- read.csv("otutab.csv", row.names = 1)
metadata <- read.csv("treatment.csv")
rownames(metadata) <- metadata$Treatment 
otutab <- otutab[, rownames(metadata)]    
dist_matrix <- vegdist(t(otutab), method = "bray")
perm_model_interaction <- adonis2(dist_matrix ~ (Year + Plant + Inoculant + Aggregate)^2, 
                                  data = metadata, 
                                  strata = metadata$Random,
                                  permutations = 999)
print(perm_model_interaction)
perm_results <- as.data.frame(perm_model_interaction)
write.csv(perm_results, "permanova_results.csv", row.names = TRUE)

################ Control and microbial inoculant
otu <- read.csv('otutab.csv', row.names = 1)
otu <- data.frame(t(otu))
otu.distance <- vegdist(otu)
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p

group <- read.csv('design.csv')
colnames(group) <- c("samples","group")
df <- merge(pc12,group,by="samples")
color=c("#F7AF34","#448DCD")
p1<-ggplot(data=df,aes(x=V1,y=V2,
                       color=group,shape=group))+
  theme_bw()+
  geom_point(size=4, shape=19)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+
  #guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"))+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#F7AF34","#448DCD"))+
  theme(axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank())
p1
p1 + stat_ellipse(data=df,geom = "polygon",level=0.9,linetype = 2,size=0.5,aes(fill=group),alpha=0.2,show.legend = T)
write.table(pc12, 'pc12.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

##################################### GLM
library(car)
library(fitdistrplus)
library(e1071)
library(lmtest)
tbl.GLM <- read.csv("GLM.csv")

#################################  BacterialRichness_model  
table(tbl.GLM$BacterialRichness)
### test mean 5485.078
mean(tbl.GLM$BacterialRichness)
### test variance 125413.2
var(tbl.GLM$BacterialRichness)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonBacterialRichness <- glm(BacterialRichness ~ Inoculants * Aggregate * Plant * Time, 
                                    data = tbl.GLM, family = poisson)
### The discrete factor was calculated  7.49899
deviance(glm_poissonBacterialRichness) / df.residual(glm_poissonBacterialRichness)
### Distribution of response variables
hist(tbl.GLM$BacterialRichness, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$BacterialRichness, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiBacterialRichness <- glm(BacterialRichness ~ Inoculants * Aggregate * Plant * Time, 
                                  data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiBacterialRichness), residuals(glm_quasiBacterialRichness), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiBacterialRichness, type = "II")

glm_quasiBacterialRichness_simplified <- glm(BacterialRichness ~ Inoculants + Aggregate + Plant + Time +
                                               Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                               Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiBacterialRichness_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacterialRichness
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              0.272  1   0.602149    
# Aggregate             123.126  2  < 2.2e-16 ***
# Plant                  32.007  1  1.537e-08 ***
# Time                   26.121  1  3.207e-07 ***
# Inoculants:Aggregate   10.407  2   0.005496 ** 
# Inoculants:Plant        0.003  1   0.955039    
# Aggregate:Plant         3.146  2   0.207375    
# Inoculants:Time         0.000  1   0.989132    
# Aggregate:Time          7.436  2   0.024278 *  
# Plant:Time                     0             


################################   BacterialShannon_model  
table(tbl.GLM$BacterialShannon)
### test mean 7.173791
mean(tbl.GLM$BacterialShannon)
### test variance 0.02037439
var(tbl.GLM$BacterialShannon)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$BacterialShannon, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$BacterialShannon) ## -0.3825057
kurtosis(tbl.GLM$BacterialShannon) ## -0.4248509
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianBacterialShannon <- glm(BacterialShannon ~ Inoculants * Aggregate * Plant * Time, 
                                    data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianBacterialShannon))
qqline(residuals(glm_gaussianBacterialShannon), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianBacterialShannon))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianBacterialShannon)
# W = 0.99269, p-value = 0.8612

## Test Homoscedasticity
plot(fitted(glm_gaussianBacterialShannon), residuals(glm_gaussianBacterialShannon), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianBacterialShannon)
# studentized Breusch-Pagan test
# data:  glm_gaussianBacterialShannon
# BP = 30.556, df = 17, p-value = 0.0226

#### log transformation
tbl.GLM$BacterialShannon_log <- log(tbl.GLM$BacterialShannon)
hist(tbl.GLM$BacterialShannon_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedBacterialShannon_log <- glm(BacterialShannon_log ~ Inoculants * Aggregate * Plant * Time, 
                                           data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedBacterialShannon_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedBacterialShannon_log)
# W = 0.99228, p-value = 0.8323
## Test Homoscedasticity
plot(fitted(glm_transformedBacterialShannon_log), residuals(glm_transformedBacterialShannon_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedBacterialShannon_log)
# studentized Breusch-Pagan test
# data:  glm_transformedBacterialShannon_log
# BP = 31.431, df = 17, p-value = 0.01769

#### sqrt transformation
tbl.GLM$BacterialShannon_sqrt <- sqrt(tbl.GLM$BacterialShannon)
hist(tbl.GLM$BacterialShannon_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedBacterialShannon_sqrt <- glm(BacterialShannon_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                            data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedBacterialShannon_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedBacterialShannon_sqrt)
# W = 0.99251, p-value = 0.8488
## Test Homoscedasticity
plot(fitted(glm_transformedBacterialShannon_sqrt), residuals(glm_transformedBacterialShannon_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedBacterialShannon_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedBacterialShannon_sqrt
# BP = 31.002, df = 17, p-value = 0.01996

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianBacterialShannon)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -182.9488  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedBacterialShannon_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -582.8866  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedBacterialShannon_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -524.3358  
### Choose log

###### Data transformation (log transformation)
tbl.GLM$BacterialShannon_log <- log(tbl.GLM$BacterialShannon)
glm_transformedBacterialShannon <- glm(BacterialShannon_log ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = gaussian)
Anova(glm_transformedBacterialShannon, type = "II") 

glm_transformedBacterialShannon_simplified <- glm(BacterialShannon_log ~ Inoculants + Aggregate + Plant + Time +
                                                    Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                    Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_transformedBacterialShannon_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacterialShannon_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              1.311  1    0.25214    
# Aggregate             119.247  2  < 2.2e-16 ***
# Plant                   7.931  1    0.00486 ** 
# Time                   28.285  1  1.047e-07 ***
# Inoculants:Aggregate    3.320  2    0.19014    
# Inoculants:Plant        0.022  1    0.88314    
# Aggregate:Plant         1.496  2    0.47339    
# Inoculants:Time         1.078  1    0.29904    
# Aggregate:Time          4.266  2    0.11850    
# Plant:Time                     0               



#################################  BacterialChao_model  
table(tbl.GLM$BacterialChao)
### test mean 6866.007
mean(tbl.GLM$BacterialChao)
### test variance 215640.1
var(tbl.GLM$BacterialChao)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonBacterialChao <- glm(BacterialChao ~ Inoculants * Aggregate * Plant * Time, 
                                data = tbl.GLM, family = poisson)
### The discrete factor was calculated  15.38836
deviance(glm_poissonBacterialChao) / df.residual(glm_poissonBacterialChao)
### Distribution of response variables
hist(tbl.GLM$BacterialChao, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$BacterialChao, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiBacterialChao <- glm(BacterialChao ~ Inoculants * Aggregate * Plant * Time, 
                              data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiBacterialChao), residuals(glm_quasiBacterialChao), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiBacterialChao, type = "II")

glm_quasiBacterialChao_simplified <- glm(BacterialChao ~ Inoculants + Aggregate + Plant + Time +
                                           Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                           Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiBacterialChao_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacterialChao
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              0.215  1   0.642954    
# Aggregate              52.867  2  3.312e-12 ***
# Plant                  21.059  1  4.454e-06 ***
# Time                    8.830  1   0.002963 ** 
# Inoculants:Aggregate    9.199  2   0.010055 *  
# Inoculants:Plant        0.685  1   0.407758    
# Aggregate:Plant         2.991  2   0.224108    
# Inoculants:Time         0.724  1   0.394942    
# Aggregate:Time          6.115  2   0.047006 *  
# Plant:Time                     0   


#################################  BacterialACE_model  
table(tbl.GLM$BacterialACE)
### test mean 6937.395
mean(tbl.GLM$BacterialACE)
### test variance 224578.2
var(tbl.GLM$BacterialACE)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonBacterialACE <- glm(BacterialACE ~ Inoculants * Aggregate * Plant * Time, 
                               data = tbl.GLM, family = poisson)
### The discrete factor was calculated  16.5113
deviance(glm_poissonBacterialACE) / df.residual(glm_poissonBacterialACE)
### Distribution of response variables
hist(tbl.GLM$BacterialACE, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$BacterialACE, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiBacterialACE <- glm(BacterialACE ~ Inoculants * Aggregate * Plant * Time, 
                             data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiBacterialACE), residuals(glm_quasiBacterialACE), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiBacterialACE, type = "II")

glm_quasiBacterialACE_simplified <- glm(BacterialACE ~ Inoculants + Aggregate + Plant + Time +
                                          Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                          Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiBacterialACE_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacterialACE
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              0.403  1   0.525791    
# Aggregate              47.394  2  5.112e-11 ***
# Plant                  18.113  1  2.082e-05 ***
# Time                    6.125  1   0.013329 *  
# Inoculants:Aggregate    9.670  2   0.007948 ** 
# Inoculants:Plant        0.686  1   0.407401    
# Aggregate:Plant         2.836  2   0.242230    
# Inoculants:Time         0.630  1   0.427505    
# Aggregate:Time          5.368  2   0.068297 .  
# Plant:Time                     0               


#################################  FungalRichness_model  
table(tbl.GLM$FungalRichness)
### test mean 727.5392
mean(tbl.GLM$FungalRichness)
### test variance 5684.251
var(tbl.GLM$FungalRichness)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonFungalRichness <- glm(FungalRichness ~ Inoculants * Aggregate * Plant * Time, 
                                 data = tbl.GLM, family = poisson)
### The discrete factor was calculated  5.538844
deviance(glm_poissonFungalRichness) / df.residual(glm_poissonFungalRichness)
### Distribution of response variables
hist(tbl.GLM$FungalRichness, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$FungalRichness, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiFungalRichness <- glm(FungalRichness ~ Inoculants * Aggregate * Plant * Time, 
                               data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiFungalRichness), residuals(glm_quasiFungalRichness), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiFungalRichness, type = "II")

glm_quasiFungalRichness_simplified <- glm(FungalRichness ~ Inoculants + Aggregate + Plant + Time +
                                            Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                            Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiFungalRichness_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FungalRichness
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             1.6066  1   0.204973    
# Aggregate             31.3424  2  1.563e-07 ***
# Plant                  4.5024  1   0.033847 *  
# Time                   6.9357  1   0.008449 ** 
# Inoculants:Aggregate   0.7622  2   0.683117    
# Inoculants:Plant       0.0104  1   0.918681    
# Aggregate:Plant        3.1788  2   0.204044    
# Inoculants:Time        0.0567  1   0.811772    
# Aggregate:Time         5.9246  2   0.051699 .  
# Plant:Time                     0            



#################################  FungalShannon_model  
table(tbl.GLM$FungalShannon)
### test mean 4.085436
mean(tbl.GLM$FungalShannon)
### test variance 0.173
var(tbl.GLM$FungalShannon)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$FungalShannon, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$FungalShannon) ## -0.7305734
kurtosis(tbl.GLM$FungalShannon) ## 0.8088136
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianFungalShannon <- glm(FungalShannon ~ Inoculants * Aggregate * Plant * Time, 
                                 data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianFungalShannon))
qqline(residuals(glm_gaussianFungalShannon), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianFungalShannon))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianFungalShannon)
# W = 0.96078, p-value = 0.004072, do not meet standard

## Test Homoscedasticity
plot(fitted(glm_gaussianFungalShannon), residuals(glm_gaussianFungalShannon), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianFungalShannon)
# studentized Breusch-Pagan test
# data:  glm_gaussianFungalShannon
# BP = 21.867, df = 17, p-value = 0.1899

#### log transformation
tbl.GLM$FungalShannon_log <- log(tbl.GLM$FungalShannon)
hist(tbl.GLM$FungalShannon_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedFungalShannon_log <- glm(FungalShannon_log ~ Inoculants * Aggregate * Plant * Time, 
                                        data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedFungalShannon_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedFungalShannon_log)
# 0.93235, p-value = 5.746e-05
## Test Homoscedasticity
plot(fitted(glm_transformedFungalShannon_log), residuals(glm_transformedFungalShannon_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedFungalShannon_log)
# studentized Breusch-Pagan test
# data:  glm_transformedFungalShannon_log
# 23.041, df = 17, p-value = 0.1479

#### sqrt transformation
tbl.GLM$FungalShannon_sqrt <- sqrt(tbl.GLM$FungalShannon)
hist(tbl.GLM$FungalShannon_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedFungalShannon_sqrt <- glm(FungalShannon_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                         data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedFungalShannon_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedFungalShannon_sqrt)
# W = 0.94802, p-value = 0.0005342
## Test Homoscedasticity
plot(fitted(glm_transformedFungalShannon_sqrt), residuals(glm_transformedFungalShannon_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedFungalShannon_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedFungalShannon_sqrt
# BP = 22.54, df = 17, p-value = 0.1648

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianFungalShannon)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 83.51507  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedFungalShannon_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -187.461  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedFungalShannon_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -194.043  
### Choose log

###### Data transformation (sqrt transformation)
tbl.GLM$FungalShannon_sqrt <- sqrt(tbl.GLM$FungalShannon)
glm_transformedFungalShannon <- glm(FungalShannon_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                    data = tbl.GLM, family = gaussian)
Anova(glm_transformedFungalShannon, type = "II") 

glm_transformedFungalShannon_simplified <- glm(FungalShannon_sqrt ~ Inoculants + Aggregate + Plant + Time +
                                                 Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                 Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_transformedFungalShannon_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FungalShannon_sqrt
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             5.5316  1  0.0186765 *  
# Aggregate             12.0901  2  0.0023696 ** 
# Plant                  4.1958  1  0.0405254 *  
# Time                  10.9770  1  0.0009225 ***
# Inoculants:Aggregate   6.3203  2  0.0424188 *  
# Inoculants:Plant       2.4603  1  0.1167596    
# Aggregate:Plant        1.7814  2  0.4103668    
# Inoculants:Time        0.0434  1  0.8350072    
# Aggregate:Time         1.4839  2  0.4761891    
# Plant:Time                     0            



#################################  FungalChao_model  
table(tbl.GLM$FungalChao)
### test mean 910.7075
mean(tbl.GLM$FungalChao)
### test variance 9458.757
var(tbl.GLM$FungalChao)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonFungalChao <- glm(FungalChao ~ Inoculants * Aggregate * Plant * Time, 
                             data = tbl.GLM, family = poisson)
### The discrete factor was calculated  8.105534
deviance(glm_poissonFungalChao) / df.residual(glm_poissonFungalChao)
### Distribution of response variables
hist(tbl.GLM$FungalChao, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$FungalChao, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiFungalChao <- glm(FungalChao ~ Inoculants * Aggregate * Plant * Time, 
                           data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiFungalChao), residuals(glm_quasiFungalChao), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiFungalChao, type = "II")

glm_quasiFungalChao_simplified <- glm(FungalChao ~ Inoculants + Aggregate + Plant + Time +
                                        Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                        Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiFungalChao_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FungalChao
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             0.6622  1  0.4157980    
# Aggregate             14.2498  2  0.0008048 ***
# Plant                  4.0581  1  0.0439596 *  
# Time                   5.0883  1  0.0240874 *  
# Inoculants:Aggregate   0.9791  2  0.6129135    
# Inoculants:Plant       0.0335  1  0.8548173    
# Aggregate:Plant        3.7847  2  0.1507191    
# Inoculants:Time        0.0022  1  0.9628947    
# Aggregate:Time         4.5693  2  0.1018103    
# Plant:Time                     0               


#################################  FungalACE_model  
table(tbl.GLM$FungalACE)
### test mean 911.331
mean(tbl.GLM$FungalACE)
### test variance 8838.838
var(tbl.GLM$FungalACE)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonFungalACE <- glm(FungalACE ~ Inoculants * Aggregate * Plant * Time, 
                            data = tbl.GLM, family = poisson)
### The discrete factor was calculated  7.370136
deviance(glm_poissonFungalACE) / df.residual(glm_poissonFungalACE)
### Distribution of response variables
hist(tbl.GLM$FungalACE, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$FungalACE, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiFungalACE <- glm(FungalACE ~ Inoculants * Aggregate * Plant * Time, 
                          data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiFungalACE), residuals(glm_quasiFungalACE), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiFungalACE, type = "II")

glm_quasiFungalACE_simplified <- glm(FungalACE ~ Inoculants + Aggregate + Plant + Time +
                                       Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                       Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiFungalACE_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FungalACE
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             0.8556  1  0.3549634    
# Aggregate             14.3474  2  0.0007665 ***
# Plant                  4.3671  1  0.0366388 *  
# Time                   6.7767  1  0.0092356 ** 
# Inoculants:Aggregate   0.9748  2  0.6142326    
# Inoculants:Plant       0.0391  1  0.8433275    
# Aggregate:Plant        5.4442  2  0.0657369 .  
# Inoculants:Time        0.0589  1  0.8082947    
# Aggregate:Time         3.3092  2  0.1911703    
# Plant:Time                     0           

#################################  ProteobacteriaPhylum_model  
table(tbl.GLM$ProteobacteriaPhylum)
### test mean 25.40098
mean(tbl.GLM$ProteobacteriaPhylum)
### test variance 12.53277
var(tbl.GLM$ProteobacteriaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ProteobacteriaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ProteobacteriaPhylum) ## 0.4924322
kurtosis(tbl.GLM$ProteobacteriaPhylum) ## 1.828315
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianProteobacteriaPhylum <- glm(ProteobacteriaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                        data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianProteobacteriaPhylum))
qqline(residuals(glm_gaussianProteobacteriaPhylum), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianProteobacteriaPhylum))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianProteobacteriaPhylum)
# W = 0.97297, p-value = 0.03442

## Test Homoscedasticity
plot(fitted(glm_gaussianProteobacteriaPhylum), residuals(glm_gaussianProteobacteriaPhylum), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianProteobacteriaPhylum)
# studentized Breusch-Pagan test
# data:  glm_gaussianProteobacteriaPhylum
# BP = 29.668, df = 17, p-value = 0.02884

#### log transformation
tbl.GLM$ProteobacteriaPhylum_log <- log(tbl.GLM$ProteobacteriaPhylum)
hist(tbl.GLM$ProteobacteriaPhylum_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedProteobacteriaPhylum_log <- glm(ProteobacteriaPhylum_log ~ Inoculants * Aggregate * Plant * Time, 
                                               data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedProteobacteriaPhylum_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedProteobacteriaPhylum_log)
# W = 0.97635, p-value = 0.0639
## Test Homoscedasticity
plot(fitted(glm_transformedProteobacteriaPhylum_log), residuals(glm_transformedProteobacteriaPhylum_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedProteobacteriaPhylum_log)
# studentized Breusch-Pagan test
# data:  glm_transformedProteobacteriaPhylum_log
# BP = 29.17, df = 17, p-value = 0.03299

#### sqrt transformation
tbl.GLM$ProteobacteriaPhylum_sqrt <- sqrt(tbl.GLM$ProteobacteriaPhylum)
hist(tbl.GLM$ProteobacteriaPhylum_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedProteobacteriaPhylum_sqrt <- glm(ProteobacteriaPhylum_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                                data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedProteobacteriaPhylum_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedProteobacteriaPhylum_sqrt)
# W = 0.97354, p-value = 0.03821
## Test Homoscedasticity
plot(fitted(glm_transformedProteobacteriaPhylum_sqrt), residuals(glm_transformedProteobacteriaPhylum_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedProteobacteriaPhylum_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedProteobacteriaPhylum_sqrt
# BP = 28.654, df = 17, p-value = 0.03786

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianProteobacteriaPhylum)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 528.4517  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedProteobacteriaPhylum_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -129.3607  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedProteobacteriaPhylum_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: 56.96237 
### Choose log

###### Data transformation (log transformation)
tbl.GLM$ProteobacteriaPhylum_log <- log(tbl.GLM$ProteobacteriaPhylum + 1)
glm_transformedProteobacteriaPhylum <- glm(ProteobacteriaPhylum_log ~ Inoculants * Aggregate * Plant * Time, 
                                           data = tbl.GLM, family = gaussian)
Anova(glm_transformedProteobacteriaPhylum, type = "II") 

glm_transformedProteobacteriaPhylum_simplified <- glm(ProteobacteriaPhylum_log ~ Inoculants + Aggregate + Plant + Time +
                                                        Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                        Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_transformedProteobacteriaPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ProteobacteriaPhylum_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants            13.2578  1  0.0002714 ***
# Aggregate              9.2931  2  0.0095948 ** 
# Plant                  1.7046  1  0.1916843    
# Time                   6.0292  1  0.0140714 *  
# Inoculants:Aggregate   0.3782  2  0.8276998    
# Inoculants:Plant       2.7400  1  0.0978636 .  
# Aggregate:Plant        0.3498  2  0.8395382    
# Inoculants:Time        1.6562  1  0.1981225    
# Aggregate:Time         0.5315  2  0.7666242    
# Plant:Time                     0               


#################################  AcidobacteriaPhylum_model  
table(tbl.GLM$AcidobacteriaPhylum)
### test mean 40.47255
mean(tbl.GLM$AcidobacteriaPhylum)
### test variance 19.03389
var(tbl.GLM$AcidobacteriaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$AcidobacteriaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$AcidobacteriaPhylum) ## -0.1731708
kurtosis(tbl.GLM$AcidobacteriaPhylum) ## 0.5260452
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianAcidobacteriaPhylum <- glm(AcidobacteriaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianAcidobacteriaPhylum))
qqline(residuals(glm_gaussianAcidobacteriaPhylum), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianAcidobacteriaPhylum))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianAcidobacteriaPhylum)
# W = 0.9914, p-value = 0.7651

## Test Homoscedasticity
plot(fitted(glm_gaussianAcidobacteriaPhylum), residuals(glm_gaussianAcidobacteriaPhylum), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianAcidobacteriaPhylum)
# studentized Breusch-Pagan test
# data:  glm_gaussianAcidobacteriaPhylum
# BP = 20.854, df = 17, p-value = 0.2329

###### 
glm_gaussianAcidobacteriaPhylum <- glm(AcidobacteriaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = gaussian)
Anova(glm_gaussianAcidobacteriaPhylum, type = "II") 

glm_gaussianAcidobacteriaPhylum_simplified <- glm(AcidobacteriaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                    Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                    Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianAcidobacteriaPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: AcidobacteriaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             34.993  1  3.309e-09 ***
# Aggregate               9.032  2    0.01093 *  
# Plant                   0.511  1    0.47484    
# Time                    1.852  1    0.17356    
# Inoculants:Aggregate    0.471  2    0.79005    
# Inoculants:Plant        5.029  1    0.02492 *  
# Aggregate:Plant         1.339  2    0.51189    
# Inoculants:Time         1.428  1    0.23207    
# Aggregate:Time          1.048  2    0.59204    
# Plant:Time                     0              


#################################  ActinobacteriaPhylum_model  
table(tbl.GLM$ActinobacteriaPhylum)
### test mean 6.463431
mean(tbl.GLM$ActinobacteriaPhylum)
### test variance 2.877017
var(tbl.GLM$ActinobacteriaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ActinobacteriaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ActinobacteriaPhylum) ## 0.02788649
kurtosis(tbl.GLM$ActinobacteriaPhylum) ## -0.6516005
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianActinobacteriaPhylum <- glm(ActinobacteriaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                        data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianActinobacteriaPhylum))
qqline(residuals(glm_gaussianActinobacteriaPhylum), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianActinobacteriaPhylum))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianActinobacteriaPhylum)
# W = 0.97252, p-value = 0.03171

## Test Homoscedasticity
plot(fitted(glm_gaussianActinobacteriaPhylum), residuals(glm_gaussianActinobacteriaPhylum), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianActinobacteriaPhylum)
# studentized Breusch-Pagan test
# data:  glm_gaussianActinobacteriaPhylum
# BP = 20.436, df = 17, p-value = 0.2526

#### log transformation
tbl.GLM$ActinobacteriaPhylum_log <- log(tbl.GLM$ActinobacteriaPhylum)
hist(tbl.GLM$ActinobacteriaPhylum_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedActinobacteriaPhylum_log <- glm(ActinobacteriaPhylum_log ~ Inoculants * Aggregate * Plant * Time, 
                                               data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedActinobacteriaPhylum_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedActinobacteriaPhylum_log)
# W = 0.94278, p-value = 0.0002462
## Test Homoscedasticity
plot(fitted(glm_transformedActinobacteriaPhylum_log), residuals(glm_transformedActinobacteriaPhylum_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedActinobacteriaPhylum_log)
# studentized Breusch-Pagan test
# data:  glm_transformedActinobacteriaPhylum_log
# BP = 16.522, df = 17, p-value = 0.4872

#### sqrt transformation
tbl.GLM$ActinobacteriaPhylum_sqrt <- sqrt(tbl.GLM$ActinobacteriaPhylum)
hist(tbl.GLM$ActinobacteriaPhylum_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedActinobacteriaPhylum_sqrt <- glm(ActinobacteriaPhylum_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                                data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedActinobacteriaPhylum_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedActinobacteriaPhylum_sqrt)
# W = 0.96155, p-value = 0.00463
## Test Homoscedasticity
plot(fitted(glm_transformedActinobacteriaPhylum_sqrt), residuals(glm_transformedActinobacteriaPhylum_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedActinobacteriaPhylum_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedActinobacteriaPhylum_sqrt
# BP = 17.956, df = 17, p-value = 0.3916

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianActinobacteriaPhylum)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 353.8586 
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedActinobacteriaPhylum_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -11.90964 
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedActinobacteriaPhylum_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: 26.78649 
### Choose log

######  log Transformation
tbl.GLM$ActinobacteriaPhylum_log <- log(tbl.GLM$ActinobacteriaPhylum)
glm_transformedActinobacteriaPhylum_log <- glm(ActinobacteriaPhylum_log ~ Inoculants * Aggregate * Plant * Time, 
                                               data = tbl.GLM, family = gaussian)
Anova(glm_transformedActinobacteriaPhylum_log, type = "II") 

glm_transformedActinobacteriaPhylum_log_simplified <- glm(ActinobacteriaPhylum_log ~ Inoculants + Aggregate + Plant + Time +
                                                            Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                            Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_transformedActinobacteriaPhylum_log_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ActinobacteriaPhylum_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              6.937  1   0.008443 ** 
# Aggregate              11.934  2   0.002562 ** 
# Plant                  37.334  1  9.954e-10 ***
# Time                   43.213  1  4.909e-11 ***
# Inoculants:Aggregate    1.476  2   0.478130    
# Inoculants:Plant        0.001  1   0.971539    
# Aggregate:Plant         0.627  2   0.730879    
# Inoculants:Time         0.036  1   0.849716    
# Aggregate:Time          0.053  2   0.973979    
# Plant:Time                     0               



#################################  PlanctomycetesPhylum_model  
table(tbl.GLM$PlanctomycetesPhylum)
### test mean 5.055196
mean(tbl.GLM$PlanctomycetesPhylum)
### test variance 0.733045
var(tbl.GLM$PlanctomycetesPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PlanctomycetesPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PlanctomycetesPhylum) ## 0.3610687
kurtosis(tbl.GLM$PlanctomycetesPhylum) ## 0.5434426
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPlanctomycetesPhylum <- glm(PlanctomycetesPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                        data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPlanctomycetesPhylum))
qqline(residuals(glm_gaussianPlanctomycetesPhylum), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPlanctomycetesPhylum))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPlanctomycetesPhylum)
# W = 0.98483, p-value = 0.2957

## Test Homoscedasticity
plot(fitted(glm_gaussianPlanctomycetesPhylum), residuals(glm_gaussianPlanctomycetesPhylum), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPlanctomycetesPhylum)
# studentized Breusch-Pagan test
# data:  glm_gaussianPlanctomycetesPhylum
# BP = 25.413, df = 17, p-value = 0.08584

###### 
glm_gaussianPlanctomycetesPhylum <- glm(PlanctomycetesPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                        data = tbl.GLM, family = gaussian)
Anova(glm_gaussianPlanctomycetesPhylum, type = "II") 

glm_gaussianPlanctomycetesPhylum_simplified <- glm(PlanctomycetesPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                     Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                     Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianPlanctomycetesPhylum_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PlanctomycetesPhylum
#                      LR Chisq Df Pr(>Chisq)   
# Inoculants             2.1540  1   0.142196   
# Aggregate             11.1099  2   0.003868 **
# Plant                  3.3702  1   0.066384 . 
# Time                   0.3259  1   0.568076   
# Inoculants:Aggregate   1.6950  2   0.428479   
# Inoculants:Plant       5.2291  1   0.022212 * 
# Aggregate:Plant        1.2814  2   0.526919   
# Inoculants:Time        5.4551  1   0.019512 * 
# Aggregate:Time         8.2527  2   0.016142 * 
# Plant:Time                     0              


#################################  FirmicutesPhylum_model  
table(tbl.GLM$FirmicutesPhylum)
### test mean 3.997059
mean(tbl.GLM$FirmicutesPhylum)
### test variance 1.212197
var(tbl.GLM$FirmicutesPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$FirmicutesPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$FirmicutesPhylum) ## 2.56516
kurtosis(tbl.GLM$FirmicutesPhylum) ## 11.72154
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaFirmicutesPhylum <- glm(FirmicutesPhylum ~ Inoculants * Aggregate * Plant * Time, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaFirmicutesPhylum, type = "II")
glm_gammaFirmicutesPhylum_simplified <- glm(FirmicutesPhylum ~ Inoculants + Aggregate + Plant + Time +
                                              Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                              Plant:Time, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaFirmicutesPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FirmicutesPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              4.270  1  0.0388009 *  
# Aggregate              13.793  2  0.0010111 ** 
# Plant                   6.384  1  0.0115177 *  
# Time                    7.939  1  0.0048370 ** 
# Inoculants:Aggregate    1.465  2  0.4808033    
# Inoculants:Plant       11.716  1  0.0006196 ***
# Aggregate:Plant         0.746  2  0.6886543    
# Inoculants:Time        33.314  1  7.841e-09 ***
# Aggregate:Time          0.543  2  0.7620534    
# Plant:Time                     0               

#################################  ChloroflexiPhylum_model  
table(tbl.GLM$ChloroflexiPhylum)
### test mean 3.80402
mean(tbl.GLM$ChloroflexiPhylum)
### test variance 1.624664
var(tbl.GLM$ChloroflexiPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ChloroflexiPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ChloroflexiPhylum) ## 1.349991
kurtosis(tbl.GLM$ChloroflexiPhylum) ## 2.39876
## The data is right-skewed, 
#  so it is  appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaChloroflexiPhylum <- glm(ChloroflexiPhylum ~ Inoculants * Aggregate * Plant * Time, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaChloroflexiPhylum, type = "II")
glm_gammaChloroflexiPhylum_simplified <- glm(ChloroflexiPhylum ~ Inoculants + Aggregate + Plant + Time +
                                               Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                               Plant:Time, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaChloroflexiPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ChloroflexiPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              3.310  1    0.06887 .  
# Aggregate               5.101  2    0.07804 .  
# Plant                   0.044  1    0.83306    
# Time                   40.491  1  1.975e-10 ***
# Inoculants:Aggregate    3.148  2    0.20724    
# Inoculants:Plant        0.484  1    0.48679    
# Aggregate:Plant         0.897  2    0.63873    
# Inoculants:Time         2.690  1    0.10098    
# Aggregate:Time          1.642  2    0.44002    
# Plant:Time                     0              


#################################  ArmatimonadetesPhylum_model  
table(tbl.GLM$ArmatimonadetesPhylum)
### test mean 2.153039
mean(tbl.GLM$ArmatimonadetesPhylum)
### test variance 0.1582035
var(tbl.GLM$ArmatimonadetesPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ArmatimonadetesPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ArmatimonadetesPhylum) ## 0.8296939
kurtosis(tbl.GLM$ArmatimonadetesPhylum) ## 1.230547
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianArmatimonadetesPhylum <- glm(ArmatimonadetesPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                         data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianArmatimonadetesPhylum))
qqline(residuals(glm_gaussianArmatimonadetesPhylum), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianArmatimonadetesPhylum))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianArmatimonadetesPhylum)
# W = 0.97516, p-value = 0.0514

## Test Homoscedasticity
plot(fitted(glm_gaussianArmatimonadetesPhylum), residuals(glm_gaussianArmatimonadetesPhylum), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianArmatimonadetesPhylum)
# studentized Breusch-Pagan test
# data:  glm_gaussianArmatimonadetesPhylum
# BP = 18.7, df = 17, p-value = 0.346

###### 
glm_gaussianArmatimonadetesPhylum <- glm(ArmatimonadetesPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                         data = tbl.GLM, family = gaussian)
Anova(glm_gaussianArmatimonadetesPhylum, type = "II") 

glm_gaussianArmatimonadetesPhylum_simplified <- glm(ArmatimonadetesPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                      Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                      Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianArmatimonadetesPhylum_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ArmatimonadetesPhylum
#                      LR Chisq Df Pr(>Chisq)  
# Inoculants             0.4623  1    0.49655  
# Aggregate              8.0519  2    0.01785 *
# Plant                  3.2731  1    0.07043 .
# Time                   5.3145  1    0.02115 *
# Inoculants:Aggregate   0.8453  2    0.65532  
# Inoculants:Plant       0.1997  1    0.65500  
# Aggregate:Plant        6.5715  2    0.03741 *
# Inoculants:Time        1.8732  1    0.17111  
# Aggregate:Time         2.5444  2    0.28022  
# Plant:Time                     0             


#################################  BacteroidetesPhylum_model  
table(tbl.GLM$BacteroidetesPhylum)
### test mean 1.298373
mean(tbl.GLM$BacteroidetesPhylum)
### test variance 0.3441602
var(tbl.GLM$BacteroidetesPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist((tbl.GLM$BacteroidetesPhylum), breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$BacteroidetesPhylum) ## 1.517305
kurtosis(tbl.GLM$BacteroidetesPhylum) ## 2.744974å
## The data is  right-skewed, 
#  so it is  appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaBacteroidetesPhylum <- glm(BacteroidetesPhylum ~ Inoculants * Aggregate * Plant * Time, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaBacteroidetesPhylum, type = "II")
glm_gammaBacteroidetesPhylum_simplified <- glm(BacteroidetesPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                 Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                 Plant:Time, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaBacteroidetesPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacteroidetesPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants            11.1959  1  0.0008198 ***
# Aggregate             16.2277  2  0.0002994 ***
# Plant                  4.5395  1  0.0331222 *  
# Time                  14.3858  1  0.0001489 ***
# Inoculants:Aggregate   0.8971  2  0.6385605    
# Inoculants:Plant       4.4726  1  0.0344425 *  
# Aggregate:Plant        1.0590  2  0.5888957    
# Inoculants:Time        5.5299  1  0.0186939 *  
# Aggregate:Time         3.5742  2  0.1674477    
# Plant:Time                     0              


#################################  VerrucomicrobiaPhylum_model  
table(tbl.GLM$VerrucomicrobiaPhylum)
### test mean 2.447392
mean(tbl.GLM$VerrucomicrobiaPhylum)
### test variance 1.177053
var(tbl.GLM$VerrucomicrobiaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$VerrucomicrobiaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$VerrucomicrobiaPhylum) ## 0.8817336
kurtosis(tbl.GLM$VerrucomicrobiaPhylum) ## 1.028995
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianVerrucomicrobiaPhylum <- glm(VerrucomicrobiaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                         data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianVerrucomicrobiaPhylum))
qqline(residuals(glm_gaussianVerrucomicrobiaPhylum), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianVerrucomicrobiaPhylum))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianVerrucomicrobiaPhylum)
# W = 0.94848, p-value = 0.0005729

## Test Homoscedasticity
plot(fitted(glm_gaussianVerrucomicrobiaPhylum), residuals(glm_gaussianVerrucomicrobiaPhylum), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianVerrucomicrobiaPhylum)
# studentized Breusch-Pagan test
# data:  glm_gaussianVerrucomicrobiaPhylum
# BP = 23.994, df = 17, p-value = 0.1196

#### log transformation
tbl.GLM$VerrucomicrobiaPhylum_log <- log(tbl.GLM$VerrucomicrobiaPhylum)
hist(tbl.GLM$VerrucomicrobiaPhylum_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedVerrucomicrobiaPhylum_log <- glm(VerrucomicrobiaPhylum_log ~ Inoculants * Aggregate * Plant * Time, 
                                                data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedVerrucomicrobiaPhylum_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedVerrucomicrobiaPhylum_log)
# W = 0.98415, p-value = 0.2628
## Test Homoscedasticity
plot(fitted(glm_transformedVerrucomicrobiaPhylum_log), residuals(glm_transformedVerrucomicrobiaPhylum_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedVerrucomicrobiaPhylum_log)
# studentized Breusch-Pagan test
# data:  glm_transformedVerrucomicrobiaPhylum_log
# BP = 18.68, df = 17, p-value = 0.3472

#### sqrt transformation
tbl.GLM$VerrucomicrobiaPhylum_sqrt <- sqrt(tbl.GLM$VerrucomicrobiaPhylum)
hist(tbl.GLM$VerrucomicrobiaPhylum_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedVerrucomicrobiaPhylum_sqrt <- glm(VerrucomicrobiaPhylum_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                                 data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedVerrucomicrobiaPhylum_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedVerrucomicrobiaPhylum_sqrt)
# W = 0.97487, p-value = 0.04866
## Test Homoscedasticity
plot(fitted(glm_transformedVerrucomicrobiaPhylum_sqrt), residuals(glm_transformedVerrucomicrobiaPhylum_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedVerrucomicrobiaPhylum_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedVerrucomicrobiaPhylum_sqrt
# BP = 20.085, df = 17, p-value = 0.2699

###### Data transformation (log transformation)
tbl.GLM$VerrucomicrobiaPhylum_log <- log(tbl.GLM$VerrucomicrobiaPhylum)
glm_transformedVerrucomicrobiaPhylum <- glm(VerrucomicrobiaPhylum_log ~ Inoculants * Aggregate * Plant * Time, 
                                            data = tbl.GLM, family = gaussian)
Anova(glm_transformedVerrucomicrobiaPhylum, type = "II") 

glm_transformedVerrucomicrobiaPhylum_simplified <- glm(VerrucomicrobiaPhylum_log ~ Inoculants + Aggregate + Plant + Time +
                                                         Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                         Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_transformedVerrucomicrobiaPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: VerrucomicrobiaPhylum_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              2.766  1   0.096300 .  
# Aggregate               1.413  2   0.493336    
# Plant                   0.157  1   0.692260    
# Time                   80.986  1  < 2.2e-16 ***
# Inoculants:Aggregate    0.852  2   0.653081    
# Inoculants:Plant        3.070  1   0.079749 .  
# Aggregate:Plant         2.349  2   0.309011    
# Inoculants:Time         7.871  1   0.005023 ** 
# Aggregate:Time          1.207  2   0.546892    
# Plant:Time                     0             


#################################  ThaumarchaeotaPhylum_model  
table(tbl.GLM$ThaumarchaeotaPhylum)
### test mean 1.826412
mean(tbl.GLM$ThaumarchaeotaPhylum)
### test variance 2.128622
var(tbl.GLM$ThaumarchaeotaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonThaumarchaeotaPhylum <- glm(ThaumarchaeotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.5037115
deviance(glm_poissonThaumarchaeotaPhylum) / df.residual(glm_poissonThaumarchaeotaPhylum)
### Distribution of response variables
hist(tbl.GLM$ThaumarchaeotaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ThaumarchaeotaPhylum, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiThaumarchaeotaPhylum <- glm(ThaumarchaeotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                     data = tbl.GLM, family = quasipoisson)
### Ensure model fit
plot(fitted(glm_quasiThaumarchaeotaPhylum), residuals(glm_quasiThaumarchaeotaPhylum), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
Anova(glm_quasiThaumarchaeotaPhylum, type = "II")

glm_quasiThaumarchaeotaPhylum_simplified <- glm(ThaumarchaeotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                  Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                  Plant:Time, family = quasipoisson, data = tbl.GLM)
### Results
Anova(glm_quasiThaumarchaeotaPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ThaumarchaeotaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              0.100  1     0.7517    
# Aggregate               3.150  2     0.2070    
# Plant                   0.673  1     0.4120    
# Time                   97.876  1     <2e-16 ***
# Inoculants:Aggregate    0.622  2     0.7326    
# Inoculants:Plant        0.360  1     0.5483    
# Aggregate:Plant         4.363  2     0.1129    
# Inoculants:Time         0.559  1     0.4546    
# Aggregate:Time          0.741  2     0.6903    
# Plant:Time                     0             


#################################  candidate_division_WPS_1Phylum_model  
table(tbl.GLM$candidate_division_WPS_1Phylum)
### test mean 2.470059
mean(tbl.GLM$candidate_division_WPS_1Phylum)
### test variance 1.060323
var(tbl.GLM$candidate_division_WPS_1Phylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$candidate_division_WPS_1Phylum, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$candidate_division_WPS_1Phylum) ## -0.1445163
kurtosis(tbl.GLM$candidate_division_WPS_1Phylum) ## -0.7813791
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussiancandidate_division_WPS_1Phylum <- glm(candidate_division_WPS_1Phylum ~ Inoculants * Aggregate * Plant * Time, 
                                                  data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussiancandidate_division_WPS_1Phylum))
qqline(residuals(glm_gaussiancandidate_division_WPS_1Phylum), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussiancandidate_division_WPS_1Phylum))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussiancandidate_division_WPS_1Phylum)
# W = 0.97551, p-value = 0.05475

## Test Homoscedasticity
plot(fitted(glm_gaussiancandidate_division_WPS_1Phylum), residuals(glm_gaussiancandidate_division_WPS_1Phylum), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussiancandidate_division_WPS_1Phylum)
# studentized Breusch-Pagan test
# data:  glm_gaussiancandidate_division_WPS_1Phylum
# BP = 28.997, df = 17, p-value = 0.03456

#### log transformation
tbl.GLM$candidate_division_WPS_1Phylum_log <- log(tbl.GLM$candidate_division_WPS_1Phylum+1)
hist(tbl.GLM$candidate_division_WPS_1Phylum_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedcandidate_division_WPS_1Phylum_log <- glm(candidate_division_WPS_1Phylum_log ~ Inoculants * Aggregate * Plant * Time, 
                                                         data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedcandidate_division_WPS_1Phylum_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedcandidate_division_WPS_1Phylum_log)
# W = 0.98742, p-value = 0.4521
## Test Homoscedasticity
plot(fitted(glm_transformedcandidate_division_WPS_1Phylum_log), residuals(glm_transformedcandidate_division_WPS_1Phylum_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedcandidate_division_WPS_1Phylum_log)
# studentized Breusch-Pagan test
# data:  glm_transformedcandidate_division_WPS_1Phylum_log
# BP = 19.962, df = 17, p-value = 0.2762

#### sqrt transformation
tbl.GLM$candidate_division_WPS_1Phylum_sqrt <- sqrt(tbl.GLM$candidate_division_WPS_1Phylum)
hist(tbl.GLM$candidate_division_WPS_1Phylum_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedcandidate_division_WPS_1Phylum_sqrt <- glm(candidate_division_WPS_1Phylum_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                                          data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedcandidate_division_WPS_1Phylum_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedcandidate_division_WPS_1Phylum_sqrt)
# W = 0.98671, p-value = 0.404
## Test Homoscedasticity
plot(fitted(glm_transformedcandidate_division_WPS_1Phylum_sqrt), residuals(glm_transformedcandidate_division_WPS_1Phylum_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedcandidate_division_WPS_1Phylum_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedcandidate_division_WPS_1Phylum_sqrt
# BP = 22.381, df = 17, p-value = 0.1705

# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedcandidate_division_WPS_1Phylum_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -100.7587 
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedcandidate_division_WPS_1Phylum_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -74.73451 
### Choose log

###### Data transformation (log transformation)
tbl.GLM$candidate_division_WPS_1Phylum_log <- log(tbl.GLM$candidate_division_WPS_1Phylum)
glm_transformedcandidate_division_WPS_1Phylum <- glm(candidate_division_WPS_1Phylum_log ~ Inoculants * Aggregate * Plant * Time, 
                                                     data = tbl.GLM, family = gaussian)
Anova(glm_transformedcandidate_division_WPS_1Phylum, type = "II") 

glm_transformedcandidate_division_WPS_1Phylum_simplified <- glm(candidate_division_WPS_1Phylum_log ~ Inoculants + Aggregate + Plant + Time +
                                                                  Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                                  Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_transformedcandidate_division_WPS_1Phylum_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: candidate_division_WPS_1Phylum_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               0.36  1  0.5461348    
# Aggregate               15.48  2  0.0004343 ***
# Plant                    5.68  1  0.0171758 *  
# Time                   449.51  1  < 2.2e-16 ***
# Inoculants:Aggregate     8.74  2  0.0126549 *  
# Inoculants:Plant         5.79  1  0.0161477 *  
# Aggregate:Plant          0.96  2  0.6178160    
# Inoculants:Time          9.62  1  0.0019268 ** 
# Aggregate:Time          46.63  2  7.487e-11 ***
# Plant:Time                     0               


#################################  AlphaproteobacteriaClass_model  
table(tbl.GLM$AlphaproteobacteriaClass)
### test mean 13.89951
mean(tbl.GLM$AlphaproteobacteriaClass)
### test variance 3.312005
var(tbl.GLM$AlphaproteobacteriaClass)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$AlphaproteobacteriaClass, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$AlphaproteobacteriaClass) ## 0.05199065
kurtosis(tbl.GLM$AlphaproteobacteriaClass) ## -0.3749388
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianAlphaproteobacteriaClass <- glm(AlphaproteobacteriaClass ~ Inoculants * Aggregate * Plant * Time, 
                                            data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianAlphaproteobacteriaClass))
qqline(residuals(glm_gaussianAlphaproteobacteriaClass), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianAlphaproteobacteriaClass))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianAlphaproteobacteriaClass)
# W = 0.97219, p-value = 0.0299

## Test Homoscedasticity
plot(fitted(glm_gaussianAlphaproteobacteriaClass), residuals(glm_gaussianAlphaproteobacteriaClass), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianAlphaproteobacteriaClass)
# studentized Breusch-Pagan test
# data:  glm_gaussianAlphaproteobacteriaClass
# BP = 37.136, df = 17, p-value = 0.003226

#### log transformation
tbl.GLM$AlphaproteobacteriaClass_log <- log(tbl.GLM$AlphaproteobacteriaClass)
hist(tbl.GLM$AlphaproteobacteriaClass_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedAlphaproteobacteriaClass_log <- glm(AlphaproteobacteriaClass_log ~ Inoculants * Aggregate * Plant * Time, 
                                                   data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedAlphaproteobacteriaClass_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedAlphaproteobacteriaClass_log)
# W = 0.9674, p-value = 0.0127
## Test Homoscedasticity
plot(fitted(glm_transformedAlphaproteobacteriaClass_log), residuals(glm_transformedAlphaproteobacteriaClass_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedAlphaproteobacteriaClass_log)
# studentized Breusch-Pagan test
# data:  glm_transformedAlphaproteobacteriaClass_log
# BP = 42.564, df = 17, p-value = 0.0005557

#### sqrt transformation
tbl.GLM$AlphaproteobacteriaClass_sqrt <- sqrt(tbl.GLM$AlphaproteobacteriaClass)
hist(tbl.GLM$AlphaproteobacteriaClass_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedAlphaproteobacteriaClass_sqrt <- glm(AlphaproteobacteriaClass_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                                    data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedAlphaproteobacteriaClass_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedAlphaproteobacteriaClass_sqrt)
# W = 0.97132, p-value = 0.0255
## Test Homoscedasticity
plot(fitted(glm_transformedAlphaproteobacteriaClass_sqrt), residuals(glm_transformedAlphaproteobacteriaClass_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedAlphaproteobacteriaClass_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedAlphaproteobacteriaClass_sqrt
# BP = 39.831, df = 17, p-value = 0.001367

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianAlphaproteobacteriaClass)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 377.9999 
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedAlphaproteobacteriaClass_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -151.7439 
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedAlphaproteobacteriaClass_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -29.0434 
### Choose log

###### Data transformation (log transformation)
tbl.GLM$AlphaproteobacteriaClass_log <- log(tbl.GLM$AlphaproteobacteriaClass)
glm_transformedAlphaproteobacteriaClass <- glm(AlphaproteobacteriaClass_log ~ Inoculants * Aggregate * Plant * Time, 
                                               data = tbl.GLM, family = gaussian)
Anova(glm_transformedAlphaproteobacteriaClass, type = "II") 

glm_transformedAlphaproteobacteriaClass_simplified <- glm(AlphaproteobacteriaClass_log ~ Inoculants + Aggregate + Plant + Time +
                                                            Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                            Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_transformedAlphaproteobacteriaClass_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: AlphaproteobacteriaClass_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              6.181  1    0.01291 *  
# Aggregate              46.645  2  7.434e-11 ***
# Plant                   2.276  1    0.13136    
# Time                    0.166  1    0.68350    
# Inoculants:Aggregate    1.770  2    0.41276    
# Inoculants:Plant        0.102  1    0.74946    
# Aggregate:Plant         0.095  2    0.95379    
# Inoculants:Time         2.742  1    0.09773 .  
# Aggregate:Time          1.391  2    0.49889    
# Plant:Time                     0             


#################################  GammaproteobacteriaClass_model  
table(tbl.GLM$GammaproteobacteriaClass)
### test mean 3.751863
mean(tbl.GLM$GammaproteobacteriaClass)
### test variance 0.343849
var(tbl.GLM$GammaproteobacteriaClass)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$GammaproteobacteriaClass, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$GammaproteobacteriaClass) ## 0.6619236
kurtosis(tbl.GLM$GammaproteobacteriaClass) ## 0.7438434
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianGammaproteobacteriaClass <- glm(GammaproteobacteriaClass ~ Inoculants * Aggregate * Plant * Time, 
                                            data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianGammaproteobacteriaClass))
qqline(residuals(glm_gaussianGammaproteobacteriaClass), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianGammaproteobacteriaClass))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianGammaproteobacteriaClass)
# W = 0.97777, p-value = 0.08292

## Test Homoscedasticity
plot(fitted(glm_gaussianGammaproteobacteriaClass), residuals(glm_gaussianGammaproteobacteriaClass), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianGammaproteobacteriaClass)
# studentized Breusch-Pagan test
# data:  glm_gaussianGammaproteobacteriaClass
# BP = 14.788, df = 17, p-value = 0.6108

###### 
glm_gaussianGammaproteobacteriaClass <- glm(GammaproteobacteriaClass ~ Inoculants * Aggregate * Plant * Time, 
                                            data = tbl.GLM, family = gaussian)
Anova(glm_gaussianGammaproteobacteriaClass, type = "II") 

glm_gaussianGammaproteobacteriaClass_simplified <- glm(GammaproteobacteriaClass ~ Inoculants + Aggregate + Plant + Time +
                                                         Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                         Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianGammaproteobacteriaClass_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: GammaproteobacteriaClass
#                      LR Chisq Df Pr(>Chisq)   
# Inoculants             5.2151  1   0.022392 * 
# Aggregate              2.7282  2   0.255605   
# Plant                  7.9817  1   0.004725 **
# Time                   0.0061  1   0.937817   
# Inoculants:Aggregate   0.2057  2   0.902248   
# Inoculants:Plant       3.6749  1   0.055236 . 
# Aggregate:Plant        0.2495  2   0.882737   
# Inoculants:Time        0.0037  1   0.951472   
# Aggregate:Time         2.3759  2   0.304850   
# Plant:Time                     0           


#################################  BetaproteobacteriaClass_model  
table(tbl.GLM$BetaproteobacteriaClass)
### test mean 5.524118
mean(tbl.GLM$BetaproteobacteriaClass)
### test variance 1.090338
var(tbl.GLM$BetaproteobacteriaClass)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$BetaproteobacteriaClass, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$BetaproteobacteriaClass) ## 1.090338
kurtosis(tbl.GLM$BetaproteobacteriaClass) ## 2.489005
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianBetaproteobacteriaClass <- glm(BetaproteobacteriaClass ~ Inoculants * Aggregate * Plant * Time, 
                                           data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianBetaproteobacteriaClass))
qqline(residuals(glm_gaussianBetaproteobacteriaClass), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianBetaproteobacteriaClass))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianBetaproteobacteriaClass)
# W = 0.95541, p-value = 0.00169

## Test Homoscedasticity
plot(fitted(glm_gaussianBetaproteobacteriaClass), residuals(glm_gaussianBetaproteobacteriaClass), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianBetaproteobacteriaClass)
# studentized Breusch-Pagan test
# data:  glm_gaussianBetaproteobacteriaClass
# BP = 21.849, df = 17, p-value = 0.1906

#### log transformation
tbl.GLM$BetaproteobacteriaClass_log <- log(tbl.GLM$BetaproteobacteriaClass)
hist(tbl.GLM$BetaproteobacteriaClass_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedBetaproteobacteriaClass_log <- glm(BetaproteobacteriaClass_log ~ Inoculants * Aggregate * Plant * Time, 
                                                  data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedBetaproteobacteriaClass_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedBetaproteobacteriaClass_log)
# W = 0.97399, p-value = 0.04147
## Test Homoscedasticity
plot(fitted(glm_transformedBetaproteobacteriaClass_log), residuals(glm_transformedBetaproteobacteriaClass_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedBetaproteobacteriaClass_log)
# studentized Breusch-Pagan test
# data:  glm_transformedBetaproteobacteriaClass_log
# BP = 13.08, df = 17, p-value = 0.7308

#### sqrt transformation
tbl.GLM$BetaproteobacteriaClass_sqrt <- sqrt(tbl.GLM$BetaproteobacteriaClass)
hist(tbl.GLM$BetaproteobacteriaClass_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedBetaproteobacteriaClass_sqrt <- glm(BetaproteobacteriaClass_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                                   data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedBetaproteobacteriaClass_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedBetaproteobacteriaClass_sqrt)
# W = 0.97807, p-value = 0.08772
## Test Homoscedasticity
plot(fitted(glm_transformedBetaproteobacteriaClass_sqrt), residuals(glm_transformedBetaproteobacteriaClass_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedBetaproteobacteriaClass_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedBetaproteobacteriaClass_sqrt
# BP = 16.959, df = 17, p-value = 0.4571

###### sqrt transformation
tbl.GLM$BetaproteobacteriaClass_sqrt <- sqrt(tbl.GLM$BetaproteobacteriaClass)
glm_transformedBetaproteobacteriaClass_sqrt <- glm(BetaproteobacteriaClass_sqrt ~ Inoculants * Aggregate * Plant * Time, 
                                                   data = tbl.GLM, family = gaussian)
Anova(glm_transformedBetaproteobacteriaClass_sqrt, type = "II") 

glm_transformedBetaproteobacteriaClass_sqrt_simplified <- glm(BetaproteobacteriaClass_sqrt ~ Inoculants + Aggregate + Plant + Time +
                                                                Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                                Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_transformedBetaproteobacteriaClass_sqrt_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BetaproteobacteriaClass_sqrt
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants            11.6378  1  0.0006462 ***
# Aggregate              5.1611  2  0.0757309 .  
# Plant                  0.1827  1  0.6691023    
# Time                   0.2230  1  0.6367623    
# Inoculants:Aggregate   0.0425  2  0.9789716    
# Inoculants:Plant       4.1240  1  0.0422797 *  
# Aggregate:Plant        0.9220  2  0.6306443    
# Inoculants:Time        0.7071  1  0.4004061    
# Aggregate:Time         0.4333  2  0.8051966    
# Plant:Time                     0             


#################################  PlanctomycetaciaClass_model  
table(tbl.GLM$PlanctomycetaciaClass)
### test mean 5.048137
mean(tbl.GLM$PlanctomycetaciaClass)
### test variance 0.73135
var(tbl.GLM$PlanctomycetaciaClass)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PlanctomycetaciaClass, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PlanctomycetaciaClass) ## 0.3660289
kurtosis(tbl.GLM$PlanctomycetaciaClass) ## 0.5531545
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPlanctomycetaciaClass <- glm(PlanctomycetaciaClass ~ Inoculants * Aggregate * Plant * Time, 
                                         data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPlanctomycetaciaClass))
qqline(residuals(glm_gaussianPlanctomycetaciaClass), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPlanctomycetaciaClass))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPlanctomycetaciaClass)
# W = 0.98485, p-value = 0.2968

## Test Homoscedasticity
plot(fitted(glm_gaussianPlanctomycetaciaClass), residuals(glm_gaussianPlanctomycetaciaClass), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPlanctomycetaciaClass)
# studentized Breusch-Pagan test
# data:  glm_gaussianPlanctomycetaciaClass
# BP = 25.249, df = 17, p-value = 0.08927

###### 
glm_gaussianPlanctomycetaciaClass <- glm(PlanctomycetaciaClass ~ Inoculants * Aggregate * Plant * Time, 
                                         data = tbl.GLM, family = gaussian)
Anova(glm_gaussianPlanctomycetaciaClass, type = "II") 

glm_gaussianPlanctomycetaciaClass_simplified <- glm(PlanctomycetaciaClass ~ Inoculants + Aggregate + Plant + Time +
                                                      Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                      Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianPlanctomycetaciaClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: PlanctomycetaciaClass
#                      LR Chisq Df Pr(>Chisq)   
# Inoculants             2.1969  1   0.138288   
# Aggregate             10.9982  2   0.004091 **
# Plant                  3.4117  1   0.064736 . 
# Time                   0.3283  1   0.566654   
# Inoculants:Aggregate   1.7106  2   0.425147   
# Inoculants:Plant       5.2548  1   0.021886 * 
# Aggregate:Plant        1.2673  2   0.530640   
# Inoculants:Time        5.4850  1   0.019180 * 
# Aggregate:Time         8.2711  2   0.015994 * 
# Plant:Time                     0              


#################################  Acidobacteria_Gp1Class_model  
table(tbl.GLM$Acidobacteria_Gp1Class)
### test mean 20.09833
mean(tbl.GLM$Acidobacteria_Gp1Class)
### test variance 18.04644
var(tbl.GLM$Acidobacteria_Gp1Class)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonAcidobacteria_Gp1Class <- glm(Acidobacteria_Gp1Class ~ Inoculants * Aggregate * Plant * Time, 
                                         data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.2315212
deviance(glm_poissonAcidobacteria_Gp1Class) / df.residual(glm_poissonAcidobacteria_Gp1Class)
### Distribution of response variables
hist(tbl.GLM$Acidobacteria_Gp1Class, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$Acidobacteria_Gp1Class) ## -0.4544026
kurtosis(tbl.GLM$Acidobacteria_Gp1Class) ## -0.2934624
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianAcidobacteria_Gp1Class <- glm(Acidobacteria_Gp1Class ~ Inoculants * Aggregate * Plant * Time, 
                                          data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianAcidobacteria_Gp1Class))
qqline(residuals(glm_gaussianAcidobacteria_Gp1Class), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianAcidobacteria_Gp1Class))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianAcidobacteria_Gp1Class)
# W = 0.99095, p-value = 0.7295

## Test Homoscedasticity
plot(fitted(glm_gaussianAcidobacteria_Gp1Class), residuals(glm_gaussianAcidobacteria_Gp1Class), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianAcidobacteria_Gp1Class)
# studentized Breusch-Pagan test
# data:  glm_gaussianAcidobacteria_Gp1Class
# BP = 18.169, df = 17, p-value = 0.3783
###### 
glm_gaussianAcidobacteria_Gp1Class <- glm(Acidobacteria_Gp1Class ~ Inoculants * Aggregate * Plant * Time, 
                                          data = tbl.GLM, family = gaussian)
Anova(glm_gaussianAcidobacteria_Gp1Class, type = "II") 

glm_gaussianAcidobacteria_Gp1Class_simplified <- glm(Acidobacteria_Gp1Class ~ Inoculants + Aggregate + Plant + Time +
                                                       Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                       Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianAcidobacteria_Gp1Class_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: Acidobacteria_Gp1Class
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             14.157  1  0.0001682 ***
# Aggregate               6.933  2  0.0312318 *  
# Plant                   1.980  1  0.1594083    
# Time                  198.717  1  < 2.2e-16 ***
# Inoculants:Aggregate    1.022  2  0.5999786    
# Inoculants:Plant        0.015  1  0.9011943    
# Aggregate:Plant         0.651  2  0.7222930    
# Inoculants:Time         3.724  1  0.0536428 .  
# Aggregate:Time          0.100  2  0.9512392    
# Plant:Time                     0              


#################################  ActinobacteriaClass_model  
table(tbl.GLM$ActinobacteriaClass)
### test mean 2.298627
mean(tbl.GLM$ActinobacteriaClass)
### test variance 0.439404
var(tbl.GLM$ActinobacteriaClass)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ActinobacteriaClass, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ActinobacteriaClass) ## 0.6288197
kurtosis(tbl.GLM$ActinobacteriaClass) ## 0.2152862
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianActinobacteriaClass <- glm(ActinobacteriaClass ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianActinobacteriaClass))
qqline(residuals(glm_gaussianActinobacteriaClass), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianActinobacteriaClass))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianActinobacteriaClass)
# W = 0.97954, p-value = 0.1149

## Test Homoscedasticity
plot(fitted(glm_gaussianActinobacteriaClass), residuals(glm_gaussianActinobacteriaClass), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianActinobacteriaClass)
# studentized Breusch-Pagan test
# data:  glm_gaussianActinobacteriaClass
# BP = 25.837, df = 17, p-value = 0.07747

###### 
glm_gaussianActinobacteriaClass <- glm(ActinobacteriaClass ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = gaussian)
Anova(glm_gaussianActinobacteriaClass, type = "II") 

glm_gaussianActinobacteriaClass_simplified <- glm(ActinobacteriaClass ~ Inoculants + Aggregate + Plant + Time +
                                                    Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                    Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianActinobacteriaClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ActinobacteriaClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             6.1993  1    0.01278 *  
# Aggregate             24.3573  2  5.139e-06 ***
# Plant                 30.5437  1  3.264e-08 ***
# Time                   3.4045  1    0.06502 .  
# Inoculants:Aggregate   0.1050  2    0.94883    
# Inoculants:Plant       0.4550  1    0.49999    
# Aggregate:Plant        0.7075  2    0.70204    
# Inoculants:Time        0.1259  1    0.72277    
# Aggregate:Time         0.0620  2    0.96948    
# Plant:Time                     0             


#################################  Acidobacteria_Gp3Class_model  
table(tbl.GLM$Acidobacteria_Gp3Class)
### test mean 4.698627
mean(tbl.GLM$Acidobacteria_Gp3Class)
### test variance 0.919608
var(tbl.GLM$Acidobacteria_Gp3Class)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$Acidobacteria_Gp3Class, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$Acidobacteria_Gp3Class) ## 0.7735834
kurtosis(tbl.GLM$Acidobacteria_Gp3Class) ## 1.607416
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianAcidobacteria_Gp3Class <- glm(Acidobacteria_Gp3Class ~ Inoculants * Aggregate * Plant * Time, 
                                          data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianAcidobacteria_Gp3Class))
qqline(residuals(glm_gaussianAcidobacteria_Gp3Class), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianAcidobacteria_Gp3Class))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianAcidobacteria_Gp3Class)
# W = 0.99233, p-value = 0.8358

## Test Homoscedasticity
plot(fitted(glm_gaussianAcidobacteria_Gp3Class), residuals(glm_gaussianAcidobacteria_Gp3Class), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianAcidobacteria_Gp3Class)
# studentized Breusch-Pagan test
# data:  glm_gaussianAcidobacteria_Gp3Class
# BP = 11.746, df = 17, p-value = 0.8153

###### 
glm_gaussianAcidobacteria_Gp3Class <- glm(Acidobacteria_Gp3Class ~ Inoculants * Aggregate * Plant * Time, 
                                          data = tbl.GLM, family = gaussian)
Anova(glm_gaussianAcidobacteria_Gp3Class, type = "II") 

glm_gaussianAcidobacteria_Gp3Class_simplified <- glm(Acidobacteria_Gp3Class ~ Inoculants + Aggregate + Plant + Time +
                                                       Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                       Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianAcidobacteria_Gp3Class_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: Acidobacteria_Gp3Class
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants            22.8097  1  1.789e-06 ***
# Aggregate             24.9512  2  3.819e-06 ***
# Plant                  0.2656  1   0.606292    
# Time                  31.1688  1  2.365e-08 ***
# Inoculants:Aggregate   0.3865  2   0.824270    
# Inoculants:Plant       3.0203  1   0.082231 .  
# Aggregate:Plant        0.0136  2   0.993230    
# Inoculants:Time        1.4241  1   0.232723    
# Aggregate:Time        11.2473  2   0.003611 ** 
# Plant:Time                     0              


#################################  BacilliClass_model  
table(tbl.GLM$BacilliClass)
### test mean 0.8596863
mean(tbl.GLM$BacilliClass)
### test variance 0.6737007
var(tbl.GLM$BacilliClass)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$BacilliClass, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$BacilliClass) ## 6.217204
kurtosis(tbl.GLM$BacilliClass) ## 47.31919
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaBacilliClass <- glm(BacilliClass ~ Inoculants * Aggregate * Plant * Time, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaBacilliClass, type = "II")
glm_gammaBacilliClass_simplified <- glm(BacilliClass ~ Inoculants + Aggregate + Plant + Time +
                                          Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                          Plant:Time, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaBacilliClass_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: BacilliClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             23.309  1  1.380e-06 ***
# Aggregate               7.277  2     0.0263 *  
# Plant                  54.851  1  1.300e-13 ***
# Time                   30.079  1  4.148e-08 ***
# Inoculants:Aggregate    4.037  2     0.1328    
# Inoculants:Plant       34.630  1  3.987e-09 ***
# Aggregate:Plant         4.301  2     0.1164    
# Inoculants:Time        31.867  1  1.651e-08 ***
# Aggregate:Time          3.244  2     0.1975    
# Plant:Time                     0              



#################################  Acidobacteria_Gp2Class_model  
table(tbl.GLM$Acidobacteria_Gp2Class)
### test mean 10.7652
mean(tbl.GLM$Acidobacteria_Gp2Class)
### test variance 7.431368
var(tbl.GLM$Acidobacteria_Gp2Class)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$Acidobacteria_Gp2Class, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$Acidobacteria_Gp2Class) ## 0.05199065
kurtosis(tbl.GLM$Acidobacteria_Gp2Class) ## -0.3749388
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianAcidobacteria_Gp2Class <- glm(Acidobacteria_Gp2Class ~ Inoculants * Aggregate * Plant * Time, 
                                          data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianAcidobacteria_Gp2Class))
qqline(residuals(glm_gaussianAcidobacteria_Gp2Class), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianAcidobacteria_Gp2Class))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianAcidobacteria_Gp2Class)
# W = 0.99256, p-value = 0.852

## Test Homoscedasticity
plot(fitted(glm_gaussianAcidobacteria_Gp2Class), residuals(glm_gaussianAcidobacteria_Gp2Class), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianAcidobacteria_Gp2Class)
# studentized Breusch-Pagan test
# data:  glm_gaussianAcidobacteria_Gp2Class
# BP = 24.978, df = 17, p-value = 0.0952

###### 
glm_gaussianAcidobacteria_Gp2Class <- glm(Acidobacteria_Gp2Class ~ Inoculants * Aggregate * Plant * Time, 
                                          data = tbl.GLM, family = gaussian)
Anova(glm_gaussianAcidobacteria_Gp2Class, type = "II") 

glm_gaussianAcidobacteria_Gp2Class_simplified <- glm(Acidobacteria_Gp2Class ~ Inoculants + Aggregate + Plant + Time +
                                                       Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                       Plant:Time, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianAcidobacteria_Gp2Class_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: Acidobacteria_Gp2Class
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants            13.8246  1  0.0002007 ***
# Aggregate              9.2564  2  0.0097724 ** 
# Plant                  0.1404  1  0.7078901    
# Time                  10.1370  1  0.0014532 ** 
# Inoculants:Aggregate   0.1432  2  0.9309143    
# Inoculants:Plant       5.5679  1  0.0182922 *  
# Aggregate:Plant        0.7731  2  0.6794099    
# Inoculants:Time        0.3434  1  0.5578935    
# Aggregate:Time         3.1203  2  0.2101051    
# Plant:Time                     0            

tbl.GLM <- read.csv("FungalGLM.csv")
#################################  AscomycotaPhylum_model  
table(tbl.GLM$AscomycotaPhylum)
### test mean 53.08824
mean(tbl.GLM$AscomycotaPhylum)
### test variance 98.13511
var(tbl.GLM$AscomycotaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonAscomycotaPhylum <- glm(AscomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                   data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.629036
deviance(glm_poissonAscomycotaPhylum) / df.residual(glm_poissonAscomycotaPhylum)
### Distribution of response variables
hist(tbl.GLM$AscomycotaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$AscomycotaPhylum, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiAscomycotaPhylum <- glm(AscomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                 data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiAscomycotaPhylum), residuals(glm_quasiAscomycotaPhylum), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiAscomycotaPhylum, type = "II")

glm_quasiAscomycotaPhylum_simplified <- glm(AscomycotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                              Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                              Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiAscomycotaPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: AscomycotaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants            11.4793  1  0.0007038 ***
# Aggregate              2.9525  2  0.2284922    
# Plant                  1.2246  1  0.2684579    
# Time                   8.3829  1  0.0037876 ** 
# Inoculants:Aggregate   1.2321  2  0.5400604    
# Inoculants:Plant      11.8524  1  0.0005759 ***
# Aggregate:Plant        0.0253  2  0.9874469    
# Inoculants:Time        5.5751  1  0.0182173 *  
# Aggregate:Time         0.3571  2  0.8364812    
# Plant:Time                     0              


#################################  BasidiomycotaPhylum_model  
table(tbl.GLM$BasidiomycotaPhylum)
### test mean 28.55627
mean(tbl.GLM$BasidiomycotaPhylum)
### test variance 203.275
var(tbl.GLM$BasidiomycotaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonBasidiomycotaPhylum <- glm(BasidiomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                      data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.629036
deviance(glm_poissonBasidiomycotaPhylum) / df.residual(glm_poissonBasidiomycotaPhylum)
### Distribution of response variables
hist(tbl.GLM$BasidiomycotaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$BasidiomycotaPhylum, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiBasidiomycotaPhylum <- glm(BasidiomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                    data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiBasidiomycotaPhylum), residuals(glm_quasiBasidiomycotaPhylum), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiBasidiomycotaPhylum, type = "II")

glm_quasiBasidiomycotaPhylum_simplified <- glm(BasidiomycotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                 Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                 Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiBasidiomycotaPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BasidiomycotaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              6.213  1    0.01268 *  
# Aggregate               0.809  2    0.66719    
# Plant                   0.052  1    0.82002    
# Time                  110.441  1    < 2e-16 ***
# Inoculants:Aggregate    2.130  2    0.34481    
# Inoculants:Plant        2.542  1    0.11088    
# Aggregate:Plant         0.125  2    0.93950    
# Inoculants:Time         0.170  1    0.68044    
# Aggregate:Time          1.026  2    0.59876    
# Plant:Time                     0              


#################################  MucoromycotaPhylum_model  
table(tbl.GLM$MucoromycotaPhylum)
### test mean 1.517608
mean(tbl.GLM$MucoromycotaPhylum)
### test variance 1.301999
var(tbl.GLM$MucoromycotaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$MucoromycotaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$MucoromycotaPhylum) ## 2.174217
kurtosis(tbl.GLM$MucoromycotaPhylum) ## 8.434694
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaMucoromycotaPhylum <- glm(MucoromycotaPhylum ~ Inoculants * Aggregate * Plant * Time, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaMucoromycotaPhylum, type = "II")
glm_gammaMucoromycotaPhylum_simplified <- glm(MucoromycotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                Plant:Time, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaMucoromycotaPhylum_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: MucoromycotaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              1.501  1    0.22057    
# Aggregate               3.421  2    0.18074    
# Plant                   1.441  1    0.22999    
# Time                  106.959  1    < 2e-16 ***
# Inoculants:Aggregate    0.534  2    0.76567    
# Inoculants:Plant        0.657  1    0.41776    
# Aggregate:Plant         2.562  2    0.27777    
# Inoculants:Time         5.593  1    0.01804 *  
# Aggregate:Time          0.455  2    0.79665    
# Plant:Time                     0             



#################################  GlomeromycotaPhylum_model  
table(tbl.GLM$GlomeromycotaPhylum)
### test mean 0.4625225
mean(tbl.GLM$GlomeromycotaPhylum)
### test variance 0.2344595
var(tbl.GLM$GlomeromycotaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$GlomeromycotaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$GlomeromycotaPhylum) ## 1.974437
kurtosis(tbl.GLM$GlomeromycotaPhylum) ## 4.980416
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaGlomeromycotaPhylum <- glm(GlomeromycotaPhylum ~ Inoculants * Aggregate * Plant * Time, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaGlomeromycotaPhylum, type = "II")
glm_gammaGlomeromycotaPhylum_simplified <- glm(GlomeromycotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                 Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                 Plant:Time, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaGlomeromycotaPhylum_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: GlomeromycotaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              0.156  1  0.6931496    
# Aggregate               4.278  2  0.1178013    
# Plant                 104.273  1  < 2.2e-16 ***
# Time                   13.265  1  0.0002704 ***
# Inoculants:Aggregate    4.402  2  0.1106762    
# Inoculants:Plant        0.124  1  0.7245848    
# Aggregate:Plant         1.892  2  0.3883120    
# Inoculants:Time         1.432  1  0.2314164    
# Aggregate:Time          0.178  2  0.9147919    
# Plant:Time                     0              



#################################  MortierellomycotaPhylum_model  
table(tbl.GLM$MortierellomycotaPhylum)
### test mean 2.6225
mean(tbl.GLM$MortierellomycotaPhylum)
### test variance 7.923284
var(tbl.GLM$MortierellomycotaPhylum)
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaMortierellomycotaPhylum <- glm(MortierellomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaMortierellomycotaPhylum, type = "II")
glm_gammaMortierellomycotaPhylum_simplified <- glm(MortierellomycotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                     Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                     Plant:Time, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaMortierellomycotaPhylum_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: MortierellomycotaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             8.2756  1  0.0040182 ** 
# Aggregate              0.1563  2  0.9248225    
# Plant                 13.6547  1  0.0002197 ***
# Time                   0.5702  1  0.4501636    
# Inoculants:Aggregate   0.1280  2  0.9379913    
# Inoculants:Plant      16.4406  1   5.02e-05 ***
# Aggregate:Plant        0.8204  2  0.6635022    
# Inoculants:Time       11.1963  1  0.0008196 ***
# Aggregate:Time         1.9046  2  0.3858459    
# Plant:Time                     0               



#################################  KickxellomycotaPhylum_model  
table(tbl.GLM$KickxellomycotaPhylum)
### test mean 0.08744314
mean(tbl.GLM$KickxellomycotaPhylum)
### test variance 0.006629734
var(tbl.GLM$KickxellomycotaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$KickxellomycotaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")

### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$KickxellomycotaPhylum) ## 3.268177
kurtosis(tbl.GLM$KickxellomycotaPhylum) ## 13.68641
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaKickxellomycotaPhylum <- glm(KickxellomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaKickxellomycotaPhylum, type = "II")
glm_gammaKickxellomycotaPhylum_simplified <- glm(KickxellomycotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                   Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                   Plant:Time, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaKickxellomycotaPhylum_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: KickxellomycotaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             0.5015  1     0.4788    
# Aggregate              1.0841  2     0.5816    
# Plant                 18.0418  1  2.161e-05 ***
# Time                   2.4995  1     0.1139    
# Inoculants:Aggregate   3.1415  2     0.2079    
# Inoculants:Plant       0.2815  1     0.5957    
# Aggregate:Plant        1.8182  2     0.4029    
# Inoculants:Time        0.5194  1     0.4711    
# Aggregate:Time         0.3315  2     0.8473    
# Plant:Time                     0             



#################################  ChytridiomycotaPhylum_model  
table(tbl.GLM$ChytridiomycotaPhylum)
### test mean 0.6263982
mean(tbl.GLM$ChytridiomycotaPhylum)
### test variance 3.232445
var(tbl.GLM$ChytridiomycotaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonChytridiomycotaPhylum <- glm(ChytridiomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                        data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.6921894
deviance(glm_poissonChytridiomycotaPhylum) / df.residual(glm_poissonChytridiomycotaPhylum)
### Distribution of response variables
hist(tbl.GLM$ChytridiomycotaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ChytridiomycotaPhylum, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiChytridiomycotaPhylum <- glm(ChytridiomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                      data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiChytridiomycotaPhylum), residuals(glm_quasiChytridiomycotaPhylum), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiChytridiomycotaPhylum, type = "II")

glm_quasiChytridiomycotaPhylum_simplified <- glm(ChytridiomycotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                   Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                   Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiChytridiomycotaPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ChytridiomycotaPhylum
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             20.850  1  4.968e-06 ***
# Aggregate              19.226  2  6.686e-05 ***
# Plant                   0.339  1     0.5604    
# Time                   70.698  1  < 2.2e-16 ***
# Inoculants:Aggregate    0.602  2     0.7401    
# Inoculants:Plant        0.000  1     0.9873    
# Aggregate:Plant         0.455  2     0.7963    
# Inoculants:Time         1.806  1     0.1789    
# Aggregate:Time          2.819  2     0.2442    
# Plant:Time                     0               



#################################  RozellomycotaPhylum_model  
table(tbl.GLM$RozellomycotaPhylum)
### test mean 1.370141
mean(tbl.GLM$RozellomycotaPhylum)
### test variance 26.47138
var(tbl.GLM$RozellomycotaPhylum)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonRozellomycotaPhylum <- glm(RozellomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                      data = tbl.GLM, family = poisson)
### The discrete factor was calculated  2.933197
deviance(glm_poissonRozellomycotaPhylum) / df.residual(glm_poissonRozellomycotaPhylum)
### Distribution of response variables
hist(tbl.GLM$RozellomycotaPhylum, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$RozellomycotaPhylum, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiRozellomycotaPhylum <- glm(RozellomycotaPhylum ~ Inoculants * Aggregate * Plant * Time, 
                                    data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiRozellomycotaPhylum), residuals(glm_quasiRozellomycotaPhylum), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiRozellomycotaPhylum, type = "II")

glm_quasiRozellomycotaPhylum_simplified <- glm(RozellomycotaPhylum ~ Inoculants + Aggregate + Plant + Time +
                                                 Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                 Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiRozellomycotaPhylum_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: RozellomycotaPhylum
#                      LR Chisq Df Pr(>Chisq)   
# Inoculants             5.9485  1   0.014729 * 
# Aggregate              6.4831  2   0.039104 * 
# Plant                  7.0383  1   0.007978 **
# Time                   0.2476  1   0.618740   
# Inoculants:Aggregate   1.6983  2   0.427784   
# Inoculants:Plant       0.0462  1   0.829831   
# Aggregate:Plant        1.7618  2   0.414405   
# Inoculants:Time        2.2204  1   0.136200   
# Aggregate:Time         9.9714  2   0.006835 **
# Plant:Time                     0              



#################################  SordariomycetesClass_model  
table(tbl.GLM$SordariomycetesClass)
### test mean 18.98402
mean(tbl.GLM$SordariomycetesClass)
### test variance 79.68218
var(tbl.GLM$SordariomycetesClass)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonSordariomycetesClass <- glm(SordariomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = poisson)
### The discrete factor was calculated  2.23519
deviance(glm_poissonSordariomycetesClass) / df.residual(glm_poissonSordariomycetesClass)
### Distribution of response variables
hist(tbl.GLM$SordariomycetesClass, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$SordariomycetesClass, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiSordariomycetesClass <- glm(SordariomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                     data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiSordariomycetesClass), residuals(glm_quasiSordariomycetesClass), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiSordariomycetesClass, type = "II")

glm_quasiSordariomycetesClass_simplified <- glm(SordariomycetesClass ~ Inoculants + Aggregate + Plant + Time +
                                                  Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                  Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiSordariomycetesClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: SordariomycetesClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants            23.8707  1  1.030e-06 ***
# Aggregate              4.6573  2    0.09743 .  
# Plant                  2.2113  1    0.13700    
# Time                   0.0190  1    0.89040    
# Inoculants:Aggregate   3.0035  2    0.22274    
# Inoculants:Plant      27.3828  1  1.669e-07 ***
# Aggregate:Plant        0.7602  2    0.68380    
# Inoculants:Time       29.0704  1  6.980e-08 ***
# Aggregate:Time         7.4970  2    0.02355 *  
# Plant:Time                     0               



#################################  AgaricomycetesClass_model  
table(tbl.GLM$AgaricomycetesClass)
### test mean 18.025
mean(tbl.GLM$AgaricomycetesClass)
### test variance 131.4566
var(tbl.GLM$AgaricomycetesClass)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonAgaricomycetesClass <- glm(AgaricomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                      data = tbl.GLM, family = poisson)
### The discrete factor was calculated  4.9159
deviance(glm_poissonAgaricomycetesClass) / df.residual(glm_poissonAgaricomycetesClass)
### Distribution of response variables
hist(tbl.GLM$AgaricomycetesClass, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$AgaricomycetesClass, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiAgaricomycetesClass <- glm(AgaricomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                    data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiAgaricomycetesClass), residuals(glm_quasiAgaricomycetesClass), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiAgaricomycetesClass, type = "II")

glm_quasiAgaricomycetesClass_simplified <- glm(AgaricomycetesClass ~ Inoculants + Aggregate + Plant + Time +
                                                 Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                 Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiAgaricomycetesClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: AgaricomycetesClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              4.204  1    0.04034 *  
# Aggregate               2.584  2    0.27467    
# Plant                   7.538  1    0.00604 ** 
# Time                   37.042  1  1.156e-09 ***
# Inoculants:Aggregate    1.780  2    0.41056    
# Inoculants:Plant        0.192  1    0.66088    
# Aggregate:Plant         0.683  2    0.71071    
# Inoculants:Time         0.472  1    0.49198    
# Aggregate:Time          0.384  2    0.82527    
# Plant:Time                     0               



#################################  EurotiomycetesClass_model  
table(tbl.GLM$EurotiomycetesClass)
### test mean 18.45029
mean(tbl.GLM$EurotiomycetesClass)
### test variance 74.42203
var(tbl.GLM$EurotiomycetesClass)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonEurotiomycetesClass <- glm(EurotiomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                      data = tbl.GLM, family = poisson)
### The discrete factor was calculated  2.147091
deviance(glm_poissonEurotiomycetesClass) / df.residual(glm_poissonEurotiomycetesClass)
### Distribution of response variables
hist(tbl.GLM$EurotiomycetesClass, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$EurotiomycetesClass, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiEurotiomycetesClass <- glm(EurotiomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                    data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiEurotiomycetesClass), residuals(glm_quasiEurotiomycetesClass), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiEurotiomycetesClass, type = "II")

glm_quasiEurotiomycetesClass_simplified <- glm(EurotiomycetesClass ~ Inoculants + Aggregate + Plant + Time +
                                                 Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                 Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiEurotiomycetesClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: EurotiomycetesClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              2.546  1  0.1105551    
# Aggregate               1.768  2  0.4131421    
# Plant                  12.674  1  0.0003708 ***
# Time                   39.859  1   2.73e-10 ***
# Inoculants:Aggregate    0.367  2  0.8323195    
# Inoculants:Plant        0.029  1  0.8657161    
# Aggregate:Plant         2.006  2  0.3668264    
# Inoculants:Time         2.302  1  0.1292310    
# Aggregate:Time          0.416  2  0.8120626    
# Plant:Time                     0               



#################################  LecanoromycetesClass_model  
table(tbl.GLM$LecanoromycetesClass)
### test mean 2.436431
mean(tbl.GLM$LecanoromycetesClass)
### test variance 2.75732
var(tbl.GLM$LecanoromycetesClass)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonLecanoromycetesClass <- glm(LecanoromycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.8819533
deviance(glm_poissonLecanoromycetesClass) / df.residual(glm_poissonLecanoromycetesClass)
### Distribution of response variables
hist(tbl.GLM$LecanoromycetesClass, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$LecanoromycetesClass, "pois")
plot(fit)

##################### Quasi-Poisson
glm_poissonLecanoromycetesClass <- glm(LecanoromycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = poisson)

### Ensure model fit
plot(fitted(glm_poissonLecanoromycetesClass), residuals(glm_poissonLecanoromycetesClass), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
Anova(glm_poissonLecanoromycetesClass, type = "II")
glm_poissonLecanoromycetesClass_simplified <- glm(LecanoromycetesClass ~ Inoculants + Aggregate + Plant + Time +
                                                    Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                    Plant:Time, family = poisson, data = tbl.GLM)
### Results
Anova(glm_poissonLecanoromycetesClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: LecanoromycetesClass
#                      LR Chisq Df Pr(>Chisq)   
# Inoculants             0.9786  1   0.322535   
# Aggregate              2.6706  2   0.263082   
# Plant                  6.7787  1   0.009225 **
# Time                   6.0443  1   0.013951 * 
# Inoculants:Aggregate   1.3039  2   0.521019   
# Inoculants:Plant       0.3464  1   0.556170   
# Aggregate:Plant        4.1769  2   0.123881   
# Inoculants:Time        0.4237  1   0.515077   
# Aggregate:Time         2.4763  2   0.289927   
# Plant:Time                     0              



#################################  DothideomycetesClass_model  
table(tbl.GLM$DothideomycetesClass)
### test mean 4.742843
mean(tbl.GLM$DothideomycetesClass)
### test variance 11.16776
var(tbl.GLM$DothideomycetesClass)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonDothideomycetesClass <- glm(DothideomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                       data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.8869922
deviance(glm_poissonDothideomycetesClass) / df.residual(glm_poissonDothideomycetesClass)
### Distribution of response variables
hist(tbl.GLM$DothideomycetesClass, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$DothideomycetesClass, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiDothideomycetesClass <- glm(DothideomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                     data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiDothideomycetesClass), residuals(glm_quasiDothideomycetesClass), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiDothideomycetesClass, type = "II")

glm_quasiDothideomycetesClass_simplified <- glm(DothideomycetesClass ~ Inoculants + Aggregate + Plant + Time +
                                                  Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                  Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiDothideomycetesClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: DothideomycetesClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             15.512  1  8.196e-05 ***
# Aggregate               5.119  2    0.07734 .  
# Plant                   2.208  1    0.13734    
# Time                   73.210  1  < 2.2e-16 ***
# Inoculants:Aggregate    0.840  2    0.65716    
# Inoculants:Plant        4.790  1    0.02862 *  
# Aggregate:Plant         2.494  2    0.28730    
# Inoculants:Time         3.879  1    0.04890 *  
# Aggregate:Time          6.095  2    0.04747 *  
# Plant:Time                     0               


#################################  LeotiomycetesClass_model  
table(tbl.GLM$LeotiomycetesClass)
### test mean 4.793549
mean(tbl.GLM$LeotiomycetesClass)
### test variance 32.33345
var(tbl.GLM$LeotiomycetesClass)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonLeotiomycetesClass <- glm(LeotiomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                     data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.538728
deviance(glm_poissonLeotiomycetesClass) / df.residual(glm_poissonLeotiomycetesClass)
### Distribution of response variables
hist(tbl.GLM$LeotiomycetesClass, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$LeotiomycetesClass, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiLeotiomycetesClass <- glm(LeotiomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                   data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiLeotiomycetesClass), residuals(glm_quasiLeotiomycetesClass), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiLeotiomycetesClass, type = "II")

glm_quasiLeotiomycetesClass_simplified <- glm(LeotiomycetesClass ~ Inoculants + Aggregate + Plant + Time +
                                                Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiLeotiomycetesClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: LeotiomycetesClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants            25.2958  1  4.918e-07 ***
# Aggregate             21.2890  2  2.383e-05 ***
# Plant                 12.6259  1  0.0003804 ***
# Time                  31.2834  1  2.230e-08 ***
# Inoculants:Aggregate   0.1925  2  0.9082520    
# Inoculants:Plant       0.0097  1  0.9214659    
# Aggregate:Plant        0.4208  2  0.8102429    
# Inoculants:Time       21.0865  1  4.390e-06 ***
# Aggregate:Time         0.5969  2  0.7419799    
# Plant:Time                     0               



#################################  MortierellomycetesClass_model  
table(tbl.GLM$MortierellomycetesClass)
### test mean 2.6225
mean(tbl.GLM$MortierellomycetesClass)
### test variance 7.923284
var(tbl.GLM$MortierellomycetesClass)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonMortierellomycetesClass <- glm(MortierellomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                          data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.296267
deviance(glm_poissonMortierellomycetesClass) / df.residual(glm_poissonMortierellomycetesClass)
### Distribution of response variables
hist(tbl.GLM$MortierellomycetesClass, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$MortierellomycetesClass, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiMortierellomycetesClass <- glm(MortierellomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                        data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiMortierellomycetesClass), residuals(glm_quasiMortierellomycetesClass), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiMortierellomycetesClass, type = "II")

glm_quasiMortierellomycetesClass_simplified <- glm(MortierellomycetesClass ~ Inoculants + Aggregate + Plant + Time +
                                                     Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                     Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiMortierellomycetesClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: MortierellomycetesClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             8.7342  1   0.003123 ** 
# Aggregate              0.1682  2   0.919332    
# Plant                 13.8217  1   0.000201 ***
# Time                   0.0782  1   0.779731    
# Inoculants:Aggregate   0.1499  2   0.927806    
# Inoculants:Plant       8.2317  1   0.004117 ** 
# Aggregate:Plant        0.4018  2   0.817998    
# Inoculants:Time        7.8181  1   0.005172 ** 
# Aggregate:Time         1.7019  2   0.427003    
# Plant:Time                     0             



#################################  GeminibasidiomycetesClass_model  
table(tbl.GLM$GeminibasidiomycetesClass)
### test mean 9.320037
mean(tbl.GLM$GeminibasidiomycetesClass)
### test variance 56.55115
var(tbl.GLM$GeminibasidiomycetesClass)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonGeminibasidiomycetesClass <- glm(GeminibasidiomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                            data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.629036
deviance(glm_poissonGeminibasidiomycetesClass) / df.residual(glm_poissonGeminibasidiomycetesClass)
### Distribution of response variables
hist(tbl.GLM$GeminibasidiomycetesClass, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$GeminibasidiomycetesClass, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiGeminibasidiomycetesClass <- glm(GeminibasidiomycetesClass ~ Inoculants * Aggregate * Plant * Time, 
                                          data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiGeminibasidiomycetesClass), residuals(glm_quasiGeminibasidiomycetesClass), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiGeminibasidiomycetesClass, type = "II")

glm_quasiGeminibasidiomycetesClass_simplified <- glm(GeminibasidiomycetesClass ~ Inoculants + Aggregate + Plant + Time +
                                                       Inoculants:Aggregate + Inoculants:Plant + Aggregate:Plant + Inoculants:Time + Aggregate:Time + 
                                                       Plant:Time, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiGeminibasidiomycetesClass_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: GeminibasidiomycetesClass
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              1.638  1   0.200620    
# Aggregate               4.892  2   0.086619 .  
# Plant                  46.410  1  9.594e-12 ***
# Time                  293.999  1  < 2.2e-16 ***
# Inoculants:Aggregate    3.577  2   0.167186    
# Inoculants:Plant        9.214  1   0.002401 ** 
# Aggregate:Plant         2.804  2   0.246110    
# Inoculants:Time         0.002  1   0.962470    
# Aggregate:Time          0.264  2   0.876449    
# Plant:Time                     0               


#################################Volcano
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(openxlsx)
df <- read.csv('Level1.csv',header = T)
head(df)
##Screening threshold determination: p＜0.001，|log2FC|＞2
pvalue = 0.05
log2FC = 2

df$group <- case_when(
  df$log2FoldChange > log2FC & df$pvalue < pvalue ~ "up",
  df$log2FoldChange < -log2FC & df$pvalue < pvalue ~ "down",
  TRUE ~ 'none'
)
head(df)
df$'-log10(pvalue)' <- -log10(df$pvalue) 
df$group <- factor(df$group, levels = c("up","down","none"))
write.xlsx(df, "Level1df.xlsx", rowNames = TRUE)

###### Plot
p1 <- ggplot(data = df,
             aes(x = log2FoldChange, y = -log10(pvalue), color = group, shape = group)) + 
  geom_point(size = 6) +
  scale_shape_manual(values = c(17, 15, 16)) + 
  scale_color_manual(values = c("up" = "#448DCD", "down" = "#F7AF34", "none" = "grey50")) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black', size = 9)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1
#
p2 <- p1 +
  geom_hline(yintercept = c(-log10(pvalue)),size = 0.7,color = "black",lty = "dashed") +
  geom_vline(xintercept = c(-log2FC, log2FC),size = 0.7,color = "black",lty = "dashed")
p2
pdf("Level1.pdf",width=6,height=5)
p2
dev.off()

################################# OPLS-DA
library(tidyverse)
library(ropls)
options(scipen = 9999)
#loading and reshaping data
files.csv <- fs::dir_ls("inputfile", recurse = TRUE, glob = "*.csv")
dat.group <-  map_dfc(files.csv,read_csv)
files.xlsx <- fs::dir_ls("inputfile",recurse = TRUE,glob = "*.xlsx")
dat <- map_dfc(files.xlsx,openxlsx::read.xlsx) %>% 
  select(c('variables',dat.group$sample)) %>% 
  column_to_rownames(var ='variables') %>% 
  t()
plsda<-opls(dat, dat.group$group,
            predI = 1, orthoI = 1,#orthoI = 1:opls-da; orthoI = 0:Pls-da
            log10L = F,
            crossvalI = nrow(dat),
            scaleC="pareto",#pareto scaling
            fig.pdfC = c("none", "interactive", "outputfile/plsda.pdf")[3],
            permI=200)
#extracting data of VIP:
vip <-  plsda@vipVn%>% as.data.frame() %>% 
  rename("VIP"=".") %>% rownames_to_column(var = 'variables')

dat_vip <- map_dfc(files.xlsx,openxlsx::read.xlsx) %>% 
  select(c('variables',dat.group$sample)) %>% 
  merge(vip,by="variables")

dat_vip %>% write_csv("outputfile/opls-daResult.csv")
# data and score data
## score data
plsdaScore <- data.frame(
  t1 =plsda@scoreMN,
  to1 =plsda@orthoScoreMN 
) %>% scale(center = T,scale = T) %>% 
  as.data.frame() %>% 
  rename(
    "t1"="p1",
    "to1"="o1"
  ) %>% 
  rownames_to_column(var = 'sample') %>% 
  merge(dat.group,by="sample")

t1Weight=sprintf("%.1f%%", plsda@modelDF[1,1] * 100);t1Weight
to1Weight=sprintf("%.1f%%", plsda@modelDF[2,1] * 100);to1Weight

R2X=plsda@modelDF[1,1]+plsda@modelDF[2,1]
R2Y=plsda@modelDF[1,3]+plsda@modelDF[2,3]
Q2Y=plsda@modelDF[1,6]+plsda@modelDF[2,6]

subTitle <- paste0("R2X=",R2X,"  R2Y=",R2Y,"  Q2Y=",Q2Y)
## opls-da 
oplsdaFig <- ggplot(plsdaScore,aes(x=t1,y=to1,color=group))+
  geom_point(aes(shape=group),size=5)+
  stat_ellipse(aes(fill=group),alpha=0.2,geom = "polygon")+
  theme_classic()+
  scale_fill_manual(values=c("#ae3b2b","#004164"))+
  scale_color_manual(values=c("#ae3b2b","#004164"))+
  labs(title = "OPLS-DA", 
       subtitle = subTitle, 
       x = paste0("t1(",t1Weight,")"),
       y = paste0("to1(",to1Weight,")"))+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)), 
    plot.subtitle = element_text(hjust = 0.5, size = rel(0.6)) 
  )
#saving
ggsave("outputfile/oplsda.png", width = 4, height = 4)
print(oplsdaFig)
ggsave("outputfile/oplsda.pdf", width = 4, height = 4,onefile=F)
print(oplsdaFig)
dev.off()     
# permutation data
permutation <- plsda@suppLs[["permMN"]] %>% as.data.frame() %>% 
  select(`R2Y(cum)`,`Q2(cum)`) %>% 
  mutate("number"=1:nrow(.)) %>% 
  pivot_longer(!number,names_to = "Type",values_to = "value") %>% 
  mutate(Type = str_replace(Type, "\\(cum\\)", ""))
permutation %>% write_csv("outputfile/200permutation.csv")

R2Xlabel <- paste0("R2X:",plsda@summaryDF[1,1])
R2Ylabel <- paste0("R2Y:",plsda@summaryDF[1,2],"\np=",
                   plsda@summaryDF[1,7],"(",plsda@summaryDF[1,7]*200,"/200)")
Q2lablel <- paste0("Q2:",plsda@summaryDF[1,3],"\np=",plsda@summaryDF[1,8],
                   "(",plsda@summaryDF[1,8]*200,"/200)")
calculate_bin_counts <- function(data,#Long-format data of permutation
                                 column, #Data column in the long format of permutation
                                 bins=30#Consistent with the 'bins' parameter of geom_histogram
) {
  if (!column %in% names(data)) {
    stop("zcp:Specified column does not exist in the dataframe!!!")
  }
  data_vector <- data[[column]]
  range_data <- range(data_vector, na.rm = TRUE)
  bin_width <- (range_data[2] - range_data[1]) / bins
  bin_breaks <- seq(from = range_data[1], to = range_data[2], by = bin_width)
  binned_data <- cut(data_vector, breaks = bin_breaks, include.lowest = TRUE, right = FALSE)
  bin_counts <- table(binned_data)
  return(bin_counts)
}
bin_counts <- calculate_bin_counts(data=permutation,column="value", bins=30)
max_counts <- max(bin_counts)
labels_data <- data.frame(
  x = c(plsda@summaryDF[1,1], plsda@summaryDF[1,2], plsda@summaryDF[1,3]),
  y = c(max_counts*0.4, max_counts*0.8,max_counts*1.2), 
  label = c(R2Xlabel, R2Ylabel,Q2lablel)
)
#Plotting Permutation Histogram
PermutationFig <- ggplot()+
  geom_histogram(data = permutation,aes(x=value,fill=Type),
                 bins=30,# Keep consistent with the bins parameter in calculate_bin_counts function!
                 alpha=0.6,
                 size=0.3,
                 color="black")+
  scale_fill_manual(values=c("#ae3b2b","#004164"))+
  xlim(c(NA, max(labels_data$x) * 1.2))+
  ylim(c(NA, max(labels_data$y) * 1.2))+
  geom_label(data = labels_data, 
             aes(x = x, y = y, label = label), 
             size = 2.5, color = "black",
             fill = "white", label.padding = unit(0.1, "lines"),
             label.size = 0.25)+
  geom_segment(data = labels_data, aes(x = x, xend = x, y = y, yend = 0,color=label),
               arrow = arrow(type = "closed", length = unit(0.15, "cm")),
               alpha=0.5,
               show.legend = FALSE)+
  labs(title = "Permutation of OPLS-DA", 
       x = "Permutations",
       y ="Frequency")+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)))

#saving
ggsave("outputfile/Permutation.png", width = 6, height = 4)
print(PermutationFig)
ggsave("outputfile/Permutation.pdf", width = 6, height = 4,onefile=F)
print(PermutationFig)
dev.off()  


################################# Random Forest
library(randomForest)
library(rfPermute)
library(ggplot2)
library(A3)
library(openxlsx)
##################
otu <- read.csv('Rootexudate_Level3.csv',row.names = 1)
otu <- data.frame(t(otu))
otu$year <- c(rep("CK",11),rep("TR",11))
otu$year <-as.factor(otu$year)
#Random forest calculation (1000 decision trees by default)
set.seed(123)
otu_forest <- rfPermute(year~., data = otu, importance = TRUE, ntree = 1000,nrep = 1000, num.cores = 1)
otu_forest
otu_forest$rf
##Importance assessment of exudates
importance_year_forest <- data.frame(importance(otu_forest, scale = TRUE), check.names = FALSE)#scale = TRUE标准化
importance_year_forest
write.csv(importance_year_forest,"importance_year_forest.csv")
# First divide the training set and test set
set.seed(111) 
train_idx <- sample(seq_len(nrow(otu)), size = 1 * nrow(otu)) # Randomly select 100% of the data as the training set
train_data <- otu[train_idx, ]
test_data <- otu[-train_idx, ]

# Ten-fold cross-validation with 10 repetitions on the training set
otu_cv <- replicate(10, rfcv(train_data[-ncol(train_data)], train_data$year, cv.fold = 10, step = 1.5), simplify = FALSE)

# Extraction cross-validation error rate
otu_cv <- data.frame(sapply(otu_cv, '[[', 'error.cv'))
otu_cv$otus <- rownames(otu_cv)
otu_cv <- reshape2::melt(otu_cv, id = 'otus')
otu_cv$otus <- as.numeric(as.character(otu_cv$otus))

# Calculate the average error rate for each number of features
otu_cv.mean <- aggregate(otu_cv$value, by = list(otu_cv$otus), FUN = mean)
head(otu_cv.mean, 10)

# Saved
write.csv(otu_cv.mean, "otu_cv.mean.csv")
#
p <- ggplot(otu_cv.mean, aes(Group.1, x)) +
  geom_line() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +  
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')
p
###################Second step
importance_otu <- importance_year_forest[order(importance_year_forest[,5], decreasing = TRUE), ]# decreasing = TRUE降序的意思
otu.select <- importance_otu[1:11, ]
otu.select <- as.data.frame(otu.select)
otu.select <- otu[,row.names(otu.select)]
otu.select$year <- otu$year
#Running
set.seed(111)
otu_select_forest <- rfPermute(year~., data = otu.select, importance = TRUE, ntree = 1000,nrep = 1000, num.cores = 1)
importance_otu_select_forest <- data.frame(importance(otu_select_forest, scale = TRUE), check.names = FALSE)#scale = TRUE标准化
importance_otu_select_forest$OTU_name <- rownames(importance_otu_select_forest)
importance_otu_select_forest <- importance_otu_select_forest[order(importance_otu_select_forest[,5], decreasing = TRUE), ]
importance_otu_select_forest$OTU_name <- factor(importance_otu_select_forest$OTU_name, levels =importance_otu_select_forest$OTU_name )
otu_select_forest$rf

for (OTU in rownames(importance_otu_select_forest)) {
  if (importance_otu_select_forest[OTU,'MeanDecreaseAccuracy.pval'] >= 0.05) importance_otu_select_forest[OTU,'sig'] <- ''
  else if (importance_otu_select_forest[OTU,'MeanDecreaseAccuracy.pval'] >= 0.01 & importance_otu_select_forest[OTU,'MeanDecreaseAccuracy.pval'] < 0.05) importance_otu_select_forest[OTU,'sig'] <- '*'
  else if (importance_otu_select_forest[OTU,'MeanDecreaseAccuracy.pval'] >= 0.001 & importance_otu_select_forest[OTU,'MeanDecreaseAccuracy.pval'] < 0.01) importance_otu_select_forest[OTU,'sig'] <- '**'
  else if (importance_otu_select_forest[OTU,'MeanDecreaseAccuracy.pval'] < 0.001) importance_otu_select_forest[OTU,'sig'] <- '***'
}
p1 <- ggplot(importance_otu_select_forest, aes(OTU_name,MeanDecreaseAccuracy )) +
  geom_col(width = 0.5, fill = '#FFC068', color = NA) +#width = 0.5为柱状图宽度, fill = '#FFC068'为柱状图填充颜色
  labs(title =NULL, x = NULL, y = 'Mean Decrease Accuracy', fill = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +#angle = 45文字倾斜45度，hjust = 1 为文与坐标轴的距离，越大则距离约远
  scale_y_continuous(expand = c(0, 0), limit = c(0, 30))+#limit = c(0, 3000) 为y轴坐标，表示从0开始到3000，根据数据自己调整  不然会报错 
  annotate('text', label = sprintf('italic(R^2) == %.2f', 7.69), x = 15, y = 28, size = 3, parse = TRUE)+
  geom_text(aes(x = OTU_name,  y=MeanDecreaseAccuracy, label = sig),position="stack",stat="identity")

p1
ggsave("rf1.PDF",p1,width = 7,height = 5)
# Output
write.xlsx(importance_otu_select_forest[, c("OTU_name", "MeanDecreaseAccuracy", "sig")], 
           file = "MeanDecreaseAccuracy_results.xlsx", 
           row.names = FALSE)
# Use A3 package to test significance
set.seed(123)
otu_forest.pval <- a3(year ~ ., data = otu, model.fn = randomForest, p.acc = 0.001, model.args = list(importance = TRUE, ntree = 1000))
print(otu_forest.pval)


#################### Relationship between soil properties and microbial community
library(dplyr)
library(linkET)
library(ggplot2)
library(readr)
library(openxlsx)
library(reshape2)
# Input data
env <- read.csv(file.choose(), header = T, row.names = 1) # 替换为实际环境因子数据路径
speciese <- read.csv(file.choose(), header = T, row.names = 1) # 替换为实际物种数据路径
common_samples <- intersect(rownames(env), rownames(speciese))
env <- env[common_samples, , drop = FALSE]
speciese <- speciese[common_samples, , drop = FALSE]
# Mantel
mantel_results <- mantel_test(
  speciese, env,
  spec_select = list(
    bacteria = 1:12704,  # bacteria
    fungi = 12705:15885, # fungi
    microbial = 1:15885 # microbe
  )
) %>%
  mutate(
    rd = cut(
      r, breaks = c(-Inf, 0.2, 0.4, Inf),
      labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")
    ),
    pd = cut(
      p, breaks = c(-Inf, 0.01, 0.05, Inf),
      labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")
    )
  )
# Mantel results
print(mantel_results)
write.xlsx(mantel_results, file = "Mantel_and_Spearman_Results.xlsx", overwrite = TRUE)
env <- env %>% select_if(is.numeric)
env <- na.omit(env)
# Spearman
env_corr <- correlate(env, method = "spearman")
qcorrplot(env_corr, type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(
    aes(
      colour = pd,
      size = rd,
      linetype = factor(sign(r))
    ), 
    data = mantel_results, curvature = 0.1
  ) +
  scale_fill_gradientn(
    colours = c("#F7AF34", "#FFFFFF", "#448DCD"),
    limits = c(-1, 1),      
    name = "Spearman's r" 
  ) +
  scale_size_manual(
    values = c(0.5, 1, 2),
    name = "Mantel's r"   
  ) +
  scale_colour_manual(
    values = c("< 0.01" = "purple", "0.01 - 0.05" = "green", ">= 0.05" = "grey"), 
    name = "P value"    
  ) +
  scale_linetype_manual(
    values = c("dashed", "solid"),
    labels = c("Negative", "Positive"),
    name = "MantelR Sign"
  ) +
  guides(
    size = guide_legend(
      title = "Mantel's r",
      override.aes = list(colour = "grey35"),
      order = 2
    ),
    colour = guide_legend(
      title = "Mantel's p",
      override.aes = list(size = 3),
      order = 1
    ),
    linetype = guide_legend(
      title = "MantelR Sign",
      override.aes = list(size = 1), 
      order = 3
    ),
    fill = guide_colorbar(
      title = "Spearman's r", 
      order = 4
    )
  ) +
  labs(
    title = "Mantel Test with Environmental Factors",
    x = NULL, 
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right" 
  )

#  Spearman and significance
cor_results <- psych::corr.test(env, method = "spearman") 
cor_matrix <- cor_results$r 
p_matrix <- cor_results$p  
# BH adjusted
p_adjusted <- apply(p_matrix, 2, p.adjust, method = "BH")
cor_long <- melt(cor_matrix, varnames = c("Var1", "Var2"), value.name = "Spearman_r")
p_long <- melt(p_matrix, varnames = c("Var1", "Var2"), value.name = "Spearman_p")
p_adj_long <- melt(p_adjusted, varnames = c("Var1", "Var2"), value.name = "BH_corrected_p")
cor_data <- cor_long %>%
  left_join(p_long, by = c("Var1", "Var2")) %>%
  left_join(p_adj_long, by = c("Var1", "Var2"))
cor_data <- cor_data %>%
  mutate(
    Spearman_signif = case_when(
      Spearman_p < 0.001 ~ "***",        # p < 0.001 用 "***"
      Spearman_p < 0.01 ~ "**",          # 0.001 <= p < 0.01 用 "**"
      Spearman_p < 0.05 ~ "*",           # 0.01 <= p < 0.05 用 "*"
      TRUE ~ ""                          # p >= 0.05 无标记
    )
  )
# Output
write.xlsx(cor_data, file = "Spearman_Correlation_Results.xlsx", overwrite = TRUE)

############################################### Relationships between MWD and microbial community
library(vegan)
# Input data
bacteria_a <- read.csv("BacteriaCommunity_A.csv", row.names = 1)
bacteria_b <- read.csv("BacteriaCommunity_B.csv", row.names = 1)
bacteria_c <- read.csv("BacteriaCommunity_C.csv", row.names = 1)
fungi_a <- read.csv("FungiCommunity_A.csv", row.names = 1)
fungi_b <- read.csv("FungiCommunity_B.csv", row.names = 1)
fungi_c <- read.csv("FungiCommunity_C.csv", row.names = 1)
mwd_data <- read.csv("MWD.csv", row.names = 1)
#  Shannon and Richness
alpha_diversity <- function(data) {
  shannon <- diversity(t(data), index = "shannon")         # Shannon
  richness <- specnumber(t(data))                          # Richness
  list(shannon = shannon, richness = richness)
}
# Alpha
bacteria_a_alpha <- alpha_diversity(bacteria_a)
bacteria_b_alpha <- alpha_diversity(bacteria_b)
bacteria_c_alpha <- alpha_diversity(bacteria_c)
fungi_a_alpha <- alpha_diversity(fungi_a)
fungi_b_alpha <- alpha_diversity(fungi_b)
fungi_c_alpha <- alpha_diversity(fungi_c)
# spearman
correlation_test <- function(diversity_metric, mwd) {
  cor.test(diversity_metric, mwd, method = "spearman")
}
# Shannon and MWD
cor_shannon_bacteria_a <- correlation_test(bacteria_a_alpha$shannon, mwd_data$MWD)
cor_shannon_bacteria_b <- correlation_test(bacteria_b_alpha$shannon, mwd_data$MWD)
cor_shannon_bacteria_c <- correlation_test(bacteria_c_alpha$shannon, mwd_data$MWD)
cor_shannon_fungi_a <- correlation_test(fungi_a_alpha$shannon, mwd_data$MWD)
cor_shannon_fungi_b <- correlation_test(fungi_b_alpha$shannon, mwd_data$MWD)
cor_shannon_fungi_c <- correlation_test(fungi_c_alpha$shannon, mwd_data$MWD)
# Richness and MWD
cor_richness_bacteria_a <- correlation_test(bacteria_a_alpha$richness, mwd_data$MWD)
cor_richness_bacteria_b <- correlation_test(bacteria_b_alpha$richness, mwd_data$MWD)
cor_richness_bacteria_c <- correlation_test(bacteria_c_alpha$richness, mwd_data$MWD)
cor_richness_fungi_a <- correlation_test(fungi_a_alpha$richness, mwd_data$MWD)
cor_richness_fungi_b <- correlation_test(fungi_b_alpha$richness, mwd_data$MWD)
cor_richness_fungi_c <- correlation_test(fungi_c_alpha$richness, mwd_data$MWD)
# Merge
correlation_results <- data.frame(
  Community = c("bacteria_a", "bacteria_b", "bacteria_c", "fungi_a", "fungi_b", "fungi_c", 
                "bacteria_a", "bacteria_b", "bacteria_c", "fungi_a", "fungi_b", "fungi_c"),
  Metric = c(rep("Shannon", 6), rep("Richness", 6)),
  Spearman_Correlation = c(cor_shannon_bacteria_a$estimate, cor_shannon_bacteria_b$estimate, cor_shannon_bacteria_c$estimate, 
                           cor_shannon_fungi_a$estimate, cor_shannon_fungi_b$estimate, cor_shannon_fungi_c$estimate,
                           cor_richness_bacteria_a$estimate, cor_richness_bacteria_b$estimate, cor_richness_bacteria_c$estimate,
                           cor_richness_fungi_a$estimate, cor_richness_fungi_b$estimate, cor_richness_fungi_c$estimate),
  P_Value = c(cor_shannon_bacteria_a$p.value, cor_shannon_bacteria_b$p.value, cor_shannon_bacteria_c$p.value, 
              cor_shannon_fungi_a$p.value, cor_shannon_fungi_b$p.value, cor_shannon_fungi_c$p.value,
              cor_richness_bacteria_a$p.value, cor_richness_bacteria_b$p.value, cor_richness_bacteria_c$p.value,
              cor_richness_fungi_a$p.value, cor_richness_fungi_b$p.value, cor_richness_fungi_c$p.value)
)
print(correlation_results)
# Saved
write.csv(correlation_results, "Alpha_Diversity_Correlation_with_MWD.csv", row.names = FALSE)
cat("Results saved to 'Alpha_Diversity_Correlation_with_MWD.csv'\n")
###############
# Bray-Curtis
calculate_bray_curtis <- function(data) {
  vegdist(t(data), method = "bray")
}
bacteria_a_dist <- calculate_bray_curtis(bacteria_a)
bacteria_b_dist <- calculate_bray_curtis(bacteria_b)
bacteria_c_dist <- calculate_bray_curtis(bacteria_c)
fungi_a_dist <- calculate_bray_curtis(fungi_a)
fungi_b_dist <- calculate_bray_curtis(fungi_b)
fungi_c_dist <- calculate_bray_curtis(fungi_c)
#  MWD
mwd_dist <- dist(mwd_data$MWD, method = "euclidean")
# Mantel
mantel_bacteria_a <- mantel_test(bacteria_a_dist, mwd_dist)
mantel_bacteria_b <- mantel_test(bacteria_b_dist, mwd_dist)
mantel_bacteria_c <- mantel_test(bacteria_c_dist, mwd_dist)
mantel_fungi_a <- mantel_test(fungi_a_dist, mwd_dist)
mantel_fungi_b <- mantel_test(fungi_b_dist, mwd_dist)
mantel_fungi_c <- mantel_test(fungi_c_dist, mwd_dist)
# Merge
mantel_results <- data.frame(
  Community = c("Bacteria_A", "Bacteria_B", "Bacteria_C", "Fungi_A", "Fungi_B", "Fungi_C"),
  Mantel_R = c(
    mantel_bacteria_a$statistic,
    mantel_bacteria_b$statistic,
    mantel_bacteria_c$statistic,
    mantel_fungi_a$statistic,
    mantel_fungi_b$statistic,
    mantel_fungi_c$statistic
  ),
  P_Value = c(
    mantel_bacteria_a$signif,
    mantel_bacteria_b$signif,
    mantel_bacteria_c$signif,
    mantel_fungi_a$signif,
    mantel_fungi_b$signif,
    mantel_fungi_c$signif
  )
)
print(mantel_results)
# Saved
write.csv(mantel_results, "Mantel_Test_Results.csv", row.names = FALSE)
cat("Mantel test results saved to 'Mantel_Test_Results.csv'\n")

library(stats)
mwd_data <- read.csv("MWD.csv", row.names = 1)  #据
bacteria_a <- read.csv("BacteriaCommunity_A.csv", row.names = 1)
bacteria_b <- read.csv("BacteriaCommunity_B.csv", row.names = 1)
bacteria_c <- read.csv("BacteriaCommunity_C.csv", row.names = 1)  
fungi_a <- read.csv("FungiCommunity_A.csv", row.names = 1)
fungi_b <- read.csv("FungiCommunity_B.csv", row.names = 1) 
fungi_c <- read.csv("FungiCommunity_C.csv", row.names = 1) 
rownames(mwd_data) <- tolower(gsub(" ", "", rownames(mwd_data)))
colnames(bacteria_a) <- tolower(gsub(" ", "", colnames(bacteria_a)))
colnames(bacteria_b) <- tolower(gsub(" ", "", colnames(bacteria_b)))
colnames(bacteria_c) <- tolower(gsub(" ", "", colnames(bacteria_c)))
colnames(fungi_a) <- tolower(gsub(" ", "", colnames(fungi_a)))
colnames(fungi_b) <- tolower(gsub(" ", "", colnames(fungi_b)))
colnames(fungi_c) <- tolower(gsub(" ", "", colnames(fungi_c)))
bacteria_a <- bacteria_a[, rownames(mwd_data)]
bacteria_b <- bacteria_b[, rownames(mwd_data)]
bacteria_c <- bacteria_c[, rownames(mwd_data)]
fungi_a <- fungi_a[, rownames(mwd_data)]
fungi_b <- fungi_b[, rownames(mwd_data)]
fungi_c <- fungi_c[, rownames(mwd_data)]
bacteria_a <- t(bacteria_a)
bacteria_b <- t(bacteria_b)
bacteria_c <- t(bacteria_c)
fungi_a <- t(fungi_a)
fungi_b <- t(fungi_b)
fungi_c <- t(fungi_c)
y <- mwd_data$MWD
bacteria_a <- bacteria_a[, apply(bacteria_a, 2, var) != 0]
bacteria_b <- bacteria_b[, apply(bacteria_b, 2, var) != 0]
bacteria_c <- bacteria_c[, apply(bacteria_c, 2, var) != 0]
fungi_a <- fungi_a[, apply(fungi_a, 2, var) != 0]
fungi_b <- fungi_b[, apply(fungi_b, 2, var) != 0]
fungi_c <- fungi_c[, apply(fungi_c, 2, var) != 0]
# PCA
bacteria_a_pca <- prcomp(bacteria_a, scale. = TRUE)
bacteria_b_pca <- prcomp(bacteria_b, scale. = TRUE)
bacteria_c_pca <- prcomp(bacteria_c, scale. = TRUE)
fungi_a_pca <- prcomp(fungi_a, scale. = TRUE)
fungi_b_pca <- prcomp(fungi_b, scale. = TRUE)
fungi_c_pca <- prcomp(fungi_c, scale. = TRUE)
# PC1 and PC2
bacteria_a_pc1 <- bacteria_a_pca$x[, 1]
bacteria_b_pc1 <- bacteria_b_pca$x[, 1]
bacteria_c_pc1 <- bacteria_c_pca$x[, 1]
fungi_a_pc1 <- fungi_a_pca$x[, 1]
fungi_b_pc1 <- fungi_b_pca$x[, 1]
fungi_c_pc1 <- fungi_c_pca$x[, 1]
bacteria_a_pc2 <- bacteria_a_pca$x[, 2]
bacteria_b_pc2 <- bacteria_b_pca$x[, 2]
bacteria_c_pc2 <- bacteria_c_pca$x[, 2]
fungi_a_pc2 <- fungi_a_pca$x[, 2]
fungi_b_pc2 <- fungi_b_pca$x[, 2]
fungi_c_pc2 <- fungi_c_pca$x[, 2]
summary(bacteria_a_pca)
summary(bacteria_b_pca)
summary(bacteria_c_pca)
summary(fungi_a_pca)
summary(fungi_b_pca)
summary(fungi_c_pca)
# Spearman
cor_bacteria_a_pc1 <- cor.test(bacteria_a_pc1, y, method = "spearman")
cor_bacteria_b_pc1 <- cor.test(bacteria_b_pc1, y, method = "spearman")
cor_bacteria_c_pc1 <- cor.test(bacteria_c_pc1, y, method = "spearman")
cor_fungi_a_pc1 <- cor.test(fungi_a_pc1, y, method = "spearman")
cor_fungi_b_pc1 <- cor.test(fungi_b_pc1, y, method = "spearman")
cor_fungi_c_pc1 <- cor.test(fungi_c_pc1, y, method = "spearman")
cor_bacteria_a_pc2 <- cor.test(bacteria_a_pc2, y, method = "spearman")
cor_bacteria_b_pc2 <- cor.test(bacteria_b_pc2, y, method = "spearman")
cor_bacteria_c_pc2 <- cor.test(bacteria_c_pc2, y, method = "spearman")
cor_fungi_a_pc2 <- cor.test(fungi_a_pc2, y, method = "spearman")
cor_fungi_b_pc2 <- cor.test(fungi_b_pc2, y, method = "spearman")
cor_fungi_c_pc2 <- cor.test(fungi_c_pc2, y, method = "spearman")
# Merge
correlation_results <- data.frame(
  Community = c("bacteria_a", "bacteria_b", "bacteria_c", "fungi_a", "fungi_b", "fungi_c",
                "bacteria_a", "bacteria_b", "bacteria_c", "fungi_a", "fungi_b", "fungi_c"),
  Metric = c(rep("PC1", 6), rep("PC2", 6)),
  Spearman_Correlation = c(cor_bacteria_a_pc1$estimate, cor_bacteria_b_pc1$estimate, cor_bacteria_c_pc1$estimate,
                           cor_fungi_a_pc1$estimate, cor_fungi_b_pc1$estimate, cor_fungi_c_pc1$estimate,
                           cor_bacteria_a_pc2$estimate, cor_bacteria_b_pc2$estimate, cor_bacteria_c_pc2$estimate,
                           cor_fungi_a_pc2$estimate, cor_fungi_b_pc2$estimate, cor_fungi_c_pc2$estimate),
  P_Value = c(cor_bacteria_a_pc1$p.value, cor_bacteria_b_pc1$p.value, cor_bacteria_c_pc1$p.value,
              cor_fungi_a_pc1$p.value, cor_fungi_b_pc1$p.value, cor_fungi_c_pc1$p.value,
              cor_bacteria_a_pc2$p.value, cor_bacteria_b_pc2$p.value, cor_bacteria_c_pc2$p.value,
              cor_fungi_a_pc2$p.value, cor_fungi_b_pc2$p.value, cor_fungi_c_pc2$p.value)
)
print(correlation_results)
# Saved
write.csv(correlation_results, "PCA_Correlation_with_MWD.csv", row.names = FALSE)
cat("PCA correlation results saved to 'PCA_Correlation_with_MWD.csv'\n")

library(ggplot2)
library(reshape2)
correlation_results <- data.frame(
  Community = c("bacteria_a", "bacteria_b", "bacteria_c", "fungi_a", "fungi_b", "fungi_c",
                "bacteria_a", "bacteria_b", "bacteria_c", "fungi_a", "fungi_b", "fungi_c"),
  Metric = c(rep("PC1", 6), rep("PC2", 6)),
  Spearman_Correlation = c(cor_bacteria_a_pc1$estimate, cor_bacteria_b_pc1$estimate, cor_bacteria_c_pc1$estimate,
                           cor_fungi_a_pc1$estimate, cor_fungi_b_pc1$estimate, cor_fungi_c_pc1$estimate,
                           cor_bacteria_a_pc2$estimate, cor_bacteria_b_pc2$estimate, cor_bacteria_c_pc2$estimate,
                           cor_fungi_a_pc2$estimate, cor_fungi_b_pc2$estimate, cor_fungi_c_pc2$estimate),
  P_Value = c(cor_bacteria_a_pc1$p.value, cor_bacteria_b_pc1$p.value, cor_bacteria_c_pc1$p.value,
              cor_fungi_a_pc1$p.value, cor_fungi_b_pc1$p.value, cor_fungi_c_pc1$p.value,
              cor_bacteria_a_pc2$p.value, cor_bacteria_b_pc2$p.value, cor_bacteria_c_pc2$p.value,
              cor_fungi_a_pc2$p.value, cor_fungi_b_pc2$p.value, cor_fungi_c_pc2$p.value)
)
# Mantel
mantel_results <- data.frame(
  Community = c("bacteria_a", "bacteria_b", "bacteria_c", "fungi_a", "fungi_b", "fungi_c"),
  Metric = "Mantel",
  Spearman_Correlation = c(mantel_bacteria_a$statistic, mantel_bacteria_b$statistic, mantel_bacteria_c$statistic,
                           mantel_fungi_a$statistic, mantel_fungi_b$statistic, mantel_fungi_c$statistic),
  P_Value = c(mantel_bacteria_a$signif, mantel_bacteria_b$signif, mantel_bacteria_c$signif,
              mantel_fungi_a$signif, mantel_fungi_b$signif, mantel_fungi_c$signif)
)
# Correlation
diversity_results <- data.frame(
  Community = c("bacteria_a", "bacteria_b", "bacteria_c", "fungi_a", "fungi_b", "fungi_c"),
  Metric = c("Shannon", "Shannon", "Shannon", "Shannon", "Shannon", "Shannon"),  # 这里可以替换为 Shannon 或 Richness
  Spearman_Correlation = c(cor_shannon_bacteria_a$estimate, cor_shannon_bacteria_b$estimate, cor_shannon_bacteria_c$estimate,
                           cor_shannon_fungi_a$estimate, cor_shannon_fungi_b$estimate, cor_shannon_fungi_c$estimate),
  P_Value = c(cor_shannon_bacteria_a$p.value, cor_shannon_bacteria_b$p.value, cor_shannon_bacteria_c$p.value,
              cor_shannon_fungi_a$p.value, cor_shannon_fungi_b$p.value, cor_shannon_fungi_c$p.value)
)

# Merge
all_results <- rbind(correlation_results, mantel_results, diversity_results)
ggplot(all_results, aes(x = Community, y = Spearman_Correlation, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Correlation of Different Metrics with MWD",
       x = "Community", y = "Spearman Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("PC1" = "blue", "PC2" = "green", "Mantel" = "purple", "Shannon" = "orange")) +
  theme(legend.title = element_blank()) +
  geom_text(aes(label = paste0("P = ", format(P_Value, digits = 3))),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5)

############################################### Relationship between exudate and MWD
library(ggplot2)
library(reshape2)
mwd_data <- read.csv("MWD.csv", row.names = 1)
level1_data <- read.csv("Rootexudate_Level1.csv", row.names = 1)
level2_data <- read.csv("Rootexudate_Level2.csv", row.names = 1)
level3_data <- read.csv("Rootexudate_Level3.csv", row.names = 1)
combined_results <- data.frame(
  Metric = c("Mantel", "Mantel", "Mantel", 
             "Total Sum", "Total Sum", "Total Sum", 
             "PCA_PC1", "PCA_PC1", "PCA_PC1", 
             "PCA_PC2", "PCA_PC2", "PCA_PC2"),
  Level = rep(c("Level1", "Level2", "Level3"), 4),
  Correlation = c(
    mantel_results$Mantel_R,
    c(cor_sum_level1$estimate, cor_sum_level2$estimate, cor_sum_level3$estimate), 
    c(cor_level1_pc1$estimate, cor_level2_pc1$estimate, cor_level3_pc1$estimate),
    c(cor_level1_pc2$estimate, cor_level2_pc2$estimate, cor_level3_pc2$estimate)
  )
)
heatmap_data <- melt(combined_results, id.vars = c("Metric", "Level"), variable.name = "MetricType", value.name = "Correlation")
ggplot(heatmap_data, aes(x = Metric, y = Level, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#F7AF34", high = "#448DCD", mid = "white", midpoint = 0,
                       name = "Correlation") +
  labs(title = "Combined Heatmap of Mantel, Sum, and PCA Correlations",
       x = "Metric", y = "Levels") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###################################### Relationship between exudate and microbial community
library(vegan)
library(openxlsx)
library(ggplot2)
level1_data <- read.xlsx("Rootexudate_Level1.xlsx", sheet = 1, rowNames = TRUE)
level2_data <- read.xlsx("Rootexudate_Level2.xlsx", sheet = 1, rowNames = TRUE)
level3_data <- read.xlsx("Rootexudate_Level3.xlsx", sheet = 1, rowNames = TRUE)
level1_data <- t(level1_data)
level2_data <- t(level2_data)
level3_data <- t(level3_data)
level1_data[is.na(level1_data)] <- 0
level2_data[is.na(level2_data)] <- 0
level3_data[is.na(level3_data)] <- 0

bacterial_a_data <- read.xlsx("Bacterial_A.xlsx", sheet = 1, rowNames = TRUE)
bacterial_b_data <- read.xlsx("Bacterial_B.xlsx", sheet = 1, rowNames = TRUE)
bacterial_c_data <- read.xlsx("Bacterial_C.xlsx", sheet = 1, rowNames = TRUE)
fungal_a_data <- read.xlsx("Fungal_A.xlsx", sheet = 1, rowNames = TRUE)
fungal_b_data <- read.xlsx("Fungal_B.xlsx", sheet = 1, rowNames = TRUE)
fungal_c_data <- read.xlsx("Fungal_C.xlsx", sheet = 1, rowNames = TRUE)
bacterial_a_data <- t(bacterial_a_data)
bacterial_b_data <- t(bacterial_b_data)
bacterial_c_data <- t(bacterial_c_data)
fungal_a_data <- t(fungal_a_data)
fungal_b_data <- t(fungal_b_data)
fungal_c_data <- t(fungal_c_data)
bacterial_a_data[is.na(bacterial_a_data)] <- 0
bacterial_b_data[is.na(bacterial_b_data)] <- 0
bacterial_c_data[is.na(bacterial_c_data)] <- 0
fungal_a_data[is.na(fungal_a_data)] <- 0
fungal_b_data[is.na(fungal_b_data)] <- 0
fungal_c_data[is.na(fungal_c_data)] <- 0

# bray
level1_dist <- vegdist(level1_data, method = "bray")
level2_dist <- vegdist(level2_data, method = "bray")
level3_dist <- vegdist(level3_data, method = "bray")

# bray
bacterial_a_dist <- vegdist(bacterial_a_data, method = "bray")
bacterial_b_dist <- vegdist(bacterial_b_data, method = "bray")
bacterial_c_dist <- vegdist(bacterial_c_data, method = "bray")
fungal_a_dist <- vegdist(fungal_a_data, method = "bray")
fungal_b_dist <- vegdist(fungal_b_data, method = "bray")
fungal_c_dist <- vegdist(fungal_c_data, method = "bray")

# Mantel
compute_mantel <- function(metabolite_dist, microbial_dist, label) {
  mantel_result <- mantel(metabolite_dist, microbial_dist, method = "spearman", permutations = 999)
  return(data.frame(Group = label,
                    Mantel_R = mantel_result$statistic,
                    P_Value = mantel_result$signif))
}
# Mantel
results <- rbind(
  compute_mantel(level1_dist, bacterial_a_dist, "Level1-Bacterial_A"),
  compute_mantel(level1_dist, bacterial_b_dist, "Level1-Bacterial_B"),
  compute_mantel(level1_dist, bacterial_c_dist, "Level1-Bacterial_C"),
  compute_mantel(level1_dist, fungal_a_dist, "Level1-Fungal_A"),
  compute_mantel(level1_dist, fungal_b_dist, "Level1-Fungal_B"),
  compute_mantel(level1_dist, fungal_c_dist, "Level1-Fungal_C"),
  compute_mantel(level2_dist, bacterial_a_dist, "Level2-Bacterial_A"),
  compute_mantel(level2_dist, bacterial_b_dist, "Level2-Bacterial_B"),
  compute_mantel(level2_dist, bacterial_c_dist, "Level2-Bacterial_C"),
  compute_mantel(level2_dist, fungal_a_dist, "Level2-Fungal_A"),
  compute_mantel(level2_dist, fungal_b_dist, "Level2-Fungal_B"),
  compute_mantel(level2_dist, fungal_c_dist, "Level2-Fungal_C"),
  compute_mantel(level3_dist, bacterial_a_dist, "Level3-Bacterial_A"),
  compute_mantel(level3_dist, bacterial_b_dist, "Level3-Bacterial_B"),
  compute_mantel(level3_dist, bacterial_c_dist, "Level3-Bacterial_C"),
  compute_mantel(level3_dist, fungal_a_dist, "Level3-Fungal_A"),
  compute_mantel(level3_dist, fungal_b_dist, "Level3-Fungal_B"),
  compute_mantel(level3_dist, fungal_c_dist, "Level3-Fungal_C")
)
# Output
write.xlsx(results, file = "Mantel_Test_Results_All_Combinations.xlsx")
ggplot(results, aes(x = Group, y = Mantel_R, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(P_Value, 3)), vjust = -0.5) +
  labs(title = "Mantel Test Results for Different Levels and Groups", 
       y = "Mantel Statistic R", x = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################ Relationship between exudate and SOC
library(vegan)
library(openxlsx)
library(ggplot2)
# SOC
SOC_A <- as.matrix(read.xlsx("SOC_A.xlsx", rowNames = TRUE))
SOC_B <- as.matrix(read.xlsx("SOC_B.xlsx", rowNames = TRUE))
SOC_C <- as.matrix(read.xlsx("SOC_C.xlsx", rowNames = TRUE))
# Rootexudate
Rootexudate_Level1 <- read.csv("Rootexudate_Level1.csv", row.names = 1)
Rootexudate_Level2 <- read.csv("Rootexudate_Level2.csv", row.names = 1)
Rootexudate_Level3 <- read.csv("Rootexudate_Level3.csv", row.names = 1)
Rootexudate_Level1 <- t(as.matrix(Rootexudate_Level1))
Rootexudate_Level2 <- t(as.matrix(Rootexudate_Level2))
Rootexudate_Level3 <- t(as.matrix(Rootexudate_Level3))

# Mantel
# Bray-Curtis
  SOC_dist <- dist(SOC) 
  Rootexudate_dist <- vegdist(Rootexudate, method = "bray")  # Bray-Curtis
  # Mantel
  mantel_result <- mantel(SOC_dist, Rootexudate_dist, method = "pearson")
  return(mantel_result)
# Mantel
mantel_results <- list(
  SOC_A_Level1 = perform_mantel(SOC_A, Rootexudate_Level1),
  SOC_A_Level2 = perform_mantel(SOC_A, Rootexudate_Level2),
  SOC_A_Level3 = perform_mantel(SOC_A, Rootexudate_Level3),
  
  SOC_B_Level1 = perform_mantel(SOC_B, Rootexudate_Level1),
  SOC_B_Level2 = perform_mantel(SOC_B, Rootexudate_Level2),
  SOC_B_Level3 = perform_mantel(SOC_B, Rootexudate_Level3),
  
  SOC_C_Level1 = perform_mantel(SOC_C, Rootexudate_Level1),
  SOC_C_Level2 = perform_mantel(SOC_C, Rootexudate_Level2),
  SOC_C_Level3 = perform_mantel(SOC_C, Rootexudate_Level3)
)
mantel_results
mantel_df <- data.frame(
  Comparison = names(mantel_results),
  Statistic = sapply(mantel_results, function(x) x$statistic),
  P_value = sapply(mantel_results, function(x) x$signif)
)
write.xlsx(mantel_df, file = "Mantel_Results.xlsx", rowNames = FALSE)
mantel_plot <- ggplot(mantel_df, aes(x = Comparison, y = Statistic)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = sprintf("P=%.3f", P_value)), vjust = -0.5, size = 3.5) + 
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Mantel Test Results",
    x = "Comparison",
    y = "Mantel Statistic"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(mantel_plot)


######################## Relationship between microbial function and MWD
library(vegan)
library(ggplot2)
library(openxlsx)
mwd_data <- read.csv("MWD.csv", row.names = 1)
fapro_a_data <- read.xlsx("picrust2_KO_A.xlsx", sheet = 1, rowNames = TRUE)
fapro_b_data <- read.xlsx("picrust2_KO_B.xlsx", sheet = 1, rowNames = TRUE)
fapro_c_data <- read.xlsx("picrust2_KO_C.xlsx", sheet = 1, rowNames = TRUE)
# CPCoA
compute_pcoa <- function(data) {
  dist_matrix <- vegdist(t(data), method = "bray")
  pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2) 
  CPCoA1 <- pcoa_result$points[, 1]
  CPCoA2 <- pcoa_result$points[, 2]
  return(list(CPCoA1 = CPCoA1, CPCoA2 = CPCoA2, DistMatrix = dist_matrix))
}
pcoa_a <- compute_pcoa(fapro_a_data)
pcoa_b <- compute_pcoa(fapro_b_data)
pcoa_c <- compute_pcoa(fapro_c_data)
# Shannon and Richness
compute_indices <- function(data) {
  otu_data <- t(data)
  shannon_index <- diversity(otu_data, index = "shannon")  # Shannon
  richness_index <- rowSums(otu_data > 0)  # Richness
  return(list(Shannon = shannon_index, Richness = richness_index))
}
# Shannon and Richness
indices_a <- compute_indices(fapro_a_data)
indices_b <- compute_indices(fapro_b_data)
indices_c <- compute_indices(fapro_c_data)
# Mantel
compute_mantel <- function(mwd_dist, func_dist, label) {
  mantel_result <- mantel(mwd_dist, func_dist, method = "spearman", permutations = 999)
  return(data.frame(Group = label,
                    R_Value = mantel_result$statistic,
                    P_Value = mantel_result$signif))
}
# bray
mwd_dist <- vegdist(mwd_data, method = "bray")
# Mantel
mantel_results <- rbind(
  compute_mantel(mwd_dist, pcoa_a$DistMatrix, "MWD vs KO_A"),
  compute_mantel(mwd_dist, pcoa_b$DistMatrix, "MWD vs KO_B"),
  compute_mantel(mwd_dist, pcoa_c$DistMatrix, "MWD vs KO_C")
)
# Spearman
compute_spearman <- function(var1, var2, label) {
  correlation <- cor.test(var1, var2, method = "spearman")
  return(data.frame(Group = label,
                    R_Value = correlation$estimate,
                    P_Value = correlation$p.value))
}

spearman_results <- rbind(
  compute_spearman(pcoa_a$CPCoA1, rowMeans(mwd_data), "CPCoA1_A"),
  compute_spearman(pcoa_a$CPCoA2, rowMeans(mwd_data), "CPCoA2_A"),
  compute_spearman(indices_a$Shannon, rowMeans(mwd_data), "Shannon_A"),
  compute_spearman(indices_a$Richness, rowMeans(mwd_data), "Richness_A"),
  
  compute_spearman(pcoa_b$CPCoA1, rowMeans(mwd_data), "CPCoA1_B"),
  compute_spearman(pcoa_b$CPCoA2, rowMeans(mwd_data), "CPCoA2_B"),
  compute_spearman(indices_b$Shannon, rowMeans(mwd_data), "Shannon_B"),
  compute_spearman(indices_b$Richness, rowMeans(mwd_data), "Richness_B"),
  
  compute_spearman(pcoa_c$CPCoA1, rowMeans(mwd_data), "CPCoA1_C"),
  compute_spearman(pcoa_c$CPCoA2, rowMeans(mwd_data), "CPCoA2_C"),
  compute_spearman(indices_c$Shannon, rowMeans(mwd_data), "Shannon_C"),
  compute_spearman(indices_c$Richness, rowMeans(mwd_data), "Richness_C")
)
# Merge
all_results <- rbind(
  data.frame(mantel_results, Type = "Mantel"),
  data.frame(spearman_results, Type = "Spearman")
)
# Plot
p <- ggplot(all_results, aes(x = Group, y = R_Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = paste0("P=", round(P_Value, 3))), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
  labs(title = "R and P Values for Mantel and Spearman Tests",
       x = "Comparison Groups",
       y = "R Value") +
  scale_fill_manual(values = c("Mantel" = "#3f88c5", "Spearman" = "#f3ac66")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
print(p)
write.csv(all_results, "All_Results.csv", row.names = FALSE)


################## Relationship between microbial functions and exudate
level1_data <- read.xlsx("Rootexudate_Level1.xlsx", sheet = 1, rowNames = TRUE)
level2_data <- read.xlsx("Rootexudate_Level2.xlsx", sheet = 1, rowNames = TRUE)
level3_data <- read.xlsx("Rootexudate_Level3.xlsx", sheet = 1, rowNames = TRUE)
ko_a_data <- read.xlsx("picrust2_KO_A.xlsx", sheet = 1, rowNames = TRUE)
ko_b_data <- read.xlsx("picrust2_KO_B.xlsx", sheet = 1, rowNames = TRUE)
ko_c_data <- read.xlsx("picrust2_KO_C.xlsx", sheet = 1, rowNames = TRUE)
process_data <- function(data) {
  data <- t(data)
  data[is.na(data)] <- 0
  return(data)
}
level1_data <- process_data(level1_data)
level2_data <- process_data(level2_data)
level3_data <- process_data(level3_data)
ko_a_data <- process_data(ko_a_data)
ko_b_data <- process_data(ko_b_data)
ko_c_data <- process_data(ko_c_data)

# Root Exudate bray
level1_dist <- vegdist(level1_data, method = "bray")
level2_dist <- vegdist(level2_data, method = "bray")
level3_dist <- vegdist(level3_data, method = "bray")
# Picrust2 bray
ko_a_dist <- vegdist(ko_a_data, method = "bray")
ko_b_dist <- vegdist(ko_b_data, method = "bray")
ko_c_dist <- vegdist(ko_c_data, method = "bray")
# Mantel
compute_mantel <- function(exudate_dist, ko_dist, label) {
  mantel_result <- mantel(exudate_dist, ko_dist, method = "spearman", permutations = 999)
  return(data.frame(Group = label,
                    Mantel_R = mantel_result$statistic,
                    P_Value = mantel_result$signif))
}
# Mantel
results <- rbind(
  compute_mantel(level1_dist, ko_a_dist, "Level1 vs KO_A"),
  compute_mantel(level1_dist, ko_b_dist, "Level1 vs KO_B"),
  compute_mantel(level1_dist, ko_c_dist, "Level1 vs KO_C"),
  
  compute_mantel(level2_dist, ko_a_dist, "Level2 vs KO_A"),
  compute_mantel(level2_dist, ko_b_dist, "Level2 vs KO_B"),
  compute_mantel(level2_dist, ko_c_dist, "Level2 vs KO_C"),
  
  compute_mantel(level3_dist, ko_a_dist, "Level3 vs KO_A"),
  compute_mantel(level3_dist, ko_b_dist, "Level3 vs KO_B"),
  compute_mantel(level3_dist, ko_c_dist, "Level3 vs KO_C")
)

# Output
write.xlsx(results, file = "Mantel_Test_Results_Picrust2_RootExudate.xlsx")
# Plot
p <- ggplot(results, aes(x = Group, y = Mantel_R, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0("P=", round(P_Value, 3))), vjust = -0.5, size = 3.5) +
  labs(title = "Mantel Test Results for Picrust2 and Root Exudates",
       y = "Mantel Statistic R", x = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
print(p)

############################### Cohesion
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}
# Input data
b <- read.csv("ACK.csv", header = TRUE, row.names = 1, sep = ",")
print("Initial matrix dimensions:")
print(dim(b))
print("First few rows of data:")
print(head(b))
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]
print("Matrix dimensions after removing empty rows and columns:")
print(dim(c))
zero.cutoff <- ceiling(0.1 * dim(c)[1])
print(paste("Zero cutoff value:", zero.cutoff))
filtered_columns <- apply(c, 2, zero) < (dim(c)[1] - zero.cutoff)
print("Filtered columns (TRUE means retained):")
print(filtered_columns)
d <- c[, filtered_columns]
print("Matrix dimensions after filtering columns:")
print(dim(d))
if(ncol(d) == 0 || nrow(d) == 0) stop("Error: Filtered matrix 'd' is empty.")
d <- d[rowSums(d) > 0, ]
print("Matrix dimensions after filtering rows:")
print(dim(d))
rel.d <- d
print("Relative abundance matrix dimensions:")
print(dim(rel.d))
cor.mat.true <- cor(rel.d)
print("Correlation matrix dimensions:")
print(dim(cor.mat.true))
# Null model
use.custom.cors <- FALSE
tax.shuffle <- TRUE
iter <- 200
med.tax.cors <- vector()
if(use.custom.cors == FALSE) {
  for(which.taxon in 1:dim(rel.d)[2]){
    print(paste("Processing taxon:", which.taxon))
    perm.cor.vec.mat <- vector()
    for(i in 1:iter){
      perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
      rownames(perm.rel.d) <- rownames(rel.d)
      colnames(perm.rel.d) <- colnames(rel.d)
      for(j in 1:dim(rel.d)[2]){ 
        perm.rel.d[, j] <- sample(rel.d[, j])
      }
      perm.rel.d[, which.taxon] <- rel.d[, which.taxon]
      cor.mat.null <- cor(perm.rel.d)
      perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
    }
    if(ncol(perm.cor.vec.mat) > 0) {
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
    } else {
      stop(paste("Error: perm.cor.vec.mat is empty for taxon", which.taxon))
    }
    if(which.taxon %% 20 == 0) {print(which.taxon)}
  }
}
if(use.custom.cors) {
  obs.exp.cors.mat <- custom.cor.mat.sub
} else {
  obs.exp.cors.mat <- cor.mat.true - med.tax.cors
}
diag(obs.exp.cors.mat) <- 0
# connectedness and cohesion
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg
# Output
output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness", "Positive Connectedness", "Negative Cohesion", "Positive Cohesion")
print(output)
max_length <- max(length(output[[1]]), length(output[[2]]), length(output[[3]]), length(output[[4]]))
output_aligned <- lapply(output, function(x) {
  c(x, rep(NA, max_length - length(x))) 
})
output_df <- data.frame(
  Negative_Connectedness = output_aligned[[1]],
  Positive_Connectedness = output_aligned[[2]],
  Negative_Cohesion = output_aligned[[3]],
  Positive_Cohesion = output_aligned[[4]]
)
write.csv(output_df, file = "output_aligned.csv", row.names = FALSE)

#################################### Robustness
otutab <- read.table("ACK.txt", header = TRUE, row.names = 1, sep = "\t")
otutab[is.na(otutab)] <- 0
counts <- rowSums(otutab > 0)
otutab <- otutab[counts >= 1, ]
comm <- t(otutab)
cormatrix <- matrix(0, ncol(comm), ncol(comm))

for (i in 1:ncol(comm)) {
  for (j in i:ncol(comm)) {
    speciesi <- sapply(1:nrow(comm), function(k) {
      ifelse(comm[k, i] > 0, comm[k, i], ifelse(comm[k, j] > 0, 0.01, NA))
    })
    speciesj <- sapply(1:nrow(comm), function(k) {
      ifelse(comm[k, j] > 0, comm[k, j], ifelse(comm[k, i] > 0, 0.01, NA))
    })
    corij <- cor(log(speciesi)[!is.na(speciesi)], log(speciesj)[!is.na(speciesj)])
    cormatrix[i, j] <- cormatrix[j, i] <- corij
  }
}
row.names(cormatrix) <- colnames(cormatrix) <- colnames(comm)
cormatrix2 <- cormatrix * (abs(cormatrix) >= 0.80)
cormatrix2[is.na(cormatrix2)] <- 0
diag(cormatrix2) <- 0
print(sum(abs(cormatrix2) > 0) / 2)
print(sum(colSums(abs(cormatrix2)) > 0))
network.raw <- cormatrix2[colSums(abs(cormatrix2)) > 0, colSums(abs(cormatrix2)) > 0]
sp.ra2 <- colSums(comm)[colSums(abs(cormatrix2)) > 0]
stopifnot(sum(row.names(network.raw) == names(sp.ra2)) == length(sp.ra2))
cor.cutoff <- 0.3
cor.adj <- ifelse(abs(network.raw) >= cor.cutoff, 1, 0)
write.table(data.frame(cor.adj, check.names = FALSE), 'CK.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)
rand.remov.once <- function(netRaw, rm.percent, sp.ra, abundance.weighted = TRUE) {
  id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
  net.Raw <- netRaw 
  net.Raw[id.rm, ] <- 0
  net.Raw[, id.rm] <- 0
  if (abundance.weighted) {
    net.strength <- net.Raw * sp.ra
  } else {
    net.strength <- net.Raw
  }
  sp.meanInteraction <- colMeans(net.strength)
  id.rm2 <- which(sp.meanInteraction <= 0)
  remain.percent <- (nrow(netRaw) - length(id.rm2)) / nrow(netRaw)
  remain.percent
}
rmsimu <- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm = 100) {
  t(sapply(rm.p.list, function(x) {
    remains <- sapply(1:nperm, function(i) {
      rand.remov.once(netRaw = netRaw, rm.percent = x, sp.ra = sp.ra, abundance.weighted = abundance.weighted)
    })
    remain.mean <- mean(remains)
    remain.sd <- sd(remains)
    remain.se <- sd(remains) / (nperm^0.5)
    result <- c(remain.mean, remain.sd, remain.se)
    names(result) <- c("remain.mean", "remain.sd", "remain.se")
    result
  }))
}
Weighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = TRUE, nperm = 100)
Unweighted.simu <- rmsimu(netRaw = network.raw, rm.p.list = seq(0.05, 1, by = 0.05), sp.ra = sp.ra2, abundance.weighted = FALSE, nperm = 100)
dat1 <- data.frame(
  Proportion.removed = rep(seq(0.05, 1, by = 0.05), 2),
  rbind(Weighted.simu, Unweighted.simu),
  weighted = rep(c("weighted", "unweighted"), each = 20),
  year = rep(2014, 20),
  treat = rep("CK", 20)
)
write.csv(dat1, "random_removal_result_CK.csv")

################## Vulnerability
library(igraph)
source("/Users/chong/Library/CloudStorage/OneDrive-Personal/1 未完成的论文/The effects of microbial inoculants on soil aggregate/Analysis/2. Greenhouse experiment/Amplicon analysis/Network(不同粒径下CKTR比）/AggregateA/CK/Vulnerability/info.centrality.R")
setwd("/Users/chong/Library/CloudStorage/OneDrive-Personal/1 未完成的论文/The effects of microbial inoculants on soil aggregate/Analysis/2. Greenhouse experiment/Amplicon analysis/Network(不同粒径下CKTR比）/AggregateA/CK/Vulnerability")
#### get graph ####
## construct a graph from OTU table
otutab<-read.table("ACK.txt",header = T,row.names=1,sep="\t")
otutab[is.na(otutab)]<-0
comm<-t(otutab)
cormatrix=matrix(0,ncol(comm),ncol(comm))
for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i]>0,comm[k,i],ifelse(comm[k,j]>0,0.01,NA))
    })
    speciesj<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j]>0,comm[k,j],ifelse(comm[k,i]>0,0.01,NA))
    })
    corij<-cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    cormatrix[i,j]<-cormatrix[j,i]<-corij
  }}
row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.
cormatrix2<-cormatrix*(abs(cormatrix)>=0.80)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self   
cormatrix2[abs(cormatrix2)>0]<-1 # adjacency matrix
g = graph_from_adjacency_matrix(as.matrix(cormatrix2), mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL) # note: this graph contains isolated nodes.
## construct a graph from OTU table
# read in the correlation matrix downloaded from MENAP to construct the graph
# cormatrix.input  <- matrix(0, 1596, 1596)
# cormatrix.input[row(cormatrix.input) >= col(cormatrix.input)] <- scan("CKtra.txt")
# cormatrix <- t(cormatrix.input)
# cormatrix[abs(cormatrix)<0.8]<-0
# cormatrix[abs(cormatrix)>0]<-1 # adjacency matrix
# g = graph_from_adjacency_matrix(as.matrix(cormatrix), mode="upper", weighted = NULL, diag = FALSE, add.colnames = NULL) # note: this graph contains isolated nodes.
# ## End 3) read in the correlation matrix downloaded from MENAP to construct the graph
# ### End get graph
# remove isolated nodes
iso_node_id = which(degree(g)==0)
g2 = delete.vertices(g, iso_node_id) # graph without isolated nodes
#check node number and links
length(V(g2));length(E(g2))   
# calculate vulnerability of each node
node.vul<-info.centrality.vertex(g2)
max(node.vul)



############################################################################# Incubation experiment

########################## CPCoA
cpcoa <- function (otutab, metadata, dis = "bray", groupID = "Group", 
                   ellipse = T, label = F) 
{
  p_list = c("ggplot2", "vegan", "ggrepel")
  for (p in p_list) {
    if (!requireNamespace(p)) {
      install.packages(p)
    }
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
  idx = rownames(metadata) %in% colnames(otutab)
  metadata = metadata[idx, , drop = F]
  otutab = otutab[, rownames(metadata)]
  sampFile = as.data.frame(metadata[, groupID], row.names = row.names(metadata))
  colnames(sampFile)[1] = "group"
  if (length(unique(sampFile$group)) > 2) {
    variability_table = function(cca) {
      chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
      variability_table = cbind(chi, chi/chi[1])
      colnames(variability_table) = c("inertia", "proportion")
      rownames(variability_table) = c("total", "constrained", 
                                      "unconstrained")
      return(variability_table)
    }
    capscale.gen = capscale(t(otutab) ~ group, data = sampFile, 
                            add = F, sqrt.dist = T, distance = dis)
    set.seed(123)
    perm_anova.gen = anova.cca(capscale.gen, permutations = 1000, 
                               parallel = 4)
    var_tbl.gen = variability_table(capscale.gen)
    eig = capscale.gen$CCA$eig
    variance = var_tbl.gen["constrained", "proportion"]
    p.val = perm_anova.gen[1, 4]
    points = as.data.frame(capscale.gen$CCA$wa)
    points = cbind(sampFile, points[rownames(points), ])
    
  } 
  return(list(eig = eig, points = points, var = variance, p = p.val))
}
library(ggplot2)
library(vegan)
library(ape)
library(plyr)
library(openxlsx)
otu <- read.csv('genes.csv', row.names = 1)
metadata <- read.csv('design.csv')
colnames(metadata) <- c("samples", "Treatment", "Group","Inter" )
otu <- otu[, metadata$samples]
metadata <- metadata[metadata$samples %in% colnames(otu), ]
sampFile <- as.data.frame(metadata$Treatment, row.names = metadata$samples)
colnames(sampFile)[1] <- "treatment"
cpcoa_result <- capscale(t(otu) ~ treatment, data = sampFile, add = FALSE, sqrt.dist = TRUE, distance = "bray")
cppoints <- as.data.frame(cpcoa_result$CCA$wa)
cppoints$samples <- row.names(cppoints)
colnames(cppoints) <- c("CPCoA1", "CPCoA2", "CPCoA3", "CPCoA4", "samples")
cppoints <- cppoints[, c("CPCoA1", "CPCoA2", "samples")]
eigenvalues <- cpcoa_result$CCA$eig
pca1 <- round(100 * eigenvalues[1] / sum(eigenvalues), 2)
pca2 <- round(100 * eigenvalues[2] / sum(eigenvalues), 2)
variability_table <- function(cca) {
  chi <- c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi/chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained")
  return(variability_table)
}
var_tbl <- variability_table(cpcoa_result)
variance <- var_tbl["constrained", "proportion"]
p_val <- anova.cca(cpcoa_result, permutations = 1000)[1, 4]
df <- merge(cppoints, metadata, by = "samples")
color <- c('black','black',"#F7AF34", "#448DCD", "#ffcd85", "#a1b8d5", "grey",'black')
shapes <- c(17, 15, 16, 18, 19, 20, 21, 22)
# Plot
p1 <- ggplot(df, aes(CPCoA1, CPCoA2, shape = Group)) +
  geom_point(aes(color = Treatment), size = 5) +  
  xlab(paste("CPCoA1 (", pca1, "%)", sep = "")) + 
  ylab(paste("CPCoA2 (", pca2, "%)", sep = "")) +
  stat_ellipse(aes(group = Inter, color = Inter), level = 0.95, linetype = "dashed", linewidth = 1) + 
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = color) + 
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black', size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(paste("Variance explained: ", round(variance * 100, 2), "%; P = ", format(p_val, digits = 2), sep = ""))

p1
pdf("CPCoA.pdf", width = 7, height = 6)
print(p1)
dev.off()
coordinates <- df[, c("samples", "CPCoA1", "CPCoA2")]
print(coordinates)
# Saved
write.csv(coordinates, "CPCoA_Coordinates.csv", row.names = FALSE)

########################## PERMANOVA
library(vegan)
otutab <- read.csv("otutab.csv", row.names = 1)
metadata <- read.csv("treatment.csv")
rownames(metadata) <- metadata$Treatment 
otutab <- otutab[, rownames(metadata)] 
# Bray-Curtis
dist_matrix <- vegdist(t(otutab), method = "bray")
# 4. PERMANOVA
# Year、Duration、Inoculant and Aggregate 
perm_model_interaction <- adonis2(dist_matrix ~ (Year + Duration + Inoculant + Aggregate)^2, 
                                  data = metadata, 
                                  strata = metadata$Random,
                                  permutations = 999)
# PERMANOVA result
print(perm_model_interaction)
# Output
perm_results <- as.data.frame(perm_model_interaction)
# Saved
write.csv(perm_results, "permanova_results.csv", row.names = TRUE)
write.table(perm_results, "permanova_results.txt", sep = "\t", row.names = TRUE, quote = FALSE)

################# Relationship between Alpha diversities and SOC
library(corrplot)
library(ggplot2)
library(reshape2)
data <- read.csv("Correlation.csv")
group_1 <- c("Bacterialrichness_A", "Bacterialchao1_A", "BacterialACE_A", "Bacterialshannon_A",
             "Fungalrichness_A", "Fungalchao1_A", "FungalACE_A", "Fungalshannon_A", 
             "Bacterialrichness_B", "Bacterialchao1_B", "BacterialACE_B", "Bacterialshannon_B", 
             "Fungalrichness_B", "Fungalchao1_B", "FungalACE_B", "Fungalshannon_B", 
             "Bacterialrichness_C", "Bacterialchao1_C", "BacterialACE_C", "Bacterialshannon_C", 
             "Fungalrichness_C", "Fungalchao1_C", "FungalACE_C", "Fungalshannon_C", 
             "Bacterialrichness_D", "Bacterialchao1_D", "BacterialACE_D", "Bacterialshannon_D", 
             "Fungalrichness_D", "Fungalchao1_D", "FungalACE_D", "Fungalshannon_D")

group_2 <- c("POC", "MAOC", "Aggregate_A", "Aggregate_B", "Aggregate_C", "Aggregate_D", "MWD")
result_cor <- matrix(NA, nrow = length(group_1), ncol = length(group_2))
result_p <- matrix(NA, nrow = length(group_1), ncol = length(group_2))
for (i in 1:length(group_1)) {
  for (j in 1:length(group_2)) {
    if (!is.null(data[[group_1[i]]]) && !is.null(data[[group_2[j]]])) {
      test <- cor.test(data[[group_1[i]]], data[[group_2[j]]], method = "spearman", use = "complete.obs")
      result_cor[i, j] <- test$estimate
      result_p[i, j] <- test$p.value
    }
  }
}
# BH adjusted
adjusted_p <- p.adjust(as.vector(result_p), method = "BH")
result_p_adjusted <- matrix(adjusted_p, nrow = length(group_1), ncol = length(group_2))
rownames(result_cor) <- group_1
colnames(result_cor) <- group_2
rownames(result_p) <- group_1
colnames(result_p) <- group_2
rownames(result_p_adjusted) <- group_1
colnames(result_p_adjusted) <- group_2
write.csv(result_cor, "correlation_results.csv")
write.csv(result_p, "raw_p_values.csv")
write.csv(result_p_adjusted, "adjusted_p_values.csv")
cor_data <- melt(result_cor)
colnames(cor_data) <- c("Group1", "Group2", "Correlation")
p_data <- melt(result_p_adjusted)
colnames(p_data) <- c("Group1", "Group2", "AdjustedP")
plot_data <- merge(cor_data, p_data, by = c("Group1", "Group2"))
plot_data$Significance <- ifelse(plot_data$AdjustedP < 0.05, "*", "")
# Plot
ggplot(plot_data, aes(x = Group2, y = Group1, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = Significance), color = "black", size = 3) +
  scale_fill_gradient2(low = "#F7AF34", high = "#448DCD", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle("Correlation Heatmap with BH-adjusted Significance")

############################### Relationship between microbial community and MAOC
library(vegan)
bacteria_a <- read.csv("Bacteria_A.csv", row.names = 1)
bacteria_b <- read.csv("Bacteria_B.csv", row.names = 1)
bacteria_c <- read.csv("Bacteria_C.csv", row.names = 1)
bacteria_d <- read.csv("Bacteria_D.csv", row.names = 1)
fungi_a <- read.csv("Fungi_A.csv", row.names = 1)
fungi_b <- read.csv("Fungi_B.csv", row.names = 1)
fungi_c <- read.csv("Fungi_C.csv", row.names = 1)
fungi_d <- read.csv("Fungi_D.csv", row.names = 1)
data <- read.csv("SOCMWD.csv", row.names = 1)
# Bray-Curtis
calculate_bray_curtis <- function(data) {
  vegdist(t(data), method = "bray")
}
bacteria_a_dist <- calculate_bray_curtis(bacteria_a)
bacteria_b_dist <- calculate_bray_curtis(bacteria_b)
bacteria_c_dist <- calculate_bray_curtis(bacteria_c)
bacteria_d_dist <- calculate_bray_curtis(bacteria_d)
fungi_a_dist <- calculate_bray_curtis(fungi_a)
fungi_b_dist <- calculate_bray_curtis(fungi_b)
fungi_c_dist <- calculate_bray_curtis(fungi_c)
fungi_d_dist <- calculate_bray_curtis(fungi_d)
MWD_dist <- dist(data$MWD, method = "euclidean")
# Mantel
MWD_mantel_bacteria_a <- mantel(bacteria_a_dist, MWD_dist, method = "spearman", permutations = 999)
MWD_mantel_bacteria_b <- mantel(bacteria_b_dist, MWD_dist, method = "spearman", permutations = 999)
MWD_mantel_bacteria_c <- mantel(bacteria_c_dist, MWD_dist, method = "spearman", permutations = 999)
MWD_mantel_bacteria_d <- mantel(bacteria_d_dist, MWD_dist, method = "spearman", permutations = 999)
MWD_mantel_fungi_a <- mantel(fungi_a_dist, MWD_dist, method = "spearman", permutations = 999)
MWD_mantel_fungi_b <- mantel(fungi_b_dist, MWD_dist, method = "spearman", permutations = 999)
MWD_mantel_fungi_c <- mantel(fungi_c_dist, MWD_dist, method = "spearman", permutations = 999)
MWD_mantel_fungi_d <- mantel(fungi_d_dist, MWD_dist, method = "spearman", permutations = 999)
# Merge
MWD_mantel_results <- data.frame(
  c("Bacteria_A", "Bacteria_B", "Bacteria_C", "Bacteria_D", "Fungi_A", "Fungi_B", "Fungi_C", "Fungi_D"),
  Mantel_R = c(
    MWD_mantel_bacteria_a$statistic,
    MWD_mantel_bacteria_b$statistic,
    MWD_mantel_bacteria_c$statistic,
    MWD_mantel_bacteria_d$statistic,
    MWD_mantel_fungi_a$statistic,
    MWD_mantel_fungi_b$statistic,
    MWD_mantel_fungi_c$statistic,
    MWD_mantel_fungi_d$statistic
  ),
  P_Value = c(
    MWD_mantel_bacteria_a$signif,
    MWD_mantel_bacteria_b$signif,
    MWD_mantel_bacteria_c$signif,
    MWD_mantel_bacteria_d$signif,
    MWD_mantel_fungi_a$signif,
    MWD_mantel_fungi_b$signif,
    MWD_mantel_fungi_c$signif,
    MWD_mantel_fungi_d$signif
  )
)
print(MWD_mantel_results)
# Saved
write.csv(MWD_mantel_results, "MWD_Mantel_Test_Results.csv", row.names = FALSE)

# MAOC euclidean
MAOC_dist <- dist(data$MAOC, method = "euclidean")
# Mantel
MAOC_mantel_bacteria_a <- mantel(bacteria_a_dist, MAOC_dist, method = "spearman", permutations = 999)
MAOC_mantel_bacteria_b <- mantel(bacteria_b_dist, MAOC_dist, method = "spearman", permutations = 999)
MAOC_mantel_bacteria_c <- mantel(bacteria_c_dist, MAOC_dist, method = "spearman", permutations = 999)
MAOC_mantel_bacteria_d <- mantel(bacteria_d_dist, MAOC_dist, method = "spearman", permutations = 999)
MAOC_mantel_fungi_a <- mantel(fungi_a_dist, MAOC_dist, method = "spearman", permutations = 999)
MAOC_mantel_fungi_b <- mantel(fungi_b_dist, MAOC_dist, method = "spearman", permutations = 999)
MAOC_mantel_fungi_c <- mantel(fungi_c_dist, MAOC_dist, method = "spearman", permutations = 999)
MAOC_mantel_fungi_d <- mantel(fungi_d_dist, MAOC_dist, method = "spearman", permutations = 999)
# Merge
MAOC_mantel_results <- data.frame(
  c("Bacteria_A", "Bacteria_B", "Bacteria_C", "Bacteria_D", "Fungi_A", "Fungi_B", "Fungi_C", "Fungi_D"),
  Mantel_R = c(
    MAOC_mantel_bacteria_a$statistic,
    MAOC_mantel_bacteria_b$statistic,
    MAOC_mantel_bacteria_c$statistic,
    MAOC_mantel_bacteria_d$statistic,
    MAOC_mantel_fungi_a$statistic,
    MAOC_mantel_fungi_b$statistic,
    MAOC_mantel_fungi_c$statistic,
    MAOC_mantel_fungi_d$statistic
  ),
  P_Value = c(
    MAOC_mantel_bacteria_a$signif,
    MAOC_mantel_bacteria_b$signif,
    MAOC_mantel_bacteria_c$signif,
    MAOC_mantel_bacteria_d$signif,
    MAOC_mantel_fungi_a$signif,
    MAOC_mantel_fungi_b$signif,
    MAOC_mantel_fungi_c$signif,
    MAOC_mantel_fungi_d$signif
  )
)
print(MAOC_mantel_results)
# Saved
write.csv(MAOC_mantel_results, "MAOC_Mantel_Test_Results.csv", row.names = FALSE)

# POC euclidean
POC_dist <- dist(data$POC, method = "euclidean")
# Mantel
POC_mantel_bacteria_a <- mantel(bacteria_a_dist, POC_dist, method = "spearman", permutations = 999)
POC_mantel_bacteria_b <- mantel(bacteria_b_dist, POC_dist, method = "spearman", permutations = 999)
POC_mantel_bacteria_c <- mantel(bacteria_c_dist, POC_dist, method = "spearman", permutations = 999)
POC_mantel_bacteria_d <- mantel(bacteria_d_dist, POC_dist, method = "spearman", permutations = 999)
POC_mantel_fungi_a <- mantel(fungi_a_dist, POC_dist, method = "spearman", permutations = 999)
POC_mantel_fungi_b <- mantel(fungi_b_dist, POC_dist, method = "spearman", permutations = 999)
POC_mantel_fungi_c <- mantel(fungi_c_dist, POC_dist, method = "spearman", permutations = 999)
POC_mantel_fungi_d <- mantel(fungi_d_dist, POC_dist, method = "spearman", permutations = 999)
# Merge
POC_mantel_results <- data.frame(
  c("Bacteria_A", "Bacteria_B", "Bacteria_C", "Bacteria_D", "Fungi_A", "Fungi_B", "Fungi_C", "Fungi_D"),
  Mantel_R = c(
    POC_mantel_bacteria_a$statistic,
    POC_mantel_bacteria_b$statistic,
    POC_mantel_bacteria_c$statistic,
    POC_mantel_bacteria_d$statistic,
    POC_mantel_fungi_a$statistic,
    POC_mantel_fungi_b$statistic,
    POC_mantel_fungi_c$statistic,
    POC_mantel_fungi_d$statistic
  ),
  P_Value = c(
    POC_mantel_bacteria_a$signif,
    POC_mantel_bacteria_b$signif,
    POC_mantel_bacteria_c$signif,
    POC_mantel_bacteria_d$signif,
    POC_mantel_fungi_a$signif,
    POC_mantel_fungi_b$signif,
    POC_mantel_fungi_c$signif,
    POC_mantel_fungi_d$signif
  )
)
print(POC_mantel_results)
# Saved
write.csv(POC_mantel_results, "POC_Mantel_Test_Results.csv", row.names = FALSE)

MWD_mantel_results$Variable <- "MWD"
MAOC_mantel_results$Variable <- "MAOC"
POC_mantel_results$Variable <- "POC"
all_results <- rbind(MWD_mantel_results, MAOC_mantel_results, POC_mantel_results)
write.csv(all_results, "all_results.csv", row.names = FALSE)
all_results <- read.csv("all_results.csv")
# Plot
library(ggplot2)
ggplot(all_results, aes(x = Name, y = Mantel_R, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mantel R Values for Different Communities",
       x = "Community",
       y = "Mantel R Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############################################################# GLM
library(car)
library(fitdistrplus)
library(e1071)
library(lmtest)
tbl.GLM <- read.csv("GLM.csv")

#################################  BacterialRichness_model  
table(tbl.GLM$BacterialRichness)
### test mean 6685.653
mean(tbl.GLM$BacterialRichness)
### test variance 1898683
var(tbl.GLM$BacterialRichness)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonBacterialRichness <- glm(BacterialRichness ~ Inoculants * Aggregate * Duration * Year, 
                                    data = tbl.GLM, family = poisson)
### The discrete factor was calculated  8.768145
deviance(glm_poissonBacterialRichness) / df.residual(glm_poissonBacterialRichness)
### Distribution of response variables
hist(tbl.GLM$BacterialRichness, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$BacterialRichness, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiBacterialRichness <- glm(BacterialRichness ~ Inoculants * Aggregate * Duration * Year, 
                                  data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiBacterialRichness), residuals(glm_quasiBacterialRichness), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiBacterialRichness, type = "II")

glm_quasiBacterialRichness_simplified <- glm(BacterialRichness ~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiBacterialRichness_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacterialRichness
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants                3.8  1  0.0515361 .  
# Aggregate                21.0  3  0.0001053 ***
# Duration              10604.1  3  < 2.2e-16 ***
# Year                     18.2  1  2.017e-05 ***
# Inoculants:Aggregate      1.1  3  0.7823108    
# Inoculants:Duration      27.5  3  4.734e-06 ***
# Aggregate:Duration        8.7  9  0.4606989    
# Inoculants:Year           2.7  1  0.0989049 .  
# Aggregate:Year            0.6  3  0.8899464    
# Duration:Year            54.9  3  7.261e-12 ***



################################   BacterialShannon_model  
table(tbl.GLM$BacterialShannon)
### test mean 6.817174
mean(tbl.GLM$BacterialShannon)
### test variance 0.1258915
var(tbl.GLM$BacterialShannon)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$BacterialShannon, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$BacterialShannon) ## 0.9616302
kurtosis(tbl.GLM$BacterialShannon) ## -0.4421884
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianBacterialShannon <- glm(BacterialShannon ~ Inoculants * Aggregate * Duration * Year, 
                                    data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianBacterialShannon))
qqline(residuals(glm_gaussianBacterialShannon), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianBacterialShannon))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianBacterialShannon)
# W = 0.95297, p-value = 3.579e-09

## Test Homoscedasticity
plot(fitted(glm_gaussianBacterialShannon), residuals(glm_gaussianBacterialShannon), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianBacterialShannon)
# studentized Breusch-Pagan test
# data:  glm_gaussianBacterialShannon
# BP = 111.01, df = 63, p-value = 0.0001809

#### log transformation
tbl.GLM$BacterialShannon_log <- log(tbl.GLM$BacterialShannon)
hist(tbl.GLM$BacterialShannon_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedBacterialShannon_log <- glm(BacterialShannon_log ~ Inoculants * Aggregate * Duration * Year, 
                                           data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedBacterialShannon_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedBacterialShannon_log)
# W = 0.95451, p-value = 5.706e-09
## Test Homoscedasticity
plot(fitted(glm_transformedBacterialShannon_log), residuals(glm_transformedBacterialShannon_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedBacterialShannon_log)
# studentized Breusch-Pagan test
# data:  glm_transformedBacterialShannon_log
# BP = 104.76, df = 63, p-value = 0.0007505

#### sqrt transformation
tbl.GLM$BacterialShannon_sqrt <- sqrt(tbl.GLM$BacterialShannon)
hist(tbl.GLM$BacterialShannon_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedBacterialShannon_sqrt <- glm(BacterialShannon_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                            data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedBacterialShannon_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedBacterialShannon_sqrt)
# W = 0.954, p-value = 4.887e-09
## Test Homoscedasticity
plot(fitted(glm_transformedBacterialShannon_sqrt), residuals(glm_transformedBacterialShannon_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedBacterialShannon_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedBacterialShannon_sqrt
# BP = 107.74, df = 63, p-value = 0.0003853

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianBacterialShannon)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -557.1599  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedBacterialShannon_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -1917.831  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedBacterialShannon_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -1726.008  
### Choose log

###### Data transformation (log transformation)
tbl.GLM$BacterialShannon_log <- log(tbl.GLM$BacterialShannon)
glm_transformedBacterialShannon <- glm(BacterialShannon_log ~ Inoculants * Aggregate * Duration * Year, 
                                       data = tbl.GLM, family = gaussian)
Anova(glm_transformedBacterialShannon, type = "II") 

glm_transformedBacterialShannon_simplified <- glm(BacterialShannon_log ~ Inoculants + Aggregate + Duration + Year +
                                                    Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                    Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedBacterialShannon_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacterialShannon_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants                1.5  1   0.213802    
# Aggregate                71.7  3  1.824e-15 ***
# Duration               3930.0  3  < 2.2e-16 ***
# Year                     24.7  1  6.645e-07 ***
# Inoculants:Aggregate      2.7  3   0.445134    
# Inoculants:Duration      12.7  3   0.005455 ** 
# Aggregate:Duration       15.7  9   0.072431 .  
# Inoculants:Year           9.5  1   0.002087 ** 
# Aggregate:Year            2.1  3   0.553077    
# Duration:Year            64.3  3  7.091e-14 ***    



#################################  BacterialChao_model  
table(tbl.GLM$BacterialChao)
### test mean 9771.292
mean(tbl.GLM$BacterialChao)
### test variance 2224076
var(tbl.GLM$BacterialChao)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonBacterialChao <- glm(BacterialChao ~ Inoculants * Aggregate * Duration * Year, 
                                data = tbl.GLM, family = poisson)
### The discrete factor was calculated  23.80857
deviance(glm_poissonBacterialChao) / df.residual(glm_poissonBacterialChao)
### Distribution of response variables
hist(tbl.GLM$BacterialChao, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$BacterialChao, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiBacterialChao <- glm(BacterialChao ~ Inoculants * Aggregate * Duration * Year, 
                              data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiBacterialChao), residuals(glm_quasiBacterialChao), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiBacterialChao, type = "II")

glm_quasiBacterialChao_simplified <- glm(BacterialChao ~ Inoculants + Aggregate + Duration + Year +
                                           Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                           Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiBacterialChao_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacterialChao
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               1.35  1    0.24559    
# Aggregate                9.06  3    0.02851 *  
# Duration              2830.54  3  < 2.2e-16 ***
# Year                    17.47  1  2.918e-05 ***
# Inoculants:Aggregate     2.85  3    0.41580    
# Inoculants:Duration     26.98  3  5.952e-06 ***
# Aggregate:Duration       6.18  9    0.72206    
# Inoculants:Year          0.02  1    0.89003    
# Aggregate:Year           0.06  3    0.99613    
# Duration:Year           30.27  3  1.208e-06 ***


#################################  BacterialACE_model  
table(tbl.GLM$BacterialACE)
### test mean 10232.71
mean(tbl.GLM$BacterialACE)
### test variance 2215698
var(tbl.GLM$BacterialACE)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonBacterialACE <- glm(BacterialACE ~ Inoculants * Aggregate * Duration * Year, 
                               data = tbl.GLM, family = poisson)
### The discrete factor was calculated  26.26661
deviance(glm_poissonBacterialACE) / df.residual(glm_poissonBacterialACE)
### Distribution of response variables
hist(tbl.GLM$BacterialACE, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$BacterialACE, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiBacterialACE <- glm(BacterialACE ~ Inoculants * Aggregate * Duration * Year, 
                             data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiBacterialACE), residuals(glm_quasiBacterialACE), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiBacterialACE, type = "II")

glm_quasiBacterialACE_simplified <- glm(BacterialACE ~ Inoculants + Aggregate + Duration + Year +
                                          Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                          Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiBacterialACE_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: BacterialACE
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               0.62  1    0.43071    
# Aggregate                8.24  3    0.04136 *  
# Duration              2409.12  3  < 2.2e-16 ***
# Year                    15.67  1  7.541e-05 ***
# Inoculants:Aggregate     2.41  3    0.49117    
# Inoculants:Duration     26.83  3  6.388e-06 ***
# Aggregate:Duration       7.92  9    0.54261    
# Inoculants:Year          0.00  1    0.96588    
# Aggregate:Year           0.23  3    0.97245    
# Duration:Year           26.97  3  5.960e-06 ***   


#################################  FungalRichness_model  
table(tbl.GLM$FungalRichness)
### test mean 676.8409
mean(tbl.GLM$FungalRichness)
### test variance 39068.57
var(tbl.GLM$FungalRichness)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonFungalRichness <- glm(FungalRichness ~ Inoculants * Aggregate * Duration * Year, 
                                 data = tbl.GLM, family = poisson)
### The discrete factor was calculated  9.57342
deviance(glm_poissonFungalRichness) / df.residual(glm_poissonFungalRichness)
### Distribution of response variables
hist(tbl.GLM$FungalRichness, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$FungalRichness, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiFungalRichness <- glm(FungalRichness ~ Inoculants * Aggregate * Duration * Year, 
                               data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiFungalRichness), residuals(glm_quasiFungalRichness), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiFungalRichness, type = "II")

glm_quasiFungalRichness_simplified <- glm(FungalRichness ~ Inoculants + Aggregate + Duration + Year +
                                            Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                            Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiFungalRichness_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FungalRichness
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               0.13  1   0.722308    
# Aggregate               12.85  3   0.004981 ** 
# Duration              1601.17  3  < 2.2e-16 ***
# Year                     4.24  1   0.039580 *  
# Inoculants:Aggregate     1.90  3   0.592417    
# Inoculants:Duration     51.05  3  4.766e-11 ***
# Aggregate:Duration      12.68  9   0.177515    
# Inoculants:Year         16.43  1  5.060e-05 ***
# Aggregate:Year           5.61  3   0.131967    
# Duration:Year           58.48  3  1.243e-12 ***


#################################  FungalShannon_model  
table(tbl.GLM$FungalShannon)
### test mean 4.047026
mean(tbl.GLM$FungalShannon)
### test variance 0.3135073
var(tbl.GLM$FungalShannon)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$FungalShannon, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$FungalShannon) ## -1.44322
kurtosis(tbl.GLM$FungalShannon) ## 3.243447
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianFungalShannon <- glm(FungalShannon ~ Inoculants * Aggregate * Duration * Year, 
                                 data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianFungalShannon))
qqline(residuals(glm_gaussianFungalShannon), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianFungalShannon))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianFungalShannon)
# W = .9023, p-value = 2.814e-14, do not meet standard

## Test Homoscedasticity
plot(fitted(glm_gaussianFungalShannon), residuals(glm_gaussianFungalShannon), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianFungalShannon)
# studentized Breusch-Pagan test
# data:  glm_gaussianFungalShannon
# BP = 89.472, df = 63, p-value = 0.01583

#### log transformation
tbl.GLM$FungalShannon_log <- log(tbl.GLM$FungalShannon)
hist(tbl.GLM$FungalShannon_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedFungalShannon_log <- glm(FungalShannon_log ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedFungalShannon_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedFungalShannon_log)
# W = 0.81765, p-value < 2.2e-16
## Test Homoscedasticity
plot(fitted(glm_transformedFungalShannon_log), residuals(glm_transformedFungalShannon_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedFungalShannon_log)
# studentized Breusch-Pagan test
# data:  glm_transformedFungalShannon_log
# BP = 83.039, df = 63, p-value = 0.0462

#### sqrt transformation
tbl.GLM$FungalShannon_sqrt <- sqrt(tbl.GLM$FungalShannon)
hist(tbl.GLM$FungalShannon_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedFungalShannon_sqrt <- glm(FungalShannon_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                         data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedFungalShannon_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedFungalShannon_sqrt)
# W = 0.86624, p-value < 2.2e-16
## Test Homoscedasticity
plot(fitted(glm_transformedFungalShannon_sqrt), residuals(glm_transformedFungalShannon_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedFungalShannon_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedFungalShannon_sqrt
# BP = 85.798, df = 63, p-value = 0.02969

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianFungalShannon)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 551.1229  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedFungalShannon_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -280.8132  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedFungalShannon_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -362.4064  
### Choose log

###### Data transformation (sqrt transformation)
tbl.GLM$FungalShannon_sqrt <- sqrt(tbl.GLM$FungalShannon)
glm_transformedFungalShannon <- glm(FungalShannon_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                    data = tbl.GLM, family = gaussian)
Anova(glm_transformedFungalShannon, type = "II") 

glm_transformedFungalShannon_simplified <- glm(FungalShannon_sqrt ~ Inoculants + Aggregate + Duration + Year +
                                                 Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                 Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedFungalShannon_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FungalShannon_sqrt
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              1.672  1   0.195965    
# Aggregate              12.429  3   0.006050 ** 
# Duration               51.551  3  3.733e-11 ***
# Year                    5.393  1   0.020220 *  
# Inoculants:Aggregate    1.737  3   0.628774    
# Inoculants:Duration     4.396  3   0.221718    
# Aggregate:Duration     12.243  9   0.199980    
# Inoculants:Year         8.630  1   0.003306 ** 
# Aggregate:Year          1.159  3   0.762785    
# Duration:Year           5.984  3   0.112393   



#################################  FungalChao_model  
table(tbl.GLM$FungalChao)
### test mean 967.5156
mean(tbl.GLM$FungalChao)
### test variance 84546.6
var(tbl.GLM$FungalChao)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonFungalChao <- glm(FungalChao ~ Inoculants * Aggregate * Duration * Year, 
                             data = tbl.GLM, family = poisson)
### The discrete factor was calculated  15.82207
deviance(glm_poissonFungalChao) / df.residual(glm_poissonFungalChao)
### Distribution of response variables
hist(tbl.GLM$FungalChao, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$FungalChao, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiFungalChao <- glm(FungalChao ~ Inoculants * Aggregate * Duration * Year, 
                           data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiFungalChao), residuals(glm_quasiFungalChao), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiFungalChao, type = "II")

glm_quasiFungalChao_simplified <- glm(FungalChao ~ Inoculants + Aggregate + Duration + Year +
                                        Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                        Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiFungalChao_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FungalChao
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               0.16  1  0.6865965    
# Aggregate                8.59  3  0.0353280 *  
# Duration              1551.44  3  < 2.2e-16 ***
# Year                     7.67  1  0.0056081 ** 
# Inoculants:Aggregate     0.69  3  0.8753103    
# Inoculants:Duration     40.88  3  6.919e-09 ***
# Aggregate:Duration      17.24  9  0.0450460 *  
# Inoculants:Year         12.87  1  0.0003347 ***
# Aggregate:Year           2.09  3  0.5534441    
# Duration:Year           55.79  3  4.657e-12 ***    


#################################  FungalACE_model  
table(tbl.GLM$FungalACE)
### test mean 982.3985
mean(tbl.GLM$FungalACE)
### test variance 89753.89
var(tbl.GLM$FungalACE)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonFungalACE <- glm(FungalACE ~ Inoculants * Aggregate * Duration * Year, 
                            data = tbl.GLM, family = poisson)
### The discrete factor was calculated  14.85503
deviance(glm_poissonFungalACE) / df.residual(glm_poissonFungalACE)
### Distribution of response variables
hist(tbl.GLM$FungalACE, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$FungalACE, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiFungalACE <- glm(FungalACE ~ Inoculants * Aggregate * Duration * Year, 
                          data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiFungalACE), residuals(glm_quasiFungalACE), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiFungalACE, type = "II")

glm_quasiFungalACE_simplified <- glm(FungalACE ~ Inoculants + Aggregate + Duration + Year +
                                       Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                       Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiFungalACE_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: FungalACE
#                      LR Chisq Df Pr(>Chisq)
# Inoculants               0.03  1   0.868775
# Aggregate                9.54  3   0.022940 *
# Duration              1765.05  3  < 2.2e-16 ***
# Year                     9.93  1   0.001624 **
# Inoculants:Aggregate     1.18  3   0.758917
# Inoculants:Duration     41.27  3  5.741e-09 ***
# Aggregate:Duration      18.31  9   0.031702 *
# Inoculants:Year         16.51  1  4.833e-05 ***
# Aggregate:Year           2.61  3   0.455576
# Duration:Year           68.46  3  9.116e-15 ***


### Generalized Linear Model（GLM）
################################   PhylumAcidobacteria_model  
table(tbl.GLM$PhylumAcidobacteria)
### test mean 17.94148
mean(tbl.GLM$PhylumAcidobacteria)
### test variance 6.031494
var(tbl.GLM$PhylumAcidobacteria)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumAcidobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumAcidobacteria) ## 0.8228475
kurtosis(tbl.GLM$PhylumAcidobacteria) ## 1.120897
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPhylumAcidobacteria <- glm(PhylumAcidobacteria ~ Inoculants * Aggregate * Duration * Year, 
                                       data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPhylumAcidobacteria))
qqline(residuals(glm_gaussianPhylumAcidobacteria), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPhylumAcidobacteria))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPhylumAcidobacteria)
# W = 0.96694, p-value = 3.646e-07

## Test Homoscedasticity
plot(fitted(glm_gaussianPhylumAcidobacteria), residuals(glm_gaussianPhylumAcidobacteria), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPhylumAcidobacteria)
# studentized Breusch-Pagan test
# data:  glm_gaussianPhylumAcidobacteria
# BP = 99.282, df = 63, p-value = 0.00241

#### log transformation
tbl.GLM$PhylumAcidobacteria_log <- log(tbl.GLM$PhylumAcidobacteria)
hist(tbl.GLM$PhylumAcidobacteria_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumAcidobacteria_log <- glm(PhylumAcidobacteria_log ~ Inoculants * Aggregate * Duration * Year, 
                                              data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumAcidobacteria_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumAcidobacteria_log)
# W = 0.98206, p-value = 0.0002258
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumAcidobacteria_log), residuals(glm_transformedPhylumAcidobacteria_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumAcidobacteria_log)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumAcidobacteria_log
# BP = 96.648, df = 63, p-value = 0.004106

#### sqrt transformation
tbl.GLM$PhylumAcidobacteria_sqrt <- sqrt(tbl.GLM$PhylumAcidobacteria)
hist(tbl.GLM$PhylumAcidobacteria_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumAcidobacteria_sqrt <- glm(PhylumAcidobacteria_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                               data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumAcidobacteria_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumAcidobacteria_sqrt)
# W = 0.97607, p-value = 1.4e-05
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumAcidobacteria_sqrt), residuals(glm_transformedPhylumAcidobacteria_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumAcidobacteria_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumAcidobacteria_sqrt
# BP = 98.18, df = 63, p-value = 0.003019

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianPhylumAcidobacteria)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 1518.818  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedPhylumAcidobacteria_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -541.6198  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedPhylumAcidobacteria_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -2.589632  
### Choose log

###### Data transformation (log transformation)
tbl.GLM$PhylumAcidobacteria_log <- log(tbl.GLM$PhylumAcidobacteria)
glm_transformedPhylumAcidobacteria <- glm(PhylumAcidobacteria_log ~ Inoculants * Aggregate * Duration * Year, 
                                          data = tbl.GLM, family = gaussian)
Anova(glm_transformedPhylumAcidobacteria, type = "II") 

glm_transformedPhylumAcidobacteria_simplified <- glm(PhylumAcidobacteria_log ~ Inoculants + Aggregate + Duration + Year +
                                                       Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                       Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedPhylumAcidobacteria_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: PhylumAcidobacteria_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              0.425  1   0.514233    
# Aggregate              16.625  3   0.000844 ***
# Duration              128.371  3  < 2.2e-16 ***
# Year                    0.019  1   0.889583    
# Inoculants:Aggregate    1.118  3   0.772666    
# Inoculants:Duration     1.657  3   0.646462    
# Aggregate:Duration     10.939  9   0.279905    
# Inoculants:Year         0.523  1   0.469724    
# Aggregate:Year          1.484  3   0.685879    
# Duration:Year          38.243  3  2.511e-08 ***   


#################################  PhylumProteobacteria_model  
table(tbl.GLM$PhylumProteobacteria)
### test mean 37.03949
mean(tbl.GLM$PhylumProteobacteria)
### test variance 24.46513
var(tbl.GLM$PhylumProteobacteria)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumProteobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumProteobacteria) ## 0.4920181
kurtosis(tbl.GLM$PhylumProteobacteria) ## 0.8544188
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPhylumProteobacteria <- glm(PhylumProteobacteria ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPhylumProteobacteria))
qqline(residuals(glm_gaussianPhylumProteobacteria), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPhylumProteobacteria))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPhylumProteobacteria)
# W = 0.97075, p-value = 1.551e-06

## Test Homoscedasticity
plot(fitted(glm_gaussianPhylumProteobacteria), residuals(glm_gaussianPhylumProteobacteria), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPhylumProteobacteria)
# studentized Breusch-Pagan test
# data:  glm_gaussianPhylumProteobacteria
# BP = 98.065, df = 63, p-value = 0.00309

#### log transformation
tbl.GLM$PhylumProteobacteria_log <- log(tbl.GLM$PhylumProteobacteria)
hist(tbl.GLM$PhylumProteobacteria_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumProteobacteria_log <- glm(PhylumProteobacteria_log ~ Inoculants * Aggregate * Duration * Year, 
                                               data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumProteobacteria_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumProteobacteria_log)
# W = 0.98665, p-value = 0.002488
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumProteobacteria_log), residuals(glm_transformedPhylumProteobacteria_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumProteobacteria_log)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumProteobacteria_log
# BP = 89.536, df = 63, p-value = 0.01566

#### sqrt transformation
tbl.GLM$PhylumProteobacteria_sqrt <- sqrt(tbl.GLM$PhylumProteobacteria)
hist(tbl.GLM$PhylumProteobacteria_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumProteobacteria_sqrt <- glm(PhylumProteobacteria_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumProteobacteria_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumProteobacteria_sqrt)
# W = 0.97984, p-value = 7.745e-05
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumProteobacteria_sqrt), residuals(glm_transformedPhylumProteobacteria_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumProteobacteria_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumProteobacteria_sqrt
# BP = 94.014, df = 63, p-value = 0.00686

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianPhylumProteobacteria)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 2010.666  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedPhylumProteobacteria_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -560.8524  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedPhylumProteobacteria_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: 234.2449
### Choose log

###### Data transformation (log transformation)
tbl.GLM$PhylumProteobacteria_log <- log(tbl.GLM$PhylumProteobacteria)
glm_transformedPhylumProteobacteria <- glm(PhylumProteobacteria_log ~ Inoculants * Aggregate * Duration * Year, 
                                           data = tbl.GLM, family = gaussian)
Anova(glm_transformedPhylumProteobacteria, type = "II") 

glm_transformedPhylumProteobacteria_simplified <- glm(PhylumProteobacteria_log ~ Inoculants + Aggregate + Duration + Year +
                                                        Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                        Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedPhylumProteobacteria_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: PhylumProteobacteria_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             13.934  1  0.0001893 ***
# Aggregate               9.747  3  0.0208417 *  
# Duration               88.746  3  < 2.2e-16 ***
# Year                   10.983  1  0.0009193 ***
# Inoculants:Aggregate    2.216  3  0.5288720    
# Inoculants:Duration    10.608  3  0.0140439 *  
# Aggregate:Duration     30.416  9  0.0003726 ***
# Inoculants:Year        12.340  1  0.0004433 ***
# Aggregate:Year          0.935  3  0.8170153    
# Duration:Year          52.664  3  2.162e-11 ***


################################   PhylumActinobacteria_model  
table(tbl.GLM$PhylumActinobacteria)
### test mean 8.014687
mean(tbl.GLM$PhylumActinobacteria)
### test variance 3.99493
var(tbl.GLM$PhylumActinobacteria)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumActinobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumActinobacteria) ## 0.7607177
kurtosis(tbl.GLM$PhylumActinobacteria) ## 0.8690163
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPhylumActinobacteria <- glm(PhylumActinobacteria ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPhylumActinobacteria))
qqline(residuals(glm_gaussianPhylumActinobacteria), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPhylumActinobacteria))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPhylumActinobacteria)
# W = 0.98151, p-value = 0.0001722

## Test Homoscedasticity
plot(fitted(glm_gaussianPhylumActinobacteria), residuals(glm_gaussianPhylumActinobacteria), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPhylumActinobacteria)
# studentized Breusch-Pagan test
# data:  glm_gaussianPhylumActinobacteria
# BP = 112.12, df = 63, p-value = 0.0001391

#### log transformation
tbl.GLM$PhylumActinobacteria_log <- log(tbl.GLM$PhylumActinobacteria)
hist(tbl.GLM$PhylumActinobacteria_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumActinobacteria_log <- glm(PhylumActinobacteria_log ~ Inoculants * Aggregate * Duration * Year, 
                                               data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumActinobacteria_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumActinobacteria_log)
# W = 0.99378, p-value = 0.1584
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumActinobacteria_log), residuals(glm_transformedPhylumActinobacteria_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumActinobacteria_log)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumActinobacteria_log
# BP = 95.962, df = 63, p-value = 0.004702

#### sqrt transformation
tbl.GLM$PhylumActinobacteria_sqrt <- sqrt(tbl.GLM$PhylumActinobacteria)
hist(tbl.GLM$PhylumActinobacteria_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumActinobacteria_sqrt <- glm(PhylumActinobacteria_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumActinobacteria_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumActinobacteria_sqrt)
# W = 0.99095, p-value = 0.02937
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumActinobacteria_sqrt), residuals(glm_transformedPhylumActinobacteria_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumActinobacteria_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumActinobacteria_sqrt
# BP = 101.16, df = 63, p-value = 0.001628

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianPhylumActinobacteria)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 1148.523  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedPhylumActinobacteria_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -340.7782  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedPhylumActinobacteria_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -94.27729  
### Choose log

###### Data transformation (log transformation)
tbl.GLM$PhylumActinobacteria_log <- log(tbl.GLM$PhylumActinobacteria)
glm_transformedPhylumActinobacteria <- glm(PhylumActinobacteria_log ~ Inoculants * Aggregate * Duration * Year, 
                                           data = tbl.GLM, family = gaussian)
Anova(glm_transformedPhylumActinobacteria, type = "II") 

glm_transformedPhylumActinobacteria_simplified <- glm(PhylumActinobacteria_log ~ Inoculants + Aggregate + Duration + Year +
                                                        Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                        Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedPhylumActinobacteria_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: PhylumActinobacteria_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               9.34  1  0.0022434 ** 
# Aggregate              111.37  3  < 2.2e-16 ***
# Duration               408.91  3  < 2.2e-16 ***
# Year                    89.61  1  < 2.2e-16 ***
# Inoculants:Aggregate     2.99  3  0.3924395    
# Inoculants:Duration     36.67  3  5.392e-08 ***
# Aggregate:Duration      64.00  9  2.262e-10 ***
# Inoculants:Year          6.23  1  0.0125626 *  
# Aggregate:Year           1.58  3  0.6636143    
# Duration:Year           17.59  3  0.0005337 ***

################################   PhylumPlanctomycetes_model  
table(tbl.GLM$PhylumPlanctomycetes)
### test mean 6.04679
mean(tbl.GLM$PhylumPlanctomycetes)
### test variance 5.665538
var(tbl.GLM$PhylumPlanctomycetes)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumPlanctomycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumPlanctomycetes) ## 0.239093
kurtosis(tbl.GLM$PhylumPlanctomycetes) ## -0.7269954
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPhylumPlanctomycetes <- glm(PhylumPlanctomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPhylumPlanctomycetes))
qqline(residuals(glm_gaussianPhylumPlanctomycetes), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPhylumPlanctomycetes))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPhylumPlanctomycetes)
# W = 0.98733, p-value = 0.003631

## Test Homoscedasticity
plot(fitted(glm_gaussianPhylumPlanctomycetes), residuals(glm_gaussianPhylumPlanctomycetes), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPhylumPlanctomycetes)
# studentized Breusch-Pagan test
# data:  glm_gaussianPhylumPlanctomycetes
# BP = 91.569, df = 63, p-value = 0.01085

#### log transformation
tbl.GLM$PhylumPlanctomycetes_log <- log(tbl.GLM$PhylumPlanctomycetes)
hist(tbl.GLM$PhylumPlanctomycetes_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumPlanctomycetes_log <- glm(PhylumPlanctomycetes_log ~ Inoculants * Aggregate * Duration * Year, 
                                               data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumPlanctomycetes_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumPlanctomycetes_log)
# W = 0.9607, p-value = 4.11e-08
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumPlanctomycetes_log), residuals(glm_transformedPhylumPlanctomycetes_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumPlanctomycetes_log)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumPlanctomycetes_log
# BP = 94.566, df = 63, p-value = 0.006171

#### sqrt transformation
tbl.GLM$PhylumPlanctomycetes_sqrt <- sqrt(tbl.GLM$PhylumPlanctomycetes)
hist(tbl.GLM$PhylumPlanctomycetes_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumPlanctomycetes_sqrt <- glm(PhylumPlanctomycetes_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumPlanctomycetes_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumPlanctomycetes_sqrt)
# W = 0.97811, p-value = 3.475e-05
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumPlanctomycetes_sqrt), residuals(glm_transformedPhylumPlanctomycetes_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumPlanctomycetes_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumPlanctomycetes_sqrt
# BP = 82.828, df = 63, p-value = 0.04774

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianPhylumPlanctomycetes)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 1163.821  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedPhylumPlanctomycetes_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -10.01764  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedPhylumPlanctomycetes_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: 59.80279
### Choose log

###### Data transformation (log transformation)
tbl.GLM$PhylumPlanctomycetes_log <- log(tbl.GLM$PhylumPlanctomycetes)
glm_transformedPhylumPlanctomycetes <- glm(PhylumPlanctomycetes_log ~ Inoculants * Aggregate * Duration * Year, 
                                           data = tbl.GLM, family = gaussian)
Anova(glm_transformedPhylumPlanctomycetes, type = "II") 

glm_transformedPhylumPlanctomycetes_simplified <- glm(PhylumPlanctomycetes_log ~ Inoculants + Aggregate + Duration + Year +
                                                        Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                        Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedPhylumPlanctomycetes_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: PhylumPlanctomycetes_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               5.82  1    0.01583 *  
# Aggregate                2.61  3    0.45578    
# Duration               856.82  3    < 2e-16 ***
# Year                     1.47  1    0.22597    
# Inoculants:Aggregate     4.97  3    0.17378    
# Inoculants:Duration      5.29  3    0.15168    
# Aggregate:Duration       6.33  9    0.70697    
# Inoculants:Year          0.04  1    0.85104    
# Aggregate:Year           1.18  3    0.75831    
# Duration:Year            7.21  3    0.06562 . 


################################   PhylumFirmicutes_model  
table(tbl.GLM$PhylumFirmicutes)
### test mean 3.078466
mean(tbl.GLM$PhylumFirmicutes)
### test variance 0.4909179
var(tbl.GLM$PhylumFirmicutes)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumFirmicutes, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumFirmicutes) ## 0.8142599
kurtosis(tbl.GLM$PhylumFirmicutes) ## 0.3572686
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPhylumFirmicutes <- glm(PhylumFirmicutes ~ Inoculants * Aggregate * Duration * Year, 
                                    data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPhylumFirmicutes))
qqline(residuals(glm_gaussianPhylumFirmicutes), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPhylumFirmicutes))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPhylumFirmicutes)
# W = 0.96545, p-value = 2.126e-07

## Test Homoscedasticity
plot(fitted(glm_gaussianPhylumFirmicutes), residuals(glm_gaussianPhylumFirmicutes), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPhylumFirmicutes)
# studentized Breusch-Pagan test
# data:  glm_gaussianPhylumFirmicutes
# BP = 114.47, df = 63, p-value = 7.92e-05

#### log transformation
tbl.GLM$PhylumFirmicutes_log <- log(tbl.GLM$PhylumFirmicutes)
hist(tbl.GLM$PhylumFirmicutes_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumFirmicutes_log <- glm(PhylumFirmicutes_log ~ Inoculants * Aggregate * Duration * Year, 
                                           data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumFirmicutes_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumFirmicutes_log)
# W = 0.98902, p-value = 0.009501
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumFirmicutes_log), residuals(glm_transformedPhylumFirmicutes_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumFirmicutes_log)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumFirmicutes_log
# BP = 103, df = 63, p-value = 0.0011

#### sqrt transformation
tbl.GLM$PhylumFirmicutes_sqrt <- sqrt(tbl.GLM$PhylumFirmicutes)
hist(tbl.GLM$PhylumFirmicutes_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumFirmicutes_sqrt <- glm(PhylumFirmicutes_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                            data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumFirmicutes_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumFirmicutes_sqrt)
# W = 0.98123, p-value = 0.0001508
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumFirmicutes_sqrt), residuals(glm_transformedPhylumFirmicutes_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumFirmicutes_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumFirmicutes_sqrt
# BP = 109.87, df = 63, p-value = 0.0002361

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianPhylumFirmicutes)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 406.7474  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedPhylumFirmicutes_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -458.8467  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedPhylumFirmicutes_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -521.7871  
### Choose log

###### Data transformation (sqrt transformation)
tbl.GLM$PhylumFirmicutes_sqrt <- log(tbl.GLM$PhylumFirmicutes)
glm_transformedPhylumFirmicutes <- glm(PhylumFirmicutes_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                       data = tbl.GLM, family = gaussian)
Anova(glm_transformedPhylumFirmicutes, type = "II") 

glm_transformedPhylumFirmicutes_simplified <- glm(PhylumFirmicutes_sqrt ~ Inoculants + Aggregate + Duration + Year +
                                                    Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                    Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedPhylumFirmicutes_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: PhylumFirmicutes_sqrt
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               0.16  1   0.692748    
# Aggregate                4.04  3   0.257454    
# Duration               754.20  3  < 2.2e-16 ***
# Year                    17.84  1  2.407e-05 ***
# Inoculants:Aggregate     1.10  3   0.776630    
# Inoculants:Duration      5.59  3   0.133572    
# Aggregate:Duration      25.92  9   0.002104 ** 
# Inoculants:Year          0.71  1   0.399389    
# Aggregate:Year           0.35  3   0.951005    
# Duration:Year           28.95  3  2.291e-06 ***


################################   PhylumChloroflexi_model  
table(tbl.GLM$PhylumChloroflexi)
### test mean 1.782287
mean(tbl.GLM$PhylumChloroflexi)
### test variance 0.5204205
var(tbl.GLM$PhylumChloroflexi)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumChloroflexi, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumChloroflexi) ## 2.432886
kurtosis(tbl.GLM$PhylumChloroflexi) ## 11.88186
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaPhylumChloroflexi <- glm(PhylumChloroflexi ~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(PhylumChloroflexi ~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumChloroflexi
#                      LR Chisq Df Pr(>Chisq)   
# Inoculants              2.324  1   0.127402    
# Aggregate              99.132  3  < 2.2e-16 ***
# Duration              212.052  3  < 2.2e-16 ***
# Year                    1.110  1   0.291983    
# Inoculants:Aggregate    2.857  3   0.414207    
# Inoculants:Duration    15.085  3   0.001746 ** 
# Aggregate:Duration      8.493  9   0.485316    
# Inoculants:Year         0.058  1   0.808893    
# Aggregate:Year          1.708  3   0.635212    
# Duration:Year           7.742  3   0.051656 .  

################################   Phylumcandidate_division_WPS_1_model  
table(tbl.GLM$Phylumcandidate_division_WPS_1)
### test mean 2.506389
mean(tbl.GLM$Phylumcandidate_division_WPS_1)
### test variance 0.5585793
var(tbl.GLM$Phylumcandidate_division_WPS_1)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$Phylumcandidate_division_WPS_1, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$Phylumcandidate_division_WPS_1) ## 0.348791
kurtosis(tbl.GLM$Phylumcandidate_division_WPS_1) ## 0.1570266
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPhylumcandidate_division_WPS_1 <- glm(Phylumcandidate_division_WPS_1 ~ Inoculants * Aggregate * Duration * Year, 
                                                  data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPhylumcandidate_division_WPS_1))
qqline(residuals(glm_gaussianPhylumcandidate_division_WPS_1), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPhylumcandidate_division_WPS_1))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPhylumcandidate_division_WPS_1)
# W = 0.99295, p-value = 0.09682

## Test Homoscedasticity
plot(fitted(glm_gaussianPhylumcandidate_division_WPS_1), residuals(glm_gaussianPhylumcandidate_division_WPS_1), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPhylumcandidate_division_WPS_1)
# studentized Breusch-Pagan test
# data:  glm_gaussianPhylumcandidate_division_WPS_1
# BP = 92.44, df = 63, p-value = 0.009232

#### log transformation
tbl.GLM$Phylumcandidate_division_WPS_1_log <- log(tbl.GLM$Phylumcandidate_division_WPS_1)
hist(tbl.GLM$Phylumcandidate_division_WPS_1_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumcandidate_division_WPS_1_log <- glm(Phylumcandidate_division_WPS_1_log ~ Inoculants * Aggregate * Duration * Year, 
                                                         data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumcandidate_division_WPS_1_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumcandidate_division_WPS_1_log)
# W = 0.96698, p-value = 3.704e-07
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumcandidate_division_WPS_1_log), residuals(glm_transformedPhylumcandidate_division_WPS_1_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumcandidate_division_WPS_1_log)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumcandidate_division_WPS_1_log
# BP = 82.809, df = 63, p-value = 0.04789

#### sqrt transformation
tbl.GLM$Phylumcandidate_division_WPS_1_sqrt <- sqrt(tbl.GLM$Phylumcandidate_division_WPS_1)
hist(tbl.GLM$Phylumcandidate_division_WPS_1_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumcandidate_division_WPS_1_sqrt <- glm(Phylumcandidate_division_WPS_1_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                          data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumcandidate_division_WPS_1_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumcandidate_division_WPS_1_sqrt)
# W = 0.98582, p-value = 0.001587
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumcandidate_division_WPS_1_sqrt), residuals(glm_transformedPhylumcandidate_division_WPS_1_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumcandidate_division_WPS_1_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumcandidate_division_WPS_1_sqrt
# 85.409, df = 63, p-value = 0.03166

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianPhylumcandidate_division_WPS_1)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 650.6452  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedPhylumcandidate_division_WPS_1_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 72.51221  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedPhylumcandidate_division_WPS_1_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -141.9845  
### Choose log

###### Data transformation (log transformation)
tbl.GLM$Phylumcandidate_division_WPS_1_sqrt <- log(tbl.GLM$Phylumcandidate_division_WPS_1)
glm_transformedPhylumcandidate_division_WPS_1 <- glm(Phylumcandidate_division_WPS_1_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                     data = tbl.GLM, family = gaussian)
Anova(glm_transformedPhylumcandidate_division_WPS_1, type = "II") 

glm_transformedPhylumcandidate_division_WPS_1_simplified <- glm(Phylumcandidate_division_WPS_1_sqrt ~ Inoculants + Aggregate + Duration + Year +
                                                                  Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                                  Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedPhylumcandidate_division_WPS_1_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: Phylumcandidate_division_WPS_1_sqrt
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             20.430  1  6.185e-06 ***
# Aggregate               7.409  3   0.059939 .  
# Duration              141.504  3  < 2.2e-16 ***
# Year                   67.608  1  < 2.2e-16 ***
# Inoculants:Aggregate    2.887  3   0.409402    
# Inoculants:Duration    12.569  3   0.005669 ** 
# Aggregate:Duration     10.731  9   0.294607    
# Inoculants:Year         0.000  1   0.984996    
# Aggregate:Year          2.066  3   0.558829    
# Duration:Year           1.078  3   0.782499 


################################   PhylumVerrucomicrobia_model  
table(tbl.GLM$PhylumVerrucomicrobia)
### test mean 3.320486
mean(tbl.GLM$PhylumVerrucomicrobia)
### test variance 3.662034
var(tbl.GLM$PhylumVerrucomicrobia)
## The Poisson distribution requires the mean and variance to be close. 
## Construct Poisson GLM model
glm_gaussianPhylumVerrucomicrobia <- glm(PhylumVerrucomicrobia ~ Inoculants * Aggregate * Duration * Year, 
                                         family = poisson(link = "log"), data = tbl.GLM)
Anova(glm_gaussianPhylumVerrucomicrobia, type = "II") 

glm_gaussianPhylumVerrucomicrobia_simplified <- glm(PhylumVerrucomicrobia ~ Inoculants + Aggregate + Duration + Year +
                                                      Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                      Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_gaussianPhylumVerrucomicrobia_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: PhylumVerrucomicrobia
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              13.28  1  0.0002685 ***
# Aggregate               52.94  3  1.892e-11 ***
# Duration               671.56  3  < 2.2e-16 ***
# Year                    21.84  1  2.958e-06 ***
# Inoculants:Aggregate     3.43  3  0.3295673    
# Inoculants:Duration     20.01  3  0.0001686 ***
# Aggregate:Duration      26.20  9  0.0018972 ** 
# Inoculants:Year         20.34  1  6.480e-06 ***
# Aggregate:Year           0.57  3  0.9031828    
# Duration:Year           19.70  3  0.0001960 ***


################################   PhylumArmatimonadetes_model  
table(tbl.GLM$PhylumArmatimonadetes)
### test mean 1.49279
mean(tbl.GLM$PhylumArmatimonadetes)
### test variance 0.08420095
var(tbl.GLM$PhylumArmatimonadetes)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumArmatimonadetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumArmatimonadetes) ## 0.1731245
kurtosis(tbl.GLM$PhylumArmatimonadetes) ## -0.1319247
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPhylumArmatimonadetes <- glm(PhylumArmatimonadetes ~ Inoculants * Aggregate * Duration * Year, 
                                         data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPhylumArmatimonadetes))
qqline(residuals(glm_gaussianPhylumArmatimonadetes), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPhylumArmatimonadetes))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPhylumArmatimonadetes)
# W = 0.98421, p-value = 0.0006732

## Test Homoscedasticity
plot(fitted(glm_gaussianPhylumArmatimonadetes), residuals(glm_gaussianPhylumArmatimonadetes), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPhylumArmatimonadetes)
# studentized Breusch-Pagan test
# data:  glm_gaussianPhylumArmatimonadetes
# BP = 117.85, df = 63, p-value = 3.452e-05

#### log transformation
tbl.GLM$PhylumArmatimonadetes_log <- log(tbl.GLM$PhylumArmatimonadetes)
hist(tbl.GLM$PhylumArmatimonadetes_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumArmatimonadetes_log <- glm(PhylumArmatimonadetes_log ~ Inoculants * Aggregate * Duration * Year, 
                                                data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumArmatimonadetes_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumArmatimonadetes_log)
# W = 0.97865, p-value = 4.445e-05
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumArmatimonadetes_log), residuals(glm_transformedPhylumArmatimonadetes_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumArmatimonadetes_log)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumArmatimonadetes_log
# BP = 129.98, df = 63, p-value = 1.46e-06

#### sqrt transformation
tbl.GLM$PhylumArmatimonadetes_sqrt <- sqrt(tbl.GLM$PhylumArmatimonadetes)
hist(tbl.GLM$PhylumArmatimonadetes_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumArmatimonadetes_sqrt <- glm(PhylumArmatimonadetes_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                 data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumArmatimonadetes_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumArmatimonadetes_sqrt)
# W = 0.98337, p-value = 0.0004365
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumArmatimonadetes_sqrt), residuals(glm_transformedPhylumArmatimonadetes_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumArmatimonadetes_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumArmatimonadetes_sqrt
# BP = 123.85, df = 63, p-value = 7.466e-06

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianPhylumArmatimonadetes)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 20.86521  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedPhylumArmatimonadetes_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -241.5604  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedPhylumArmatimonadetes_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -604.2457  
### Choose log

###### Data transformation (sqrt transformation)
tbl.GLM$PhylumArmatimonadetes_sqrt <- log(tbl.GLM$PhylumArmatimonadetes)
glm_transformedPhylumArmatimonadetes <- glm(PhylumArmatimonadetes_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                            data = tbl.GLM, family = gaussian)
Anova(glm_transformedPhylumArmatimonadetes, type = "II") 

glm_transformedPhylumArmatimonadetes_simplified <- glm(PhylumArmatimonadetes_sqrt ~ Inoculants + Aggregate + Duration + Year +
                                                         Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                         Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedPhylumArmatimonadetes_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: PhylumArmatimonadetes_sqrt
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              9.494  1   0.002061 ** 
# Aggregate               6.888  3   0.075546 .  
# Duration               93.848  3  < 2.2e-16 ***
# Year                   21.732  1  3.135e-06 ***
# Inoculants:Aggregate    1.036  3   0.792613    
# Inoculants:Duration     4.149  3   0.245844    
# Aggregate:Duration      9.903  9   0.358397    
# Inoculants:Year         2.243  1   0.134263    
# Aggregate:Year          0.743  3   0.862969    
# Duration:Year          29.222  3  2.011e-06 ***


################################   PhylumThaumarchaeota_model  
table(tbl.GLM$PhylumThaumarchaeota)
### test mean 1.017309
mean(tbl.GLM$PhylumThaumarchaeota)
### test variance 0.8213671
var(tbl.GLM$PhylumThaumarchaeota)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumThaumarchaeota, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumThaumarchaeota) ## 1.767573
kurtosis(tbl.GLM$PhylumThaumarchaeota) ## 4.956235
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaPhylumThaumarchaeota<- glm(PhylumThaumarchaeota~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(PhylumThaumarchaeota~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumThaumarchaeota
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              46.17  1  1.082e-11 ***
# Aggregate               31.72  3  5.995e-07 ***
# Duration               867.37  3  < 2.2e-16 ***
# Year                   124.21  1  < 2.2e-16 ***
# Inoculants:Aggregate     2.02  3  0.5684885    
# Inoculants:Duration     20.12  3  0.0001604 ***
# Aggregate:Duration      16.81  9  0.0518389 .  
# Inoculants:Year          5.88  1  0.0152716 *  
# Aggregate:Year           1.63  3  0.6528828    
# Duration:Year            4.95  3  0.1758349    


################################   PhylumBacteroidetes_model  
table(tbl.GLM$PhylumBacteroidetes)
### test mean 6.712892
mean(tbl.GLM$PhylumBacteroidetes)
### test variance 15.78237
var(tbl.GLM$PhylumBacteroidetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonPhylumBacteroidetes <- glm(PhylumBacteroidetes ~ Inoculants * Aggregate * Duration * Year, 
                                      data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.3365201
deviance(glm_poissonPhylumBacteroidetes) / df.residual(glm_poissonPhylumBacteroidetes)
### Distribution of response variables
hist(tbl.GLM$PhylumBacteroidetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$PhylumBacteroidetes, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiPhylumBacteroidetes <- glm(PhylumBacteroidetes ~ Inoculants * Aggregate * Duration * Year, 
                                    data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiPhylumBacteroidetes), residuals(glm_quasiPhylumBacteroidetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiPhylumBacteroidetes, type = "II")

glm_quasiPhylumBacteroidetes_simplified <- glm(PhylumBacteroidetes ~ Inoculants + Aggregate + Duration + Year +
                                                 Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                 Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiPhylumBacteroidetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumBacteroidetes
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              42.61  1  6.669e-11 ***
# Aggregate                8.49  3   0.036833 *  
# Duration              1020.13  3  < 2.2e-16 ***
# Year                    54.66  1  1.432e-13 ***
# Inoculants:Aggregate     2.17  3   0.537662    
# Inoculants:Duration     58.06  3  1.530e-12 ***
# Aggregate:Duration      25.28  9   0.002673 ** 
# Inoculants:Year          0.04  1   0.849933    
# Aggregate:Year           0.10  3   0.991811    
# Duration:Year          115.70  3  < 2.2e-16 ***



################################   Phylumcandidate_division_WPS_2_model  
table(tbl.GLM$Phylumcandidate_division_WPS_2)
### test mean 7.424034
mean(tbl.GLM$Phylumcandidate_division_WPS_2)
### test variance 3.86082
var(tbl.GLM$Phylumcandidate_division_WPS_2)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$Phylumcandidate_division_WPS_2, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$Phylumcandidate_division_WPS_2) ## 0.4682367
kurtosis(tbl.GLM$Phylumcandidate_division_WPS_2) ## 0.5031752
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianPhylumcandidate_division_WPS_2 <- glm(Phylumcandidate_division_WPS_2 ~ Inoculants * Aggregate * Duration * Year, 
                                                  data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianPhylumcandidate_division_WPS_2))
qqline(residuals(glm_gaussianPhylumcandidate_division_WPS_2), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianPhylumcandidate_division_WPS_2))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianPhylumcandidate_division_WPS_2)
# W = 0.96715, p-value = 3.941e-07

## Test Homoscedasticity
plot(fitted(glm_gaussianPhylumcandidate_division_WPS_2), residuals(glm_gaussianPhylumcandidate_division_WPS_2), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianPhylumcandidate_division_WPS_2)
# studentized Breusch-Pagan test
# data:  glm_gaussianPhylumcandidate_division_WPS_2
# BP = 110, df = 63, p-value = 0.000229

#### log transformation
tbl.GLM$Phylumcandidate_division_WPS_2_log <- log(tbl.GLM$Phylumcandidate_division_WPS_2)
hist(tbl.GLM$Phylumcandidate_division_WPS_2_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumcandidate_division_WPS_2_log <- glm(Phylumcandidate_division_WPS_2_log ~ Inoculants * Aggregate * Duration * Year, 
                                                         data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumcandidate_division_WPS_2_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumcandidate_division_WPS_2_log)
# W = 0.97775, p-value = 2.945e-05
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumcandidate_division_WPS_2_log), residuals(glm_transformedPhylumcandidate_division_WPS_2_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumcandidate_division_WPS_2_log)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumcandidate_division_WPS_2_log
# BP = 103.94, df = 63, p-value = 0.0008975

#### sqrt transformation
tbl.GLM$Phylumcandidate_division_WPS_2_sqrt <- sqrt(tbl.GLM$Phylumcandidate_division_WPS_2)
hist(tbl.GLM$Phylumcandidate_division_WPS_2_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedPhylumcandidate_division_WPS_2_sqrt <- glm(Phylumcandidate_division_WPS_2_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                          data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedPhylumcandidate_division_WPS_2_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedPhylumcandidate_division_WPS_2_sqrt)
# W = 0.97867, p-value = 4.476e-05
## Test Homoscedasticity
plot(fitted(glm_transformedPhylumcandidate_division_WPS_2_sqrt), residuals(glm_transformedPhylumcandidate_division_WPS_2_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedPhylumcandidate_division_WPS_2_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedPhylumcandidate_division_WPS_2_sqrt
# BP = 107.17, df = 63, p-value = 0.0004382

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianPhylumcandidate_division_WPS_2)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 1339.88  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedPhylumcandidate_division_WPS_2_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -61.78741  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedPhylumcandidate_division_WPS_2_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: 137.9651 
### Choose log

###### Data transformation (log transformation)
tbl.GLM$Phylumcandidate_division_WPS_2_log <- log(tbl.GLM$Phylumcandidate_division_WPS_2)
glm_transformedPhylumcandidate_division_WPS_2 <- glm(Phylumcandidate_division_WPS_2_log ~ Inoculants * Aggregate * Duration * Year, 
                                                     data = tbl.GLM, family = gaussian)
Anova(glm_transformedPhylumcandidate_division_WPS_2, type = "II") 

glm_transformedPhylumcandidate_division_WPS_2_simplified <- glm(Phylumcandidate_division_WPS_2_log ~ Inoculants + Aggregate + Duration + Year +
                                                                  Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                                  Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedPhylumcandidate_division_WPS_2_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: Phylumcandidate_division_WPS_2_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             23.630  1  1.168e-06 ***
# Aggregate              50.837  3  5.300e-11 ***
# Duration              104.804  3  < 2.2e-16 ***
# Year                    5.015  1   0.025124 *  
# Inoculants:Aggregate    1.288  3   0.732000    
# Inoculants:Duration     8.604  3   0.035055 *  
# Aggregate:Duration      7.998  9   0.534395    
# Inoculants:Year        10.152  1   0.001442 ** 
# Aggregate:Year          1.707  3   0.635361    
# Duration:Year          21.530  3  8.170e-05 ***


################################   ClassAlphaproteobacteria_model  
table(tbl.GLM$ClassAlphaproteobacteria)
### test mean 15.06068
mean(tbl.GLM$ClassAlphaproteobacteria)
### test variance 21.04493
var(tbl.GLM$ClassAlphaproteobacteria)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassAlphaproteobacteria <- glm(ClassAlphaproteobacteria ~ Inoculants * Aggregate * Duration * Year, 
                                           data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.2309756
deviance(glm_poissonClassAlphaproteobacteria) / df.residual(glm_poissonClassAlphaproteobacteria)
### Distribution of response variables
hist(tbl.GLM$ClassAlphaproteobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassAlphaproteobacteria, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassAlphaproteobacteria <- glm(ClassAlphaproteobacteria ~ Inoculants * Aggregate * Duration * Year, 
                                         data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassAlphaproteobacteria), residuals(glm_quasiClassAlphaproteobacteria), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassAlphaproteobacteria, type = "II")

glm_quasiClassAlphaproteobacteria_simplified <- glm(ClassAlphaproteobacteria ~ Inoculants + Aggregate + Duration + Year +
                                                      Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                      Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassAlphaproteobacteria_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassAlphaproteobacteria
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               9.96  1   0.001603 ** 
# Aggregate               28.16  3  3.361e-06 ***
# Duration              1426.08  3  < 2.2e-16 ***
# Year                     1.06  1   0.302875    
# Inoculants:Aggregate     2.88  3   0.410044    
# Inoculants:Duration      9.25  3   0.026094 *  
# Aggregate:Duration      21.21  9   0.011735 *  
# Inoculants:Year          0.80  1   0.370276    
# Aggregate:Year           2.00  3   0.572085    
# Duration:Year           21.98  3  6.575e-05 ***


################################   ClassAcidobacteria_Gp1_model  
table(tbl.GLM$ClassAcidobacteria_Gp1)
### test mean 9.466193
mean(tbl.GLM$ClassAcidobacteria_Gp1)
### test variance 7.970478
var(tbl.GLM$ClassAcidobacteria_Gp1)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassAcidobacteria_Gp1, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassAcidobacteria_Gp1) ## 0.9532282
kurtosis(tbl.GLM$ClassAcidobacteria_Gp1) ## 0.2190687
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianClassAcidobacteria_Gp1 <- glm(ClassAcidobacteria_Gp1 ~ Inoculants * Aggregate * Duration * Year, 
                                          data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianClassAcidobacteria_Gp1))
qqline(residuals(glm_gaussianClassAcidobacteria_Gp1), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianClassAcidobacteria_Gp1))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianClassAcidobacteria_Gp1)
# W = 0.96289, p-value = 8.621e-08

## Test Homoscedasticity
plot(fitted(glm_gaussianClassAcidobacteria_Gp1), residuals(glm_gaussianClassAcidobacteria_Gp1), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianClassAcidobacteria_Gp1)
# studentized Breusch-Pagan test
# data:  glm_gaussianClassAcidobacteria_Gp1
# BP = 96.131, df = 63, p-value = 0.004548

#### log transformation
tbl.GLM$ClassAcidobacteria_Gp1_log <- log(tbl.GLM$ClassAcidobacteria_Gp1)
hist(tbl.GLM$ClassAcidobacteria_Gp1_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedClassAcidobacteria_Gp1_log <- glm(ClassAcidobacteria_Gp1_log ~ Inoculants * Aggregate * Duration * Year, 
                                                 data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedClassAcidobacteria_Gp1_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedClassAcidobacteria_Gp1_log)
# W = 0.98738, p-value = 0.003734
## Test Homoscedasticity
plot(fitted(glm_transformedClassAcidobacteria_Gp1_log), residuals(glm_transformedClassAcidobacteria_Gp1_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedClassAcidobacteria_Gp1_log)
# studentized Breusch-Pagan test
# data:  glm_transformedClassAcidobacteria_Gp1_log
# BP = 95.113, df = 63, p-value = 0.005551

#### sqrt transformation
tbl.GLM$ClassAcidobacteria_Gp1_sqrt <- sqrt(tbl.GLM$ClassAcidobacteria_Gp1)
hist(tbl.GLM$ClassAcidobacteria_Gp1_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedClassAcidobacteria_Gp1_sqrt <- glm(ClassAcidobacteria_Gp1_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                  data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedClassAcidobacteria_Gp1_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedClassAcidobacteria_Gp1_sqrt)
# W = 0.97881, p-value = 4.79e-05
## Test Homoscedasticity
plot(fitted(glm_transformedClassAcidobacteria_Gp1_sqrt), residuals(glm_transformedClassAcidobacteria_Gp1_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedClassAcidobacteria_Gp1_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedClassAcidobacteria_Gp1_sqrt
# BP = 93.295, df = 63, p-value = 0.007864

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianClassAcidobacteria_Gp1)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 1351.04  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedClassAcidobacteria_Gp1_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -268.2827  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedClassAcidobacteria_Gp1_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: 41.02213 
### Choose log

###### Data transformation (log transformation)
tbl.GLM$ClassAcidobacteria_Gp1_log <- log(tbl.GLM$ClassAcidobacteria_Gp1)
glm_transformedClassAcidobacteria_Gp1 <- glm(ClassAcidobacteria_Gp1_log ~ Inoculants * Aggregate * Duration * Year, 
                                             data = tbl.GLM, family = gaussian)
Anova(glm_transformedClassAcidobacteria_Gp1, type = "II") 

glm_transformedClassAcidobacteria_Gp1_simplified <- glm(ClassAcidobacteria_Gp1_log ~ Inoculants + Aggregate + Duration + Year +
                                                          Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                          Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedClassAcidobacteria_Gp1_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ClassAcidobacteria_Gp1_log
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               1.30  1    0.25338    
# Aggregate               24.07  3  2.415e-05 ***
# Duration               709.49  3  < 2.2e-16 ***
# Year                     0.06  1    0.81448    
# Inoculants:Aggregate     0.62  3    0.89122    
# Inoculants:Duration      2.08  3    0.55680    
# Aggregate:Duration      11.42  9    0.24832    
# Inoculants:Year          6.01  1    0.01422 *  
# Aggregate:Year           1.03  3    0.79371    
# Duration:Year           34.20  3  1.795e-07 ***


################################   ClassBetaproteobacteria_model  
table(tbl.GLM$ClassBetaproteobacteria)
### test mean 18.4917
mean(tbl.GLM$ClassBetaproteobacteria)
### test variance 48.59837
var(tbl.GLM$ClassBetaproteobacteria)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassBetaproteobacteria <- glm(ClassBetaproteobacteria ~ Inoculants * Aggregate * Duration * Year, 
                                          data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.3423959
deviance(glm_poissonClassBetaproteobacteria) / df.residual(glm_poissonClassBetaproteobacteria)
### Distribution of response variables
hist(tbl.GLM$ClassBetaproteobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassBetaproteobacteria, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassBetaproteobacteria <- glm(ClassBetaproteobacteria ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassBetaproteobacteria), residuals(glm_quasiClassBetaproteobacteria), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassBetaproteobacteria, type = "II")

glm_quasiClassBetaproteobacteria_simplified <- glm(ClassBetaproteobacteria ~ Inoculants + Aggregate + Duration + Year +
                                                     Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                     Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassBetaproteobacteria_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassBetaproteobacteria
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               7.10  1   0.007715 ** 
# Aggregate                2.84  3   0.416999    
# Duration              1859.16  3  < 2.2e-16 ***
# Year                     6.10  1   0.013535 *  
# Inoculants:Aggregate     2.28  3   0.517053    
# Inoculants:Duration      7.10  3   0.068705 .  
# Aggregate:Duration      25.94  9   0.002093 ** 
# Inoculants:Year         22.25  1  2.398e-06 ***
# Aggregate:Year           1.83  3   0.608302    
# Duration:Year           50.69  3  5.703e-11 ***


################################   ClassAcidobacteria_Gp3_model  
table(tbl.GLM$ClassAcidobacteria_Gp3)
### test mean 3.828153
mean(tbl.GLM$ClassAcidobacteria_Gp3)
### test variance 0.9875273
var(tbl.GLM$ClassAcidobacteria_Gp3)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassAcidobacteria_Gp3, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassAcidobacteria_Gp3) ## 1.336569
kurtosis(tbl.GLM$ClassAcidobacteria_Gp3) ## 3.23056
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaClassAcidobacteria_Gp3<- glm(ClassAcidobacteria_Gp3~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(ClassAcidobacteria_Gp3~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassAcidobacteria_Gp3
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             17.729  1  2.548e-05 ***
# Aggregate              31.990  3  5.259e-07 ***
# Duration              311.948  3  < 2.2e-16 ***
# Year                   45.297  1  1.693e-11 ***
# Inoculants:Aggregate    0.325  3    0.95529    
# Inoculants:Duration    23.704  3  2.880e-05 ***
# Aggregate:Duration      9.165  9    0.42216    
# Inoculants:Year        10.042  1    0.00153 ** 
# Aggregate:Year          2.132  3    0.54545    
# Duration:Year          44.347  3  1.274e-09 ***


################################   ClassPlanctomycetacia_model  
table(tbl.GLM$ClassPlanctomycetacia)
### test mean 6.045824
mean(tbl.GLM$ClassPlanctomycetacia)
### test variance 5.664682
var(tbl.GLM$ClassPlanctomycetacia)
## The Poisson distribution requires the mean and variance to be close. 
glm_gammaClassPlanctomycetacia<- glm(ClassPlanctomycetacia~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = poisson(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(ClassPlanctomycetacia~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = poisson(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassPlanctomycetacia
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              1.234  1     0.2666    
# Aggregate               0.644  3     0.8864    
# Duration              244.468  3     <2e-16 ***
# Year                    0.252  1     0.6159    
# Inoculants:Aggregate    0.797  3     0.8501    
# Inoculants:Duration     1.455  3     0.6928    
# Aggregate:Duration      1.517  9     0.9970    
# Inoculants:Year         0.408  1     0.5228    
# Aggregate:Year          0.344  3     0.9515    
# Duration:Year           1.903  3     0.5928  


################################   ClassActinobacteria_model  
table(tbl.GLM$ClassActinobacteria)
### test mean 4.879801
mean(tbl.GLM$ClassActinobacteria)
### test variance 1.8744
var(tbl.GLM$ClassActinobacteria)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassActinobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassActinobacteria) ## 1.146671
kurtosis(tbl.GLM$ClassActinobacteria) ## 1.802152
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaClassActinobacteria<- glm(ClassActinobacteria~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(ClassActinobacteria~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassActinobacteria
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              9.276  1   0.002322 ** 
# Aggregate             100.148  3  < 2.2e-16 ***
# Duration              298.271  3  < 2.2e-16 ***
# Year                   57.367  1  3.617e-14 ***
# Inoculants:Aggregate    3.201  3   0.361704    
# Inoculants:Duration    38.326  3  2.411e-08 ***
# Aggregate:Duration     51.729  9  5.086e-08 ***
# Inoculants:Year         5.090  1   0.024070 *  
# Aggregate:Year          2.052  3   0.561786    
# Duration:Year          22.751  3  4.550e-05 ***


################################   ClassGammaproteobacteria_model  
table(tbl.GLM$ClassGammaproteobacteria)
### test mean 2.069631
mean(tbl.GLM$ClassGammaproteobacteria)
### test variance 0.2802264
var(tbl.GLM$ClassGammaproteobacteria)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassGammaproteobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassGammaproteobacteria) ## 1.455876
kurtosis(tbl.GLM$ClassGammaproteobacteria) ## 5.163267
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaClassGammaproteobacteria<- glm(ClassGammaproteobacteria~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(ClassGammaproteobacteria~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassGammaproteobacteria
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              2.584  1  0.1079347    
# Aggregate              73.961  3  6.050e-16 ***
# Duration              229.320  3  < 2.2e-16 ***
# Year                   27.620  1  1.476e-07 ***
# Inoculants:Aggregate    4.232  3  0.2374654    
# Inoculants:Duration     4.372  3  0.2239866    
# Aggregate:Duration     22.478  9  0.0074805 ** 
# Inoculants:Year        13.769  1  0.0002067 ***
# Aggregate:Year          3.229  3  0.3576682    
# Duration:Year          28.040  3  3.562e-06 ***


################################   ClassBacilli_model  
table(tbl.GLM$ClassBacilli)
### test mean 1.666108
mean(tbl.GLM$ClassBacilli)
### test variance 0.2061332
var(tbl.GLM$ClassBacilli)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassBacilli, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassBacilli) ## 1.192681
kurtosis(tbl.GLM$ClassBacilli) ## 1.613083
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaClassBacilli<- glm(ClassBacilli~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(ClassBacilli~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassBacilli
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               0.11  1   0.741025    
# Aggregate                7.60  3   0.055098 .  
# Duration               488.09  3  < 2.2e-16 ***
# Year                     0.46  1   0.497103    
# Inoculants:Aggregate     0.76  3   0.858717    
# Inoculants:Duration      7.07  3   0.069572 .  
# Aggregate:Duration      27.15  9   0.001323 ** 
# Inoculants:Year          0.19  1   0.661728    
# Aggregate:Year           0.38  3   0.943630    
# Duration:Year           23.71  3  2.866e-05 ***


################################   ClassChitinophagia_model  
table(tbl.GLM$ClassChitinophagia)
### test mean 4.512347
mean(tbl.GLM$ClassChitinophagia)
### test variance 6.631325
var(tbl.GLM$ClassChitinophagia)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassChitinophagia <- glm(ClassChitinophagia ~ Inoculants * Aggregate * Duration * Year, 
                                     data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.2813815
deviance(glm_poissonClassChitinophagia) / df.residual(glm_poissonClassChitinophagia)
### Distribution of response variables
hist(tbl.GLM$ClassChitinophagia, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassChitinophagia, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassChitinophagia <- glm(ClassChitinophagia ~ Inoculants * Aggregate * Duration * Year, 
                                   data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassChitinophagia), residuals(glm_quasiClassChitinophagia), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassChitinophagia, type = "II")

glm_quasiClassChitinophagia_simplified <- glm(ClassChitinophagia ~ Inoculants + Aggregate + Duration + Year +
                                                Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassChitinophagia_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassChitinophagia
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              70.82  1  < 2.2e-16 ***
# Aggregate               46.01  3  5.644e-10 ***
# Duration               614.09  3  < 2.2e-16 ***
# Year                    61.13  1  5.348e-15 ***
# Inoculants:Aggregate     2.24  3  0.5234728    
# Inoculants:Duration     30.51  3  1.077e-06 ***
# Aggregate:Duration      37.24  9  2.380e-05 ***
# Inoculants:Year         12.52  1  0.0004018 ***
# Aggregate:Year           0.48  3  0.9223857    
# Duration:Year           97.89  3  < 2.2e-16 ***


################################   ClassSphingobacteriia_model  
table(tbl.GLM$ClassSphingobacteriia)
### test mean 2.138511
mean(tbl.GLM$ClassSphingobacteriia)
### test variance 2.973544
var(tbl.GLM$ClassSphingobacteriia)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassSphingobacteriia <- glm(ClassSphingobacteriia ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.1542039
deviance(glm_poissonClassSphingobacteriia) / df.residual(glm_poissonClassSphingobacteriia)
### Distribution of response variables
hist(tbl.GLM$ClassSphingobacteriia, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassSphingobacteriia, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassSphingobacteriia <- glm(ClassSphingobacteriia ~ Inoculants * Aggregate * Duration * Year, 
                                      data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassSphingobacteriia), residuals(glm_quasiClassSphingobacteriia), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassSphingobacteriia, type = "II")

glm_quasiClassSphingobacteriia_simplified <- glm(ClassSphingobacteriia ~ Inoculants + Aggregate + Duration + Year +
                                                   Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                   Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassSphingobacteriia_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassSphingobacteriia
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               1.78  1    0.18270    
# Aggregate               25.44  3  1.252e-05 ***
# Duration              1375.18  3  < 2.2e-16 ***
# Year                    21.28  1  3.972e-06 ***
# Inoculants:Aggregate     4.38  3    0.22366    
# Inoculants:Duration     98.05  3  < 2.2e-16 ***
# Aggregate:Duration      17.46  9    0.04196 *  
# Inoculants:Year         58.87  1  1.686e-14 ***
# Aggregate:Year           1.16  3    0.76240    
# Duration:Year          118.60  3  < 2.2e-16 ***


################################   ClassDeltaproteobacteria_model  
table(tbl.GLM$ClassDeltaproteobacteria)
### test mean 1.275537
mean(tbl.GLM$ClassDeltaproteobacteria)
### test variance 0.05073425
var(tbl.GLM$ClassDeltaproteobacteria)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassDeltaproteobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassDeltaproteobacteria) ## 0.04708381
kurtosis(tbl.GLM$ClassDeltaproteobacteria) ## 0.2403432
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianClassDeltaproteobacteria <- glm(ClassDeltaproteobacteria ~ Inoculants * Aggregate * Duration * Year, 
                                            data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianClassDeltaproteobacteria))
qqline(residuals(glm_gaussianClassDeltaproteobacteria), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianClassDeltaproteobacteria))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianClassDeltaproteobacteria)
# W = 0.98644, p-value = 0.002213

## Test Homoscedasticity
plot(fitted(glm_gaussianClassDeltaproteobacteria), residuals(glm_gaussianClassDeltaproteobacteria), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianClassDeltaproteobacteria)
# studentized Breusch-Pagan test
# data:  glm_gaussianClassDeltaproteobacteria
# BP = 76.179, df = 63, p-value = 0.1231

#### log transformation
tbl.GLM$ClassDeltaproteobacteria_log <- log(tbl.GLM$ClassDeltaproteobacteria)
hist(tbl.GLM$ClassDeltaproteobacteria_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedClassDeltaproteobacteria_log <- glm(ClassDeltaproteobacteria_log ~ Inoculants * Aggregate * Duration * Year, 
                                                   data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedClassDeltaproteobacteria_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedClassDeltaproteobacteria_log)
# W = 0.98749, p-value = 0.003964
## Test Homoscedasticity
plot(fitted(glm_transformedClassDeltaproteobacteria_log), residuals(glm_transformedClassDeltaproteobacteria_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedClassDeltaproteobacteria_log)
# studentized Breusch-Pagan test
# data:  glm_transformedClassDeltaproteobacteria_log
# BP = 93.997, df = 63, p-value = 0.006883

#### sqrt transformation
tbl.GLM$ClassDeltaproteobacteria_sqrt <- sqrt(tbl.GLM$ClassDeltaproteobacteria)
hist(tbl.GLM$ClassDeltaproteobacteria_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedClassDeltaproteobacteria_sqrt <- glm(ClassDeltaproteobacteria_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                    data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedClassDeltaproteobacteria_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedClassDeltaproteobacteria_sqrt)
# W = 0.9908, p-value = 0.02696
## Test Homoscedasticity
plot(fitted(glm_transformedClassDeltaproteobacteria_sqrt), residuals(glm_transformedClassDeltaproteobacteria_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedClassDeltaproteobacteria_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedClassDeltaproteobacteria_sqrt
# BP = 84.54, df = 63, p-value = 0.03645

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianClassDeltaproteobacteria)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -166.4643 
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedClassDeltaproteobacteria_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -305.81  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedClassDeltaproteobacteria_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -730.4451  
### Choose log

###### Data transformation (sqrt transformation)
tbl.GLM$ClassDeltaproteobacteria_sqrt <- log(tbl.GLM$ClassDeltaproteobacteria)
glm_transformedClassDeltaproteobacteria <- glm(ClassDeltaproteobacteria_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                               data = tbl.GLM, family = gaussian)
Anova(glm_transformedClassDeltaproteobacteria, type = "II") 

glm_transformedClassDeltaproteobacteria_simplified <- glm(ClassDeltaproteobacteria_sqrt ~ Inoculants + Aggregate + Duration + Year +
                                                            Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                            Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedClassDeltaproteobacteria_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ClassDeltaproteobacteria_sqrt
#                      LR Chisq Df Pr(>Chisq) 
# Inoculants              6.462  1  0.0110226 *  
# Aggregate              10.078  3  0.0179142 *  
# Duration              141.473  3  < 2.2e-16 ***
# Year                   32.434  1  1.233e-08 ***
# Inoculants:Aggregate    1.107  3  0.7753482    
# Inoculants:Duration    26.662  3  6.932e-06 ***
# Aggregate:Duration     28.347  9  0.0008346 ***
# Inoculants:Year         7.465  1  0.0062901 ** 
# Aggregate:Year          0.796  3  0.8505113    
# Duration:Year          23.009  3  4.020e-05 ***


################################   ClassKtedonobacteria_model  
table(tbl.GLM$ClassKtedonobacteria)
### test mean 1.398662
mean(tbl.GLM$ClassKtedonobacteria)
### test variance 0.4564266
var(tbl.GLM$ClassKtedonobacteria)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassKtedonobacteria, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassKtedonobacteria) ## 2.949611
kurtosis(tbl.GLM$ClassKtedonobacteria) ## 16.41688
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaClassKtedonobacteria<- glm(ClassKtedonobacteria~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(ClassKtedonobacteria~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassKtedonobacteria
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              4.006  1  0.0453347 *  
# Aggregate              68.547  3  8.737e-15 ***
# Duration              245.387  3  < 2.2e-16 ***
# Year                    6.723  1  0.0095199 ** 
# Inoculants:Aggregate    1.578  3  0.6644521    
# Inoculants:Duration    18.105  3  0.0004184 ***
# Aggregate:Duration      5.428  9  0.7954747    
# Inoculants:Year         1.156  1  0.2822123    
# Aggregate:Year          2.351  3  0.5027561    
# Duration:Year          12.632  3  0.0055028 ** 


################################   ClassThermoleophilia_model  
table(tbl.GLM$ClassThermoleophilia)
### test mean 2.495662
mean(tbl.GLM$ClassThermoleophilia)
### test variance 1.829129
var(tbl.GLM$ClassThermoleophilia)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassThermoleophilia, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassThermoleophilia) ## 0.8065601
kurtosis(tbl.GLM$ClassThermoleophilia) ## 0.876809
## The data is not right-skewed, 
#  so it is not appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gaussian
### Check the following assumptions
## Construct Gaussian GLM model
glm_gaussianClassThermoleophilia <- glm(ClassThermoleophilia ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = gaussian)
##  Q-Q Plot
qqnorm(residuals(glm_gaussianClassThermoleophilia))
qqline(residuals(glm_gaussianClassThermoleophilia), col = "red")

## Test Shapiro-Wilk
shapiro.test(residuals(glm_gaussianClassThermoleophilia))
# Shapiro-Wilk normality test
# data:  residuals(glm_gaussianClassThermoleophilia)
# W = 0.87095, p-value < 2.2e-16

## Test Homoscedasticity
plot(fitted(glm_gaussianClassThermoleophilia), residuals(glm_gaussianClassThermoleophilia), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_gaussianClassThermoleophilia)
# studentized Breusch-Pagan test
# data:  glm_gaussianClassThermoleophilia
# BP = 118.4, df = 63, p-value = 3.011e-05

#### log transformation
tbl.GLM$ClassThermoleophilia_log <- log(tbl.GLM$ClassThermoleophilia)
hist(tbl.GLM$ClassThermoleophilia_log, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedClassThermoleophilia_log <- glm(ClassThermoleophilia_log ~ Inoculants * Aggregate * Duration * Year, 
                                               data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedClassThermoleophilia_log))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedClassThermoleophilia_log)
# W = 0.98914, p-value = 0.01014
## Test Homoscedasticity
plot(fitted(glm_transformedClassThermoleophilia_log), residuals(glm_transformedClassThermoleophilia_log), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedClassThermoleophilia_log)
# studentized Breusch-Pagan test
# data:  glm_transformedClassThermoleophilia_log
# BP = 86.265, df = 63, p-value = 0.02748

#### sqrt transformation
tbl.GLM$ClassThermoleophilia_sqrt <- sqrt(tbl.GLM$ClassThermoleophilia)
hist(tbl.GLM$ClassThermoleophilia_sqrt, breaks = 20, main = "Histogram of Bacterial Richness")
glm_transformedClassThermoleophilia_sqrt <- glm(ClassThermoleophilia_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                                data = tbl.GLM, family = gaussian)
shapiro.test(residuals(glm_transformedClassThermoleophilia_sqrt))
# Shapiro-Wilk normality test
# data:  residuals(glm_transformedClassThermoleophilia_sqrt)
# W = 0.95839, p-value = 1.928e-08
## Test Homoscedasticity
plot(fitted(glm_transformedClassThermoleophilia_sqrt), residuals(glm_transformedClassThermoleophilia_sqrt), 
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
bptest(glm_transformedClassThermoleophilia_sqrt)
# studentized Breusch-Pagan test
# data:  glm_transformedClassThermoleophilia_sqrt
# BP = 111.11, df = 63, p-value = 0.0001768

# Calculate AIC for initial model
AIC_log <- AIC(glm_gaussianClassThermoleophilia)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: 508.272  
# Calculate AIC for log-transformed model
AIC_log <- AIC(glm_transformedClassThermoleophilia_log)
cat("AIC for log transformation:", AIC_log, "\n")
### AIC for log transformation: -264.0966  
# Calculate AIC for sqrt-transformed model
AIC_sqrt <- AIC(glm_transformedClassThermoleophilia_sqrt)
cat("AIC for sqrt transformation:", AIC_sqrt, "\n")
### AIC for sqrt transformation: -425.4098  
### Choose log

###### Data transformation (sqrt transformation)
tbl.GLM$ClassThermoleophilia_sqrt <- log(tbl.GLM$ClassThermoleophilia)
glm_transformedClassThermoleophilia <- glm(ClassThermoleophilia_sqrt ~ Inoculants * Aggregate * Duration * Year, 
                                           data = tbl.GLM, family = gaussian)
Anova(glm_transformedClassThermoleophilia, type = "II") 

glm_transformedClassThermoleophilia_simplified <- glm(ClassThermoleophilia_sqrt ~ Inoculants + Aggregate + Duration + Year +
                                                        Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                        Duration:Year, family = gaussian, data = tbl.GLM)
Anova(glm_transformedClassThermoleophilia_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: ClassThermoleophilia_sqrt
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants                3.9  1    0.04816 *  
# Aggregate                55.9  3  4.345e-12 ***
# Duration               4509.3  3  < 2.2e-16 ***
# Year                     36.0  1  1.936e-09 ***
# Inoculants:Aggregate      5.6  3    0.13004    
# Inoculants:Duration      27.1  3  5.552e-06 ***
# Aggregate:Duration       49.1  9  1.598e-07 ***
# Inoculants:Year           3.7  1    0.05378 .  
# Aggregate:Year            1.3  3    0.72753    
# Duration:Year            41.8  3  4.382e-09 ***


################################   ClassClostridia_model  
table(tbl.GLM$ClassClostridia)
### test mean 1.212455
mean(tbl.GLM$ClassClostridia)
### test variance 0.09535232
var(tbl.GLM$ClassClostridia)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$ClassClostridia, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$ClassClostridia) ## 0.906392
kurtosis(tbl.GLM$ClassClostridia) ## 1.640688
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
glm_gammaClassClostridia<- glm(ClassClostridia~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(ClassClostridia~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassClostridia
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               0.78  1   0.376775    
# Aggregate                4.54  3   0.208740    
# Duration               621.50  3  < 2.2e-16 ***
# Year                    94.38  1  < 2.2e-16 ***
# Inoculants:Aggregate     1.28  3   0.734951    
# Inoculants:Duration     58.45  3  1.259e-12 ***
# Aggregate:Duration      27.06  9   0.001368 ** 
# Inoculants:Year          5.90  1   0.015171 *  
# Aggregate:Year           0.82  3   0.844696    
# Duration:Year           27.51  3  4.596e-06 ***


################################   PhylumAscomycota_model  
table(tbl.GLM$PhylumAscomycota)
### test mean 66.9571
mean(tbl.GLM$PhylumAscomycota)
### test variance 94.90137
var(tbl.GLM$PhylumAscomycota)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonPhylumAscomycota <- glm(PhylumAscomycota ~ Inoculants * Aggregate * Duration * Year, 
                                   data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.387224
deviance(glm_poissonPhylumAscomycota) / df.residual(glm_poissonPhylumAscomycota)
### Distribution of response variables
hist(tbl.GLM$PhylumAscomycota, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$PhylumAscomycota, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiPhylumAscomycota <- glm(PhylumAscomycota ~ Inoculants * Aggregate * Duration * Year, 
                                 data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiPhylumAscomycota), residuals(glm_quasiPhylumAscomycota), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiPhylumAscomycota, type = "II")

glm_quasiPhylumAscomycota_simplified <- glm(PhylumAscomycota ~ Inoculants + Aggregate + Duration + Year +
                                              Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                              Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiPhylumAscomycota_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumAscomycota
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             1.1604  1    0.28139    
# Aggregate              0.9608  3    0.81074    
# Duration              23.8806  3  2.646e-05 ***
# Year                   4.7186  1    0.02984 *  
# Inoculants:Aggregate   2.7407  3    0.43336    
# Inoculants:Duration    1.8693  3    0.59997    
# Aggregate:Duration    12.8152  9    0.17115    
# Inoculants:Year        1.8755  1    0.17084    
# Aggregate:Year         3.7124  3    0.29424    
# Duration:Year          9.3054  3    0.02549 * 


################################   PhylumBasidiomycota_model  
table(tbl.GLM$PhylumBasidiomycota)
### test mean 6.392685
mean(tbl.GLM$PhylumBasidiomycota)
### test variance 11.60398
var(tbl.GLM$PhylumBasidiomycota)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonPhylumBasidiomycota <- glm(PhylumBasidiomycota ~ Inoculants * Aggregate * Duration * Year, 
                                      data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.467919
deviance(glm_poissonPhylumBasidiomycota) / df.residual(glm_poissonPhylumBasidiomycota)
### Distribution of response variables
hist(tbl.GLM$PhylumBasidiomycota, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$PhylumBasidiomycota, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiPhylumBasidiomycota <- glm(PhylumBasidiomycota ~ Inoculants * Aggregate * Duration * Year, 
                                    data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiPhylumBasidiomycota), residuals(glm_quasiPhylumBasidiomycota), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiPhylumBasidiomycota, type = "II")

glm_quasiPhylumBasidiomycota_simplified <- glm(PhylumBasidiomycota ~ Inoculants + Aggregate + Duration + Year +
                                                 Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                 Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiPhylumBasidiomycota_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumBasidiomycota
#                      LR Chisq Df Pr(>Chisq)   
# Inoculants             2.2557  1   0.133120   
# Aggregate              1.5222  3   0.677160   
# Duration              14.5686  3   0.002225 **
# Year                   5.9737  1   0.014521 * 
# Inoculants:Aggregate   0.4445  3   0.930895   
# Inoculants:Duration    2.7214  3   0.436594   
# Aggregate:Duration    15.7072  9   0.073253 . 
# Inoculants:Year        2.1728  1   0.140471   
# Aggregate:Year         2.6956  3   0.440980   
# Duration:Year          4.7505  3   0.191003   


################################   PhylumMortierellomycota_model  
table(tbl.GLM$PhylumMortierellomycota)
### test mean 5.744645
mean(tbl.GLM$PhylumMortierellomycota)
### test variance 8.775734
var(tbl.GLM$PhylumMortierellomycota)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonPhylumMortierellomycota <- glm(PhylumMortierellomycota ~ Inoculants * Aggregate * Duration * Year, 
                                          data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.288307
deviance(glm_poissonPhylumMortierellomycota) / df.residual(glm_poissonPhylumMortierellomycota)
### Distribution of response variables
hist(tbl.GLM$PhylumMortierellomycota, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$PhylumMortierellomycota, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiPhylumMortierellomycota <- glm(PhylumMortierellomycota ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiPhylumMortierellomycota), residuals(glm_quasiPhylumMortierellomycota), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiPhylumMortierellomycota, type = "II")

glm_quasiPhylumMortierellomycota_simplified <- glm(PhylumMortierellomycota ~ Inoculants + Aggregate + Duration + Year +
                                                     Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                     Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiPhylumMortierellomycota_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumMortierellomycota
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             1.3643  1  0.2427917    
# Aggregate             27.8971  3  3.817e-06 ***
# Duration              13.1481  3  0.0043271 ** 
# Year                  12.2616  1  0.0004624 ***
# Inoculants:Aggregate   0.7001  3  0.8731706    
# Inoculants:Duration    8.2373  3  0.0413534 *  
# Aggregate:Duration    11.3900  9  0.2499174    
# Inoculants:Year        6.1087  1  0.0134516 *  
# Aggregate:Year         2.1674  3  0.5384063    
# Duration:Year          5.8399  3  0.1196645  


#################################  PhylumMucoromycota_model  
table(tbl.GLM$PhylumMucoromycota)
### test mean 1.869344
mean(tbl.GLM$PhylumMucoromycota)
### test variance 3.368417
var(tbl.GLM$PhylumMucoromycota)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonPhylumMucoromycota <- glm(PhylumMucoromycota ~ Inoculants * Aggregate * Duration * Year, 
                                     data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.6843681
deviance(glm_poissonPhylumMucoromycota) / df.residual(glm_poissonPhylumMucoromycota)
### Distribution of response variables
hist(tbl.GLM$PhylumMucoromycota, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$PhylumMucoromycota, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiPhylumMucoromycota <- glm(PhylumMucoromycota ~ Inoculants * Aggregate * Duration * Year, 
                                   data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiPhylumMucoromycota), residuals(glm_quasiPhylumMucoromycota), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiPhylumMucoromycota, type = "II")

glm_quasiPhylumMucoromycota_simplified <- glm(PhylumMucoromycota ~ Inoculants + Aggregate + Duration + Year +
                                                Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiPhylumMucoromycota_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumMucoromycota
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              1.349  1    0.24542    
# Aggregate              41.688  3  4.673e-09 ***
# Duration               66.927  3  1.941e-14 ***
# Year                   42.465  1  7.196e-11 ***
# Inoculants:Aggregate    6.551  3    0.08766 .  
# Inoculants:Duration     4.495  3    0.21275    
# Aggregate:Duration      8.364  9    0.49790    
# Inoculants:Year         0.004  1    0.95253    
# Aggregate:Year          8.353  3    0.03925 *  
# Duration:Year           4.998  3    0.17193


################################   PhylumRozellomycota_model  
table(tbl.GLM$PhylumRozellomycota)
### test mean 0.0220142
mean(tbl.GLM$PhylumRozellomycota)
### test variance 0.001758635
var(tbl.GLM$PhylumRozellomycota)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumRozellomycota, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumRozellomycota) ## 4.67344
kurtosis(tbl.GLM$PhylumRozellomycota) ## 28.9093
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
tbl.GLM$PhylumRozellomycota <- tbl.GLM$PhylumRozellomycota + 0.01
glm_gammaPhylumRozellomycota<- glm(PhylumRozellomycota~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(PhylumRozellomycota~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumRozellomycota
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              3.194  1    0.07390 .  
# Aggregate               4.979  3    0.17332    
# Duration              182.219  3    < 2e-16 ***
# Year                    3.393  1    0.06549 .  
# Inoculants:Aggregate    1.677  3    0.64213    
# Inoculants:Duration     3.222  3    0.35865    
# Aggregate:Duration      9.117  9    0.42651    
# Inoculants:Year         1.000  1    0.31725    
# Aggregate:Year          5.233  3    0.15553    
# Duration:Year           8.335  3    0.03957 *


################################   PhylumChytridiomycota_model  
table(tbl.GLM$PhylumChytridiomycota)
### test mean 4.104122
mean(tbl.GLM$PhylumChytridiomycota)
### test variance 29.92413
var(tbl.GLM$PhylumChytridiomycota)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonPhylumChytridiomycota <- glm(PhylumChytridiomycota ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = poisson)
### The discrete factor was calculated  3.936918
deviance(glm_poissonPhylumChytridiomycota) / df.residual(glm_poissonPhylumChytridiomycota)
### Distribution of response variables
hist(tbl.GLM$PhylumChytridiomycota, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$PhylumChytridiomycota, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiPhylumChytridiomycota <- glm(PhylumChytridiomycota ~ Inoculants * Aggregate * Duration * Year, 
                                      data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiPhylumChytridiomycota), residuals(glm_quasiPhylumChytridiomycota), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiPhylumChytridiomycota, type = "II")

glm_quasiPhylumChytridiomycota_simplified <- glm(PhylumChytridiomycota ~ Inoculants + Aggregate + Duration + Year +
                                                   Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                   Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiPhylumChytridiomycota_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumChytridiomycota
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              0.280  1   0.596772    
# Aggregate               2.111  3   0.549679    
# Duration               58.006  3  1.567e-12 ***
# Year                    5.965  1   0.014589 *  
# Inoculants:Aggregate    1.101  3   0.776776    
# Inoculants:Duration     8.017  3   0.045655 *  
# Aggregate:Duration     19.081  9   0.024509 *  
# Inoculants:Year         0.001  1   0.977647    
# Aggregate:Year          2.672  3   0.445071    
# Duration:Year          12.248  3   0.006581 ** 


################################   PhylumGlomeromycota_model  
table(tbl.GLM$PhylumGlomeromycota)
### test mean 0.1016903
mean(tbl.GLM$PhylumGlomeromycota)
### test variance 0.02839399
var(tbl.GLM$PhylumGlomeromycota)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumGlomeromycota, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumGlomeromycota) ## 3.877086
kurtosis(tbl.GLM$PhylumGlomeromycota) ## 20.08614
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
tbl.GLM$PhylumGlomeromycota <- tbl.GLM$PhylumGlomeromycota + 0.01
glm_gammaPhylumGlomeromycota<- glm(PhylumGlomeromycota~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(PhylumGlomeromycota~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumGlomeromycota
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               7.85  1    0.00507 ** 
# Aggregate                2.80  3    0.42303    
# Duration               382.61  3  < 2.2e-16 ***
# Year                     0.12  1    0.73142    
# Inoculants:Aggregate     4.96  3    0.17449    
# Inoculants:Duration      3.61  3    0.30721    
# Aggregate:Duration       7.78  9    0.55615    
# Inoculants:Year          6.14  1    0.01324 *  
# Aggregate:Year           3.21  3    0.36016    
# Duration:Year           42.34  3  3.396e-09 ***


################################   PhylumKickxellomycota_model  
table(tbl.GLM$PhylumKickxellomycota)
### test mean 0.05242898
mean(tbl.GLM$PhylumKickxellomycota)
### test variance 0.02820294
var(tbl.GLM$PhylumKickxellomycota)
## The Poisson distribution requires the mean and variance to be close. 
### If the variance of the data is much smaller than the mean, the data is less discrete, 
### which may mean that the data is more concentrated or less variable, We can choose Gaussian or Gamma.
### Distribution of response variables
hist(tbl.GLM$PhylumKickxellomycota, breaks = 20, main = "Histogram of Bacterial Richness")
### Calculate the skewness and kurtosis of the data
skewness(tbl.GLM$PhylumKickxellomycota) ## 8.836338
kurtosis(tbl.GLM$PhylumKickxellomycota) ## 95.92426
## The data is right-skewed, 
#  so it is appropriate to use the Gamma distribution, 
#  because the Gamma distribution usually applies to right-skewed data
# Therefore, we choose Gamma
tbl.GLM$PhylumKickxellomycota <- tbl.GLM$PhylumKickxellomycota + 0.01
glm_gammaPhylumKickxellomycota<- glm(PhylumKickxellomycota~ Inoculants * Aggregate * Duration * Year, data = tbl.GLM, family = Gamma(link = "log"))
Anova(glm_gammaPhylumChloroflexi, type = "II")
glm_gammaPhylumChloroflexi_simplified <- glm(PhylumKickxellomycota~ Inoculants + Aggregate + Duration + Year +
                                               Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                               Duration:Year, family = Gamma(link = "log"), data = tbl.GLM)
Anova(glm_gammaPhylumChloroflexi_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: PhylumKickxellomycota
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             0.1628  1   0.686615    
# Aggregate              0.6989  3   0.873461    
# Duration              28.7953  3  2.473e-06 ***
# Year                   0.0038  1   0.950642    
# Inoculants:Aggregate   3.1802  3   0.364662    
# Inoculants:Duration   14.4659  3   0.002335 ** 
# Aggregate:Duration    11.7692  9   0.226634    
# Inoculants:Year        8.6032  1   0.003356 ** 
# Aggregate:Year         3.7404  3   0.290898    
# Duration:Year          6.0539  3   0.109015 


################################   ClassSordariomycetes_model  
table(tbl.GLM$ClassSordariomycetes)
### test mean 32.40472
mean(tbl.GLM$ClassSordariomycetes)
### test variance 160.3405
var(tbl.GLM$ClassSordariomycetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassSordariomycetes <- glm(ClassSordariomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                       data = tbl.GLM, family = poisson)
### The discrete factor was calculated  3.890873
deviance(glm_poissonClassSordariomycetes) / df.residual(glm_poissonClassSordariomycetes)
### Distribution of response variables
hist(tbl.GLM$ClassSordariomycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassSordariomycetes, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassSordariomycetes <- glm(ClassSordariomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                     data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassSordariomycetes), residuals(glm_quasiClassSordariomycetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassSordariomycetes, type = "II")

glm_quasiClassSordariomycetes_simplified <- glm(ClassSordariomycetes ~ Inoculants + Aggregate + Duration + Year +
                                                  Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                  Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassSordariomycetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassSordariomycetes
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             4.7276  1  0.0296820 *  
# Aggregate              3.1602  3  0.3675812    
# Duration              18.2112  3  0.0003979 ***
# Year                  10.5954  1  0.0011337 ** 
# Inoculants:Aggregate   1.9147  3  0.5903006    
# Inoculants:Duration    5.1181  3  0.1633530    
# Aggregate:Duration    23.5281  9  0.0051131 ** 
# Inoculants:Year        4.6369  1  0.0312907 *  
# Aggregate:Year         3.1682  3  0.3664111    
# Duration:Year         19.2102  3  0.0002474 ***


################################   ClassEurotiomycetes_model  
table(tbl.GLM$ClassEurotiomycetes)
### test mean 21.23724
mean(tbl.GLM$ClassEurotiomycetes)
### test variance 96.51153
var(tbl.GLM$ClassEurotiomycetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassEurotiomycetes <- glm(ClassEurotiomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                      data = tbl.GLM, family = poisson)
### The discrete factor was calculated  3.519167
deviance(glm_poissonClassEurotiomycetes) / df.residual(glm_poissonClassEurotiomycetes)
### Distribution of response variables
hist(tbl.GLM$ClassEurotiomycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassEurotiomycetes, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassEurotiomycetes <- glm(ClassEurotiomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                    data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassEurotiomycetes), residuals(glm_quasiClassEurotiomycetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassEurotiomycetes, type = "II")

glm_quasiClassEurotiomycetes_simplified <- glm(ClassEurotiomycetes ~ Inoculants + Aggregate + Duration + Year +
                                                 Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                 Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassEurotiomycetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassEurotiomycetes
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              6.071  1    0.01374 *  
# Aggregate              27.177  3  5.406e-06 ***
# Duration               43.179  3  2.255e-09 ***
# Year                    0.237  1    0.62630    
# Inoculants:Aggregate    4.768  3    0.18963    
# Inoculants:Duration     3.871  3    0.27572    
# Aggregate:Duration      3.934  9    0.91571    
# Inoculants:Year         0.000  1    0.99249    
# Aggregate:Year          1.404  3    0.70462    
# Duration:Year           6.563  3    0.08721 . 


################################   ClassAgaricomycetes_model  
table(tbl.GLM$ClassAgaricomycetes)
### test mean 2.482423
mean(tbl.GLM$ClassAgaricomycetes)
### test variance 4.137048
var(tbl.GLM$ClassAgaricomycetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassAgaricomycetes <- glm(ClassAgaricomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                      data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.9565022
deviance(glm_poissonClassAgaricomycetes) / df.residual(glm_poissonClassAgaricomycetes)
### Distribution of response variables
hist(tbl.GLM$ClassAgaricomycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassAgaricomycetes, "pois")
plot(fit)

##################### Poisson
glm_poissonClassAgaricomycetes <- glm(ClassAgaricomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                      data = tbl.GLM, family = poisson)

### Ensure model fit
plot(fitted(glm_poissonClassAgaricomycetes), residuals(glm_poissonClassAgaricomycetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_poissonClassAgaricomycetes, type = "II")

glm_poissonClassAgaricomycetes_simplified <- glm(ClassAgaricomycetes ~ Inoculants + Aggregate + Duration + Year +
                                                   Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                   Duration:Year, family = poisson, data = tbl.GLM)
Anova(glm_poissonClassAgaricomycetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassAgaricomycetes
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              8.042  1   0.004571 ** 
# Aggregate              10.719  3   0.013349 *  
# Duration               15.919  3   0.001178 ** 
# Year                    0.737  1   0.390746    
# Inoculants:Aggregate    3.084  3   0.378928    
# Inoculants:Duration     6.206  3   0.101999    
# Aggregate:Duration     35.239  9  5.406e-05 ***
# Inoculants:Year         0.017  1   0.897282    
# Aggregate:Year          4.695  3   0.195560    
# Duration:Year           4.683  3   0.196507  


################################   ClassDothideomycetes_model  
table(tbl.GLM$ClassDothideomycetes)
### test mean 4.581523
mean(tbl.GLM$ClassDothideomycetes)
### test variance 20.54686
var(tbl.GLM$ClassDothideomycetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassDothideomycetes <- glm(ClassDothideomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                       data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.803757
deviance(glm_poissonClassDothideomycetes) / df.residual(glm_poissonClassDothideomycetes)
### Distribution of response variables
hist(tbl.GLM$ClassDothideomycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassDothideomycetes, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassDothideomycetes <- glm(ClassDothideomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                     data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassDothideomycetes), residuals(glm_quasiClassDothideomycetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassDothideomycetes, type = "II")

glm_quasiClassDothideomycetes_simplified <- glm(ClassDothideomycetes ~ Inoculants + Aggregate + Duration + Year +
                                                  Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                  Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassDothideomycetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassDothideomycetes
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             5.4436  1   0.019640 *  
# Aggregate             14.1779  3   0.002673 ** 
# Duration              25.9020  3  9.999e-06 ***
# Year                   0.0770  1   0.781397    
# Inoculants:Aggregate   4.9948  3   0.172180    
# Inoculants:Duration    6.6645  3   0.083398 .  
# Aggregate:Duration    16.1362  9   0.064091 .  
# Inoculants:Year        0.0086  1   0.926139    
# Aggregate:Year         5.2818  3   0.152287    
# Duration:Year         13.0065  3   0.004623 **


################################   ClassMortierellomycetes_model  
table(tbl.GLM$ClassMortierellomycetes)
### test mean 5.744645
mean(tbl.GLM$ClassMortierellomycetes)
### test variance 8.775734
var(tbl.GLM$ClassMortierellomycetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassMortierellomycetes <- glm(ClassMortierellomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                          data = tbl.GLM, family = poisson)
### The discrete factor was calculated  1.288307
deviance(glm_poissonClassMortierellomycetes) / df.residual(glm_poissonClassMortierellomycetes)
### Distribution of response variables
hist(tbl.GLM$ClassMortierellomycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassMortierellomycetes, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassMortierellomycetes <- glm(ClassMortierellomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                        data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassMortierellomycetes), residuals(glm_quasiClassMortierellomycetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassMortierellomycetes, type = "II")

glm_quasiClassMortierellomycetes_simplified <- glm(ClassMortierellomycetes ~ Inoculants + Aggregate + Duration + Year +
                                                     Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                     Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassMortierellomycetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassMortierellomycetes
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants             1.3643  1  0.2427917    
# Aggregate             27.8971  3  3.817e-06 ***
# Duration              13.1481  3  0.0043271 ** 
# Year                  12.2616  1  0.0004624 ***
# Inoculants:Aggregate   0.7001  3  0.8731706    
# Inoculants:Duration    8.2373  3  0.0413534 *  
# Aggregate:Duration    11.3900  9  0.2499174    
# Inoculants:Year        6.1087  1  0.0134516 *  
# Aggregate:Year         2.1674  3  0.5384063    
# Duration:Year          5.8399  3  0.1196645 


#################################  ClassLecanoromycetes_model  
table(tbl.GLM$ClassLecanoromycetes)
### test mean 0.5273182
mean(tbl.GLM$ClassLecanoromycetes)
### test variance 5.70306
var(tbl.GLM$ClassLecanoromycetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassLecanoromycetes <- glm(ClassLecanoromycetes ~ Inoculants * Aggregate * Duration * Year, 
                                       data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.6228756
deviance(glm_poissonClassLecanoromycetes) / df.residual(glm_poissonClassLecanoromycetes)
### Distribution of response variables
hist(tbl.GLM$ClassLecanoromycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassLecanoromycetes, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassLecanoromycetes <- glm(ClassLecanoromycetes ~ Inoculants * Aggregate * Duration * Year, 
                                     data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassLecanoromycetes), residuals(glm_quasiClassLecanoromycetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassLecanoromycetes, type = "II")

glm_quasiClassLecanoromycetes_simplified <- glm(ClassLecanoromycetes ~ Inoculants + Aggregate + Duration + Year +
                                                  Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                  Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassLecanoromycetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassLecanoromycetes
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              5.300  1   0.021322 *  
# Aggregate              38.864  3  1.855e-08 ***
# Duration               57.947  3  1.614e-12 ***
# Year                    6.697  1   0.009655 ** 
# Inoculants:Aggregate    3.086  3   0.378559    
# Inoculants:Duration     6.415  3   0.093063 .  
# Aggregate:Duration     12.990  9   0.163058    
# Inoculants:Year        15.376  1  8.810e-05 ***
# Aggregate:Year          5.235  3   0.155386    
# Duration:Year           1.689  3   0.639355


#################################  ClassLeotiomycetes_model  
table(tbl.GLM$ClassLeotiomycetes)
### test mean 6.724097
mean(tbl.GLM$ClassLeotiomycetes)
### test variance 27.43008
var(tbl.GLM$ClassLeotiomycetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassLeotiomycetes <- glm(ClassLeotiomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                     data = tbl.GLM, family = poisson)
### The discrete factor was calculated  2.519842
deviance(glm_poissonClassLeotiomycetes) / df.residual(glm_poissonClassLeotiomycetes)
### Distribution of response variables
hist(tbl.GLM$ClassLeotiomycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassLeotiomycetes, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassLeotiomycetes <- glm(ClassLeotiomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                   data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassLeotiomycetes), residuals(glm_quasiClassLeotiomycetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassLeotiomycetes, type = "II")

glm_quasiClassLeotiomycetes_simplified <- glm(ClassLeotiomycetes ~ Inoculants + Aggregate + Duration + Year +
                                                Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassLeotiomycetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassLeotiomycetes
#                      LR Chisq Df Pr(>Chisq)  
# Inoculants             1.0189  1    0.31279  
# Aggregate              1.7212  3    0.63224  
# Duration               9.3279  3    0.02523 *
# Year                   3.5703  1    0.05882 .
# Inoculants:Aggregate   2.6084  3    0.45603  
# Inoculants:Duration    1.6502  3    0.64805  
# Aggregate:Duration     6.0553  9    0.73437  
# Inoculants:Year        0.7504  1    0.38634  
# Aggregate:Year         0.9576  3    0.81152  
# Duration:Year          7.3905  3    0.06044 .


#################################  ClassEndogonomycetes_model  
table(tbl.GLM$ClassEndogonomycetes)
### test mean 0.3244261
mean(tbl.GLM$ClassEndogonomycetes)
### test variance 0.4263041
var(tbl.GLM$ClassEndogonomycetes)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.
#################### Fitting the Poisson model
glm_poissonClassEndogonomycetes <- glm(ClassEndogonomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                       data = tbl.GLM, family = poisson)
### The discrete factor was calculated  0.1959901
deviance(glm_poissonClassEndogonomycetes) / df.residual(glm_poissonClassEndogonomycetes)
### Distribution of response variables
hist(tbl.GLM$ClassEndogonomycetes, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$ClassEndogonomycetes, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiClassEndogonomycetes <- glm(ClassEndogonomycetes ~ Inoculants * Aggregate * Duration * Year, 
                                     data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiClassEndogonomycetes), residuals(glm_quasiClassEndogonomycetes), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiClassEndogonomycetes, type = "II")

glm_quasiClassEndogonomycetes_simplified <- glm(ClassEndogonomycetes ~ Inoculants + Aggregate + Duration + Year +
                                                  Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                                  Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiClassEndogonomycetes_simplified, type = "II")
# Analysis of Deviance Table (Type II tests)
# Response: ClassEndogonomycetes
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants               1.20  1  0.2733499    
# Aggregate               52.81  3  2.011e-11 ***
# Duration               373.86  3  < 2.2e-16 ***
# Year                   128.50  1  < 2.2e-16 ***
# Inoculants:Aggregate    16.81  3  0.0007722 ***
# Inoculants:Duration      4.47  3  0.2154054    
# Aggregate:Duration      21.26  9  0.0115379 *  
# Inoculants:Year          9.16  1  0.0024727 ** 
# Aggregate:Year           2.85  3  0.4159681    
# Duration:Year           46.26  3  4.982e-10 ***


#################################  Starch_model  
table(tbl.GLM$Starch)
### test mean 17936.41
mean(tbl.GLM$Starch)
### test variance 13186562
var(tbl.GLM$Starch)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonStarch <- glm(Starch ~ Inoculants * Aggregate * Duration * Year, 
                         data = tbl.GLM, family = poisson)
### The discrete factor was calculated  71.31744
deviance(glm_poissonStarch) / df.residual(glm_poissonStarch)
### Distribution of response variables
hist(tbl.GLM$Starch, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$Starch, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiStarch <- glm(Starch ~ Inoculants * Aggregate * Duration * Year, 
                       data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiStarch), residuals(glm_quasiStarch), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiStarch, type = "II")

glm_quasiStarch_simplified <- glm(Starch ~ Inoculants + Aggregate + Duration + Year +
                                    Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                    Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiStarch_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: Starch
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              27.99  1  1.218e-07 ***
# Aggregate               19.09  3  0.0002616 ***
# Duration              1635.17  3  < 2.2e-16 ***
# Year                    46.78  1  7.953e-12 ***
# Inoculants:Aggregate     0.97  3  0.8096829    
# Inoculants:Duration     60.93  3  3.713e-13 ***
# Aggregate:Duration      14.61  9  0.1022329    
# Inoculants:Year          0.77  1  0.3806322    
# Aggregate:Year           1.18  3  0.7576342    
# Duration:Year           87.76  3  < 2.2e-16 ***


#################################  Hemicellulose_model  
table(tbl.GLM$Hemicellulose)
### test mean 111218.2
mean(tbl.GLM$Hemicellulose)
### test variance 856927842
var(tbl.GLM$Hemicellulose)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonHemicellulose <- glm(Hemicellulose ~ Inoculants * Aggregate * Duration * Year, 
                                data = tbl.GLM, family = poisson)
### The discrete factor was calculated  717.4378
deviance(glm_poissonHemicellulose) / df.residual(glm_poissonHemicellulose)
### Distribution of response variables
hist(tbl.GLM$Hemicellulose, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$Hemicellulose, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiHemicellulose <- glm(Hemicellulose ~ Inoculants * Aggregate * Duration * Year, 
                              data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiHemicellulose), residuals(glm_quasiHemicellulose), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiHemicellulose, type = "II")

glm_quasiHemicellulose_simplified <- glm(Hemicellulose ~ Inoculants + Aggregate + Duration + Year +
                                           Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                           Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiHemicellulose_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: Hemicellulose
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              17.00  1  3.730e-05 ***
# Aggregate                8.05  3   0.044918 *  
# Duration              1453.46  3  < 2.2e-16 ***
# Year                    26.02  1  3.377e-07 ***
# Inoculants:Aggregate     1.28  3   0.733662    
# Inoculants:Duration     65.99  3  3.080e-14 ***
# Aggregate:Duration      27.35  9   0.001223 ** 
# Inoculants:Year          0.70  1   0.401318    
# Aggregate:Year           1.26  3   0.739772    
# Duration:Year          117.05  3  < 2.2e-16 ***


#################################  Cellulose_model  
table(tbl.GLM$Cellulose)
### test mean 142130.8
mean(tbl.GLM$Cellulose)
### test variance 853949083
var(tbl.GLM$Cellulose)
## The Poisson distribution requires the mean and variance to be close. 
## If the variance is much larger than the mean, there may be an overspread.

#################### Fitting the Poisson model
glm_poissonCellulose <- glm(Cellulose ~ Inoculants * Aggregate * Duration * Year, 
                            data = tbl.GLM, family = poisson)
### The discrete factor was calculated  434.2489
deviance(glm_poissonCellulose) / df.residual(glm_poissonCellulose)
### Distribution of response variables
hist(tbl.GLM$Cellulose, breaks = 20, main = "Histogram of Bacterial Richness")
### Fitting the theoretical Poisson distribution
fit <- fitdist(tbl.GLM$Cellulose, "pois")
plot(fit)

##################### Quasi-Poisson
glm_quasiCellulose <- glm(Cellulose ~ Inoculants * Aggregate * Duration * Year, 
                          data = tbl.GLM, family = quasipoisson)

### Ensure model fit
plot(fitted(glm_quasiCellulose), residuals(glm_quasiCellulose), 
     main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, col = "red")
### Results
Anova(glm_quasiCellulose, type = "II")

glm_quasiCellulose_simplified <- glm(Cellulose ~ Inoculants + Aggregate + Duration + Year +
                                       Inoculants:Aggregate + Inoculants:Duration + Aggregate:Duration + Inoculants:Year + Aggregate:Year + 
                                       Duration:Year, family = quasipoisson, data = tbl.GLM)
Anova(glm_quasiCellulose_simplified, type = "II") 
# Analysis of Deviance Table (Type II tests)
# Response: Cellulose
#                      LR Chisq Df Pr(>Chisq)    
# Inoculants              13.45  1  0.0002448 ***
# Aggregate               18.36  3  0.0003703 ***
# Duration              1835.26  3  < 2.2e-16 ***
# Year                    21.83  1  2.975e-06 ***
# Inoculants:Aggregate     1.24  3  0.7433987    
# Inoculants:Duration     61.55  3  2.746e-13 ***
# Aggregate:Duration      14.84  9  0.0953947 .  
# Inoculants:Year          0.67  1  0.4143485    
# Aggregate:Year           0.61  3  0.8942380    
# Duration:Year           97.71  3  < 2.2e-16 ***


############################# Relationship between MAOC and microbial functions
library(ggplot2)
library(dplyr)
data <- read.csv("FunctionandSOCMWD.csv")
variables <- c(
  "Starch_A", "Hemicellulose_A", "Cellulose_A",
  "Starch_B", "Hemicellulose_B", "Cellulose_B",
  "Starch_C", "Hemicellulose_C", "Cellulose_C", 
  "Starch_D", "Hemicellulose_D", "Cellulose_D")
targets <- c("MAOC", "MWD")
custom_colors <- c("1CK" = "#F7AF34", "1TR" = "#448DCD", "3CK" = "#ffcd85", "3TR" = "#a1b8d5", "None" = "#BEBEBE")
custom_shapes <- c("10day" = 17, "30day" = 15, "60day" = 16, "100day" = 18)
plot_and_save <- function(x_var, y_var) {
  cor_test <- cor.test(data[[x_var]], data[[y_var]], method = "spearman")
  plot <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(aes(color = Color, shape = Shape), size = 8) + 
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", fill = "gray80") +
    scale_color_manual(values = custom_colors) + 
    scale_shape_manual(values = custom_shapes) + 
    annotate(
      "text", 
      x = quantile(data[[x_var]], 0.75, na.rm = TRUE),
      y = quantile(data[[y_var]], 0.95, na.rm = TRUE), 
      label = paste("r =", round(cor_test$estimate, 3), 
                    "\np =", format.pval(cor_test$p.value, digits = 3)), 
      size = 5, hjust = 0
    ) +
    labs(
      x = x_var,
      y = y_var
    ) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank(), 
      axis.line = element_line(color = "black"),
      legend.position = "right"
    )
  # Output
  pdf_name <- paste0("Scatterplot_", x_var, "_with_", y_var, ".pdf")
  pdf(pdf_name, width = 8, height = 8)
  print(plot)
  dev.off()
}
for (var in variables) {
  for (target in targets) {
    plot_and_save(var, target)
  }
}








