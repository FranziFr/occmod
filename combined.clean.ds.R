## R code for analysising data collected for occupancy modelling
## from Ny Friesland, Spitsbergen
## and Utah, USA
##
## clean up workspace
rm(list=ls())

# load required libraries
library(RPresence)
library(gplots)
library(sp)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggpubr)

## set working directory
setwd("//kant/nhm-sfs-u2/franzif/paper3")

## loading data from Utah
## Wah Wah and Juab formations
UT.data <- read.csv("occ_mod_fielddata_Utah.csv", sep = ";", header=T)
UT.data <- cbind.data.frame(UT.data, "locality"="UT")
UT.data <- select(UT.data, site, subsite, section, ft.Log, facies, Nautiloids, Trilobites, Gastropods, Brachiopod.art, Crinoid, Echinoderm, Cystoid, sponge, area, locality)
## renaming some cols, so that they fit the Ny Friesland dataset
names(UT.data)[names(UT.data) == "Brachiopod.art"] <- "Brachiopods"
names(UT.data)[names(UT.data) == "Echinoderm"] <- "Echinodermata"
UT.data$section[UT.data$section=="HHJ"]<-"HJ"
# UT.data.old <- UT.data

UT.data <- UT.data %>% mutate(Echinoderms = rowSums(UT.data[,c("Crinoid", "Echinodermata","Cystoid")],na.rm = TRUE))



## loading data from Ny Friesland
## Valhallfonna Fm. / Olenidsletta Mb.
NF.data <- read.csv("NF_data.csv", sep = ";", header=T)
NF.PO.data <- filter(NF.data, grepl('PO', m.Log))
NF.PO.data.a <- select(NF.PO.data, site, subsite, facies, Nautiloids, Trilobites, Gastropods, Graptolites, Brachiopods, Echinodermata, sponge.spicules., area)

## renaming some cols, so that they fit the Utah dataset
names(NF.PO.data.a)[names(NF.PO.data.a) == "sponge.spicules."] <- "sponge"
names(NF.PO.data.a)[names(NF.PO.data.a) == "Echinodermata"] <- "Echinoderms"

NF.PO.data.m <- colsplit(NF.PO.data$m.Log, "_", names = c("Log","m"))[,2]
NF.PO.data <- cbind.data.frame(NF.PO.data.a, "m"=NF.PO.data.m, "locality"="NF")

## create a section input for lower and upper section in Olenidsletta Mb.
NF.PO.data <- NF.PO.data %>%
  mutate(section = ifelse(m < 70, "Low", "Up"))

## as UT.data has data from 28 sites, and I don't want a doubling of site numbers, I adapt the NF.dataset sitenumbers by increasing them by 28, so that site numbers are running
NF.PO.data$site <- NF.PO.data$site+28 

## combined dataframe from Ny Friesland and Utah
##
df <- bind_rows(UT.data,NF.PO.data)
df.taxa <- select(df, site, subsite, Nautiloids, Trilobites, Brachiopods, Gastropods, Echinoderms)
df.taxa[is.na(df.taxa)] <- 0


## create taxa-specific detection/non.detection matrices
Nautiloids <- dcast(df.taxa, site~subsite, value.var = "Nautiloids")[,-1]
Trilobites <- dcast(df.taxa, site~subsite, value.var = "Trilobites")[,-1]
Brachiopods <- dcast(df.taxa, site~subsite, value.var = "Brachiopods")[,-1]
Gastropods <- dcast(df.taxa, site~subsite, value.var = "Gastropods")[,-1] 
# Cystoid <- dcast(df.taxa, site~subsite, value.var = "Cystoid")[,-1]
Echinoderms <- dcast(df.taxa, site~subsite, value.var = "Echinoderms")[,-1]

## the above data aboundance counts
## as only detection/non.detection is of interest in the present analysis, all values >=1 are changed to 1
Nautiloids <- ifelse(Nautiloids >=1, 1, 0)
Trilobites <- ifelse(Trilobites >=1, 1, 0)
Brachiopods <- ifelse(Brachiopods >=1, 1, 0)
Gastropods <- ifelse(Gastropods >=1, 1, 0)
# Cystoid <- ifelse(Cystoid >=1, 1, 0)
Echinoderms <- ifelse(Echinoderms >=1, 1, 0)

## creating input covariate file
covariate <- unique(df[,c("site","section", "facies", "locality", "m", "area", "ft.Log")])
covariate <- covariate[-1]

# covariate <- covariate %>% mutate(area.sq=case_when(
#   area %in% "1*1" ~ 1,
#   area %in% "1*2" ~ 2,
#   area %in% "1*3" ~ 3,
#   area %in% "2*2" ~ 4,
#   area %in% "2*3" ~ 6,
#   area %in% "3*3" ~ 9
# ))
covariate <- covariate %>% mutate(area.sq=case_when(
  area %in% "1*1" ~ 625,
  area %in% "1*2" ~ 1250,
  area %in% "1*3" ~ 1875,
  area %in% "2*2" ~ 2500,
  area %in% "2*3" ~ 3750,
  area %in% "3*3" ~ 5625
))


covariate <- rbind.data.frame(cbind.data.frame(covariate, "Taxon"= "B"),
                              cbind.data.frame(covariate, "Taxon" = "E"),
                              cbind.data.frame(covariate, "Taxon" = "G"),
                              cbind.data.frame(covariate, "Taxon"= "N"),
                              cbind.data.frame(covariate, "Taxon" = "T"))

input <- rbind.data.frame(Brachiopods,
                          Echinoderms,
                          Gastropods,
                          Nautiloids,
                          Trilobites)



## for sites as transects and subsites as replicates
mydata <- createPao(data = input,
                    unitcov = covariate)
summary(mydata)


# fit psi(section*Taxon)p(section*Taxon) model
mod1 <- occMod(model = list(psi~section*Taxon,
                            p~section*Taxon),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(section+Taxon)p(section+Taxon) model
mod2 <- occMod(model = list(psi~section+Taxon,
                            p~section+Taxon),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(locality*Taxon)p(locality*Taxon) model
mod3 <- occMod(model = list(psi~locality*Taxon,
                            p~locality*Taxon),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(locality+Taxon)p(locality+Taxon) model
mod4 <- occMod(model = list(psi~locality+Taxon,
                            p~locality+Taxon),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(section*Taxon*area)p(section*Taxon*area) model
mod5 <- occMod(model = list(psi~section*Taxon*area.sq,
                            p~section*Taxon*area.sq),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(section+Taxon+area)p(section*Taxon*area) model
mod6 <- occMod(model = list(psi~section+Taxon+area.sq,
                            p~section+Taxon+area.sq),
               data = mydata, type = "so",
               randinit = 5)


# fit psi(locality*Taxon*area)p(locality*Taxon*area) model
mod7 <- occMod(model = list(psi~locality*Taxon*area.sq,
                            p~locality*Taxon*area.sq),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(locality+Taxon+area)p(locality+Taxon+area) model
mod8 <- occMod(model = list(psi~locality+Taxon+area.sq,
                            p~locality+Taxon+area.sq),
               data = mydata, type = "so",
               randinit = 5)


# fit psi(Taxon*area)p(Taxon*area) model
mod9 <- occMod(model = list(psi~Taxon*area.sq,
                            p~Taxon*area.sq),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(Taxon+area)p(Taxon+area) model
mod10 <- occMod(model = list(psi~Taxon+area.sq,
                            p~Taxon+area.sq),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(Taxon)p(Taxon)
mod11 <- occMod(model = list(psi~Taxon,
                             p~Taxon),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon*section)p(Taxon)
mod12 <- occMod(model = list(psi~Taxon*section,
                             p~Taxon),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon)p(Taxon*section)
mod13 <- occMod(model = list(psi~Taxon,
                             p~Taxon*section),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon*section)p(section)
mod14 <- occMod(model = list(psi~Taxon*section,
                             p~section),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(section)p(1)
mod15 <- occMod(model = list(psi~section,
                             p~1),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon)p(1)
mod16 <- occMod(model = list(psi~Taxon,
                             p~1),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon*section)p(1)
mod17 <- occMod(model = list(psi~Taxon*section,
                             p~1),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(1)p(Taxon*section)
mod18 <- occMod(model = list(psi~1,
                             p~Taxon*section),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(1)p(Taxon)
mod19 <- occMod(model = list(psi~1,
                             p~Taxon),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(1)p(section)
mod20 <- occMod(model = list(psi~1,
                             p~section),
                data = mydata, type = "so",
                randinit = 5)


# fit psi(Taxon+section)p(1)
mod21 <- occMod(model = list(psi~Taxon+section,
                             p~1),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(1)p(Taxon+section)
mod22 <- occMod(model = list(psi~1,
                             p~Taxon+section),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(1)p(Taxon*facies)
mod23 <- occMod(model = list(psi~1,
                             p~Taxon*facies),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(1)p(Taxon+facies)
mod24 <- occMod(model = list(psi~1,
                             p~Taxon+facies),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon*facies)p(1)
mod25 <- occMod(model = list(psi~Taxon*facies,
                             p~1),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon+facies)p(1)
mod26 <- occMod(model = list(psi~Taxon+facies,
                             p~1),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon+facies)p(Taxon+facies)
mod27 <- occMod(model = list(psi~Taxon+facies,
                             p~Taxon+facies),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon*facies)p(Taxon*facies)
mod28 <- occMod(model = list(psi~Taxon*facies,
                             p~Taxon*facies),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon*facies)p(Taxon+facies)
mod29 <- occMod(model = list(psi~Taxon*facies,
                             p~Taxon+facies),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon+facies)p(Taxon*facies)
mod30 <- occMod(model = list(psi~Taxon+facies,
                             p~Taxon*facies),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(facies*area)p(facies*area) model
mod31 <- occMod(model = list(psi~facies*area.sq,
                            p~facies*area.sq),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(locality*Taxon*facies)p(locality*Taxon*facies) model
mod32 <- occMod(model = list(psi~locality+facies,
                            p~locality+facies),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(section)p(Taxon) model
mod33 <- occMod(model = list(psi~section,
                             p~Taxon),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(section)p(Taxon+section) model
mod34 <- occMod(model = list(psi~section,
                             p~Taxon+section),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(section)p(Taxon+section) model
mod35 <- occMod(model = list(psi~section+Taxon,
                             p~section),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(area)p(area) model
mod36 <- occMod(model = list(psi~area.sq,
                             p~area.sq),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(area+Taxon)p(area) model
mod37 <- occMod(model = list(psi~area.sq+Taxon,
                             p~area.sq),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(area)p(area+Taxon) model
mod38 <- occMod(model = list(psi~area.sq,
                             p~area.sq+Taxon),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(area*Taxon)p(area) model
mod39 <- occMod(model = list(psi~area.sq*Taxon,
                             p~area.sq),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(area)p(area*Taxon) model
mod40 <- occMod(model = list(psi~area.sq,
                             p~area.sq*Taxon),
                data = mydata, type = "so",
                randinit = 5)

# fit psi()p(area) model
mod41 <- occMod(model = list(psi~1,
                             p~area.sq),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(area)p() model
mod42 <- occMod(model = list(psi~area.sq,
                             p~1),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(Taxon)p(area) model
mod43 <- occMod(model = list(psi~Taxon,
                             p~area.sq),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(area)p(Taxon) model
mod44 <- occMod(model = list(psi~area.sq,
                             p~Taxon),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(section*Taxon+area)p(section*Taxon+area) model
mod45 <- occMod(model = list(psi~section*Taxon+area.sq,
                            p~section*Taxon+area.sq),
               data = mydata, type = "so",
               randinit = 5)

# fit psi(section+Taxon+area)p(section*Taxon+area) model
mod46 <- occMod(model = list(psi~section+Taxon+area.sq,
                             p~section*Taxon+area.sq),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(section*Taxon+area)p(section+Taxon+area) model
mod47 <- occMod(model = list(psi~section*Taxon+area.sq,
                             p~section+Taxon+area.sq),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(section+Taxon*area)p(section+Taxon+area) model
mod48 <- occMod(model = list(psi~section+Taxon*area.sq,
                             p~section+Taxon+area.sq),
                data = mydata, type = "so",
                randinit = 5)

# fit psi(section+Taxon+area)p(section+Taxon*area) model
mod49 <- occMod(model = list(psi~section+Taxon+area.sq,
                             p~section+Taxon*area.sq),
                data = mydata, type = "so",
                randinit = 5)

models<-list(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41,mod42,mod43,mod44,mod45,mod46,mod47,mod48,mod49)
results<-createAicTable(models)
summary(results)
write.table(summary(results), file = "AIC.txt", sep = ",", quote = FALSE, row.names = F)
##
save.image("//kant/nhm-sfs-u2/franzif/paper3/occmod.RData")
load("//kant/nhm-sfs-u2/franzif/paper3/occmod.RData")
## 
##
##
## model results for the third ranked model are extracted
## these results are below compared to naïve occupancy estimates
# coef(mod1,"psi") #regression coefficients for occupancy
# coef(mod1,"p") #regression coefficients for detection

## real estimates
m1fit.psi <- cbind.data.frame(covariate, fitted(mod1,"psi"))
m1fit.psi.unique <- unique(m1fit.psi[,c('section', 'Taxon', 'est', 'se', 'lower_0.95', 'upper_0.95')])
m1fit.psi.unique.a <- unite(m1fit.psi.unique[,c(2,1)], "group")
m1fit.psi.unique <- cbind.data.frame(m1fit.psi.unique, m1fit.psi.unique.a)

m1fit.p <- cbind.data.frame(covariate, fitted(mod1,"p"))
m1fit.p.unique <- unique(m1fit.p[,c('section', 'Taxon', 'est', 'se', 'lower_0.95', 'upper_0.95')])
m1fit.p.unique.a <- unite(m1fit.p.unique[,c(2,1)], "group")
m1fit.p.unique <- cbind.data.frame(m1fit.p.unique, m1fit.p.unique.a)

m1fit.p <- fitted(mod1,"p")
m1fit <- rbind.data.frame(cbind.data.frame(m1fit.psi.unique, "Parameter"="psi"),
                          cbind.data.frame(m1fit.p.unique, "Parameter"="p"))


m1fit <- m1fit[!(m1fit$Taxon=="E" & m1fit$section =="Low"),]
m1fit <- m1fit[!(m1fit$Taxon=="E" & m1fit$section =="Up"),]
m1fit <- m1fit[!(m1fit$Taxon=="G" & m1fit$section =="Low"),]


# ggplot(m1fit.psi.unique, aes(x = section, y = est))+
#   geom_point()+
#   geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95))+
#   facet_grid(. ~ Taxon) +
#   theme(axis.text.x = element_text(angle = 45))
# 
# ggplot(m1fit.p.unique, aes(x = section, y = est))+
#   geom_point()+
#   geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95))+
#   facet_grid(. ~ Taxon) +
#   theme(axis.text.x = element_text(angle = 45))
# 
# tiff("mod1.tiff",
#        res = 300,
#        units = "cm",
#        width = 15,
#        height = 8,
#        pointsize = 9)

ggplot(m1fit, aes(x = section, y = est))+
  geom_point()+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95), width=.3)+
  facet_grid(Parameter ~ Taxon) +
  theme(axis.text.x = element_text(angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# dev.off()
###
###
###
## model results for the preferred model are extracted and plotted
m5fit.psi <- cbind.data.frame(covariate, fitted(mod5,"psi"))
m5fit.psi.unique <- unique(m5fit.psi[,c('section', 'Taxon', 'area.sq', 'est', 'se', 'lower_0.95', 'upper_0.95')])
m5fit.psi.unique.a <- unite(m5fit.psi.unique[,c(2,1)], "group")
m5fit.psi.unique <- cbind.data.frame(m5fit.psi.unique, m5fit.psi.unique.a)

m5fit.p <- cbind.data.frame(covariate, fitted(mod5,"p"))
m5fit.p.unique <- unique(m5fit.p[,c('section', 'Taxon', 'area.sq', 'est', 'se', 'lower_0.95', 'upper_0.95')])
m5fit.p.unique.a <- unite(m5fit.p.unique[,c(2,1)], "group")
m5fit.p.unique <- cbind.data.frame(m5fit.p.unique, m5fit.p.unique.a)

m5fit.p <- fitted(mod5,"p")
m5fit <- rbind.data.frame(cbind.data.frame(m5fit.psi.unique, "Parameter"="psi"),
                          cbind.data.frame(m5fit.p.unique, "Parameter"="p"))
m5fit$area.sq <- as.character(m5fit$area.sq)

m5fit <- m5fit[!(m5fit$Taxon=="E" & m5fit$section =="Low"),]
m5fit <- m5fit[!(m5fit$Taxon=="E" & m5fit$section =="Up"),]
m5fit <- m5fit[!(m5fit$Taxon=="G" & m5fit$section =="Low"),]

cbp2 <- c(# "#000000",
          "#E69F00",
          "#56B4E9",
          "#F0E442",
          "#009E73",
          "#0072B2",
          "#D55E00",
          "#CC79A7")

m5fit$section <- factor(m5fit$section, levels = c("Low", "Up", "HJ", "JUAB"))
m5fit$area.sq <- factor(m5fit$area.sq, levels = c("625","1250","1875","2500","3750","5625"))

plot5section <- ggplot(m5fit, aes(x = section, y = est, group = area.sq, color = area.sq))+
  geom_point(position = position_dodge(width = 0.4))+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95, color=area.sq), width=.3 , na.rm=T, position =position_dodge(width = 0.4))+
  facet_grid(Parameter ~ Taxon) +         ## first element: Rows, second element: Columns
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
        # legend.position = "bottom",
        # legend.box = "vertical")+
  # guides(colour = guide_legend(nrow = 1))+ ## defines that the legend is shown in one horizontal line
  scale_colour_manual(values=cbp2,
                      name="Area in"~cm^2)


plot5area <- ggplot(m5fit, aes(x = area.sq, y = est, group = section, color = section))+
  geom_point(position = position_dodge(width = 0.4))+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95, color=section), width=.3 , na.rm=T, position =position_dodge(width = 0.4))+
  facet_grid(Parameter ~ Taxon) +         ## first element: Rows, second element: Columns
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
        # legend.position = "bottom",
        # legend.box = "vertical")+
  # guides(colour = guide_legend(nrow = 1))+ ## defines that the legend is shown in one horizontal line
  scale_colour_manual(values=cbp2,
                      name="Section")


tiff("mod5.double.tiff",
     res = 300,
     units = "cm",
     width = 18.75,
     height = 15,
     pointsize = 7)

ggarrange(plot5section, plot5area,
          labels = c("A","B"),
          ncol = 1, nrow = 2)

dev.off()

##
###
##
###
##


## 
## model results from the second best model are extracted and plotted
m46fit.psi <- cbind.data.frame(covariate, fitted(mod46,"psi"))
m46fit.psi.unique <- unique(m46fit.psi[,c('section', 'Taxon', 'area.sq', 'est', 'se', 'lower_0.95', 'upper_0.95')])
m46fit.psi.unique.a <- unite(m46fit.psi.unique[,c(2,1)], "group")
m46fit.psi.unique <- cbind.data.frame(m46fit.psi.unique, m46fit.psi.unique.a)

m46fit.p <- cbind.data.frame(covariate, fitted(mod46,"p"))
m46fit.p.unique <- unique(m46fit.p[,c('section', 'Taxon', 'area.sq', 'est', 'se', 'lower_0.95', 'upper_0.95')])
m46fit.p.unique.a <- unite(m46fit.p.unique[,c(2,1)], "group")
m46fit.p.unique <- cbind.data.frame(m46fit.p.unique, m46fit.p.unique.a)

m46fit.p <- fitted(mod46,"p")
m46fit <- rbind.data.frame(cbind.data.frame(m46fit.psi.unique, "Parameter"="psi"),
                          cbind.data.frame(m46fit.p.unique, "Parameter"="p"))
m46fit$area.sq <- as.character(m46fit$area.sq)

m46fit <- m46fit[!(m46fit$Taxon=="E" & m46fit$section =="Low"),]
m46fit <- m46fit[!(m46fit$Taxon=="E" & m46fit$section =="Up"),]
m46fit <- m46fit[!(m46fit$Taxon=="G" & m46fit$section =="Low"),]

m46fit$section <- factor(m46fit$section, levels = c("Low", "Up", "HJ", "JUAB"))
m46fit$area.sq <- factor(m46fit$area.sq, levels = c("625","1250","1875","2500","3750","5625"))

plot46section <- ggplot(m46fit, aes(x = section, y = est, group = area.sq, color = area.sq))+
  geom_point(position = position_dodge(width = 0.4))+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95, color=area.sq), width=.3 , na.rm=T, position =position_dodge(width = 0.4))+
  facet_grid(Parameter ~ Taxon) +         ## first element: Rows, second element: Columns
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
        # legend.position = "bottom",
        # legend.box = "vertical")+
  # guides(colour = guide_legend(nrow = 1))+ ## defines that the legend is shown in one horizontal line
  scale_colour_manual(values=cbp2,
                      name="Area in"~cm^2)

plot46area <- ggplot(m46fit, aes(x = area.sq, y = est, group = section, color = section))+
  geom_point(position = position_dodge(width = 0.4))+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95, color=section), width=.3 , na.rm=T, position =position_dodge(width = 0.4))+
  facet_grid(Parameter ~ Taxon) +         ## first element: Rows, second element: Columns
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  #       legend.position = "bottom",
  #       legend.box = "vertical")+
  # guides(colour = guide_legend(nrow = 1))+ ## defines that the legend is shown in one horizontal line
  scale_colour_manual(values=cbp2,
                      name="Section")


tiff("mod46.double.tiff",
     res = 300,
     units = "cm",
     width = 18.75,
     height = 15,
     pointsize = 7)

ggarrange(plot46section, plot46area,
          labels = c("A","B"),
          ncol = 1, nrow = 2)

dev.off()






## naïve occupancy estimates
nest.df <- cbind.data.frame(input, covariate)

Zuversicht <- function(x,y){
  abc <- nest.df %>%
    filter(Taxon %in% x) %>%
    filter(section %in% y)
  no.presences <- sum(colSums(abc[1:8],na.rm = T))
  no.0 <- sum(!is.na(abc[1:8]))
  nest <- no.presences/no.0
  
}
# 
# for (j in Taxon) {
#   for (i in section){
#     
#     nest <- Zuversicht(j,i)
#     print(nest)
#     # result <- rbind(result,cbind("section"=i,"Taxon"=j,"estimate"=nest))
#     MuhKuh <- data.frame(cbind("Taxon"=j,"section"=i,"estimate"=nest))
#     print(MuhKuh)
#     nest.result <- merge(result, MuhKuh, by=intersect(names(result),names(MuhKuh)),all=TRUE)
#     
#   }
# }
# nest.result

## and their respective confidence intervals are estimated through bootstrapping
set.seed(13579)   # set a seed for consistency/reproducibility
Optimist.boot <- function(x,y){
  abc <- nest.df %>%
    filter(Taxon %in% x) %>% 
    filter(section %in% y)
  vec <- as.vector(as.matrix(abc[1:8]))
  vec <- vec[!is.na(vec)] ## vector that contains all presences/absences that were recorded from the subsites, and that we sample from with replacement
  no.0 <- sum(!is.na(abc[1:8])) ## amount of samples we want to take from our presences/absences (with replacement), equals the number of sampled subsites
  B <- 1000 ## number of bootstrap samples we like to take
  
  boot <- matrix( sample(vec, size= B*no.0, replace=TRUE), ncol=B, nrow=no.0)
  
  boot.presences <- colSums(boot,na.rm = T)
  boot.distribution <- boot.presences/no.0
  CIs <- cbind.data.frame("lCL"=quantile(boot.distribution, prob=0.025), "uCL"= quantile(boot.distribution, prob=0.975))
  
}


Taxon <- c("B", "E","G", "N", "T")
section <- c("Low", "Up", "HJ", "JUAB")
nest.result <- data.frame("section"=NA,"Taxon"=NA,"estimate"=NA)

for (j in Taxon) {
  for (i in section){
    
    nest <- Zuversicht(j,i)
    boots <- Optimist.boot(j,i)
    print(nest)
    # result <- rbind(result,cbind("section"=i,"Taxon"=j,"estimate"=nest))
    MuhKuh <- data.frame(cbind("Taxon"=j,"section"=i,"estimate"=nest, "lCL"=boots))
    print(MuhKuh)
    nest.result <- merge(nest.result, MuhKuh, by=intersect(names(nest.result),names(MuhKuh)),all=TRUE)
    
  }
}
nest.result <- nest.result[1:20,]

# tiff("nest.tiff",
#      res = 300,
#      units = "cm",
#      width = 15,
#      height = 5,
#      pointsize = 9)
# 
# ggplot(nest.result, aes(x = section, y = estimate))+
#   geom_point()+
#   geom_errorbar(aes(ymin=lCL.lCL, ymax=lCL.uCL), width=.3)+
#   facet_grid(. ~ Taxon) +
#   theme(axis.text.x = element_text(angle = 45),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank())
# 
# dev.off()


df <- m1fit %>% select(section, Taxon, est, lower_0.95, upper_0.95, Parameter)
nest.result1 <- nest.result %>% 
  rename(
     lower_0.95 = lCL.lCL,
    upper_0.95 = lCL.uCL,
    est = estimate
  ) %>%
  mutate("Parameter"="psi")

dat <- rbind.data.frame(cbind.data.frame(df, ID ="modelled"),
                        cbind.data.frame(nest.result1, ID="naïve estimate"))
dat$section <- factor(dat$section, levels = c("Low", "Up", "HJ", "JUAB"))

## estimates from the third best model and the naïve estimates are plotted in comparison.
tiff("pref.plus.naïve.tiff",
     res = 300,
     units = "cm",
     width = 15,
     height = 10,
     pointsize = 9)

ggplot(dat, aes(x = section, y = est, group = ID, color=ID))+
  geom_point(position = position_dodge(width = 0.3))+
  geom_errorbar(aes(ymin=lower_0.95, ymax=upper_0.95, color=ID), width=.3, na.rm=T, position =position_dodge(width = 0.3))+
  facet_grid(Parameter ~ Taxon) +
  theme(axis.text.x = element_text(angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical")+
  guides(colour = guide_legend(nrow = 1))+ ## defines that the legend is shown in one horizontal line
  scale_color_manual(values = c("#009E73", "#D55E00"))

dev.off()
