
################################################################################
#                                                                              #
# Purpose:       Code to reproduce the analyses and figures of PEru SARs-CoV-2 #
#                study 2019-2021 for AJTMH                                     #
#                                                                              #
# Author:        Andres Moreira-Soto                                           #
# Contact:       andres.moreira-soto@charite.de                                #
# Client:                                                                      #
#                                                                              #
# Code created:  2022-01-12                                                    #
# Last updated:  2022-01-12                                                    #
# Source:                                                                      #
#                                                                              #
# Comment:                                                                     #
#                                                                              #
################################################################################

########################### Load required libraries ############################


library(xlsx)
library(tidyverse)
library(summarytools)

#Avidity index higher than 50% were considered as high antibody avidity; between 31 to 49% as intermediate avidity, and below 30% as low avidity.
#IgG antibody avidity of SARS-CoV-2 antibody and PRNT positive samples. Dashed lines denote the ratio thresholds of >60% (positive) and <40% (negative); results between these values are considered borderline, as defined by the manufacturer EUROIMMUN (https://www.euroimmun.com). The median avidity is indicated by bar, mean by cross, quartiles by boxes, and minimum and maximum by whiskers. D) SARS-CoV-2 IgA ELISA reactivity of samples showing high, borderline or low IgG avidity.
#An avidity index for each sample was calculated by dividing the extinction of the sample with urea treatment by the extinction of the sample without the urea treatment and a final multiplication of 100

data <- read.xlsx("PATHTFILE/NAME.xlsx")

range(data$date)

data$datacategory <- ifelse (data$date < as.Date("2020-05-14", "%Y-%m-%d") & data$date > as.Date("2020-04-01", "%Y-%m-%d"), "T1",
                              ifelse (data$date > as.Date("2020-05-14", "%Y-%m-%d"), "T2", "T0"))

dataselected <- data %>% select (4,6,30,31,40,42,43,44,50,51,52,53,54) #%>% mutate_if(is.character,as.numeric) 

dataurea <- dataselected %>% select(2,3,4,5,9,10,11,12,13) %>% mutate (SARSNCPavidityindex = (SARS.CoV.2.NCP.IgGwith.Urea / SARS.CoV.2.NCP.IgGwithout.Urea) *100, SARSSavidityindex = (SARS.CoV.2..IgGwith.Urea / SARS.CoV.2..IgGwithout.Urea) *100  )

dataurea2 <- dataurea %>% select(date,prntresult,positivity,datacategory,SARSNCPavidityindex,SARSSavidityindex) %>% filter (positivity== "p")

dataurea3 <- dataurea2 %>% select (date,prntresult,datacategory,SARSNCPavidityindex,SARSSavidityindex) 

dataurea4 <- dataurea3 %>% pivot_longer (cols = !c(date, datacategory, prntresult), names_to = "test", values_to = "Avidityindex")

ggplot(dataurea4) +
  aes(x = date, y = Avidityindex, group=test, color= test) +
  geom_smooth(alpha=0.1, method = "lm") +
  #geom_point(aes(color=test), alpha= 0.8) +
  geom_point (aes(shape=prntresult), size=2)+
  scale_color_manual(values= c("gray", "darkred"))+
  scale_y_continuous(limits = c(0, 120))+
  geom_hline(yintercept = 40)+
  geom_hline(yintercept = 60)+
  scale_x_date(date_breaks = "1 month", date_labels = "%m")+
  facet_wrap(test~prntresult)+
  theme_Publication (22)+
  theme(axis.line = element_line(color = "black", 
                             linetype = "solid"))+
  theme(axis.text.x = element_text(angle = 0))


ggplot(dataurea4) +
  aes(x = prntresult, y = Avidityindex, color= test) +
  geom_boxplot() +
  geom_point(aes (color= test, shape= prntresult), position = position_jitterdodge())+
  #geom_point (aes(shape=prntresult), position = position_jitterdodge())+
  scale_color_manual(values= c("gray", "darkred"))+
  geom_hline(yintercept = 40)+
  geom_hline(yintercept = 60)+
  scale_y_continuous(limits = c(0, 120))+
  #scale_x_date(date_breaks = "1 month", date_labels = "%m-%y" , limits=as.Date(c("2020-01-01","2020-07-01")))+
  #facet_wrap(~positivity)+
  theme_Publication (22)+
  theme(axis.line = element_line(color = "black", 
                                 linetype = "solid"))+
  theme(axis.text.x = element_text(angle = 90))

#category of time

NCP1 <- filter(dataurea4, datacategory == "T0" & test== "SARSNCPavidityindex")$Avidityindex
NCP2 <- filter(dataurea4, datacategory == "T1" & test== "SARSNCPavidityindex")$Avidityindex
NCP3 <- filter(dataurea4, datacategory == "T2" & test== "SARSNCPavidityindex")$Avidityindex

S1 <- filter(dataurea4, datacategory == "T0" & test== "SARSSavidityindex")$Avidityindex
S2 <- filter(dataurea4, datacategory == "T1" & test== "SARSSavidityindex")$Avidityindex
S3 <- filter(dataurea4, datacategory == "T2" & test== "SARSSavidityindex")$Avidityindex

t.test(NCP1, NCP2)
t.test(NCP1, NCP3)
t.test(NCP2, NCP3)

t.test(S1, S2)
t.test(S1, S3)
t.test(S2, S3)


#category of PRNT positivity 
NCPn <- filter(dataurea4, prntresult == "n" & test== "SARSNCPavidityindex")$Avidityindex
NCPp <- filter(dataurea4, prntresult == "p" & test== "SARSNCPavidityindex")$Avidityindex
t.test(NCPn, NCPp)

Sn <- filter(dataurea4, prntresult == "n" & test== "SARSSavidityindex")$Avidityindex
Sp <- filter(dataurea4, prntresult == "p" & test== "SARSSavidityindex")$Avidityindex
t.test(Sn, Sp)

# some descriptive statistics 

freq(dataurea3$prntresult)
ctable(dataurea2$date, dataurea2$positivity)
ctable(dataurea2$datacategory, dataurea2$prntresult)

chisq.test (dataurea2$datacategory, dataurea2$prntresult)

###########################
#analises of the data regarding prnt negatives

data2 <- data %>% mutate (SARSNCPavidityindex = (SARS.CoV.2.NCP.IgGwith.Urea / SARS.CoV.2.NCP.IgGwithout.Urea) *100, 
                          SARSSavidityindex = (SARS.CoV.2..IgGwith.Urea / SARS.CoV.2..IgGwithout.Urea) *100  )

data2$Inhibition <- as.numeric(data2$Inhibition)

data3 <- data2 %>% select(Concentrationmagcopy, Inhibition, prntresult, date, 
                          SARSNCPavidityindex, SARSSavidityindex, datefactor ,datacategory) %>%
                  filter (datacategory == "T0") %>%
                  filter (Concentrationmagcopy > 0.99 & Inhibition > 30) %>% 
                  filter(SARSSavidityindex > 60 & SARSNCPavidityindex > 60) 


ggplot(data3) +
 aes(x = SARSSavidityindex, y = SARSNCPavidityindex, colour = Inhibition, size = Concentrationmagcopy) +
 geom_point(shape = "bullet") +
 scale_color_gradient(low = "gray", high = "darkred") +
 ggthemes::theme_par() +
 theme(legend.position = "bottom") +
 facet_wrap(vars(datefactor))

