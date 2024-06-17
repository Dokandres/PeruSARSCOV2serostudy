###############################################################################
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


#the goal of this script is to analize the whole datasets on seroprevalence data from https://serotracker.com/en/Explore

#load packages
{
library(tidyverse)
library(readxl)
library(xlsx)
library(janitor)
library(summarytools)
library (esquisse)
library(ggplot2)
library(DataExplorer)
library(plotly)
}
#### LOAD DATA ----------------------------------------------------------------
## Set working directory

datacomplete <- read.csv("SeroTracker_ Serosurveys Reporting Prevalence.csv", na.strings = "")

datacompletecleanperu <- datacomplete %>% clean_names() %>% 

filter(country %in% c("Peru"))

#clean and mofify the data 
#take all % out
datacompletecleanperu[] <- lapply(datacompletecleanperu, gsub, pattern='%', replacement='')
#convert to numeric
datacompletecleanperu <- type.convert(datacompletecleanperu, as.is = TRUE)
#convert dates 
datacompletecleanperu <- datacompletecleanperu %>% mutate_at(vars(publication_date,sampling_start_date,sampling_end_date), as.Date, format="%d. %B %Y")
datacompletecleanperu <- datacompletecleanperu %>% mutate_at(vars(date_created, last_modified_time), as.Date, format="%d.%m.%Y")

#datacompletecleanperu %>% create_report()

datacompletecleanperu$state_province <- ifelse(grepl("Lima", datacompletecleanperu$state_province, ignore.case = TRUE), "Lima", datacompletecleanperu$state_province)
datacompletecleanperu$state_province <- ifelse(datacompletecleanperu$state_province == "Coronel Portillo", "Ucayali", datacompletecleanperu$state_province )
datacompletecleanperu$state_province <- ifelse(datacompletecleanperu$state_province == "Jaen", "Cajamarca", datacompletecleanperu$state_province )

# factor for the risk of bias 

datacompletecleanperu$overall_risk_of_bias_jbi_archivefactor <- factor(datacompletecleanperu$overall_risk_of_bias_jbi_archive, levels = c("Unclear",  "Low", "Moderate", "High"))

#desired colors not to do manually
datacompletecleanperu <- datacompletecleanperu %>%
  mutate(Color = ifelse(grepl("Lima", state_province), "Lima", "Others"))


#plot the data

ggplot(datacompletecleanperu) +
 aes(x = sampling_start_date, y = serum_positive_prevalence, color= Color ) +
  geom_point(aes (x=sampling_start_date, shape =  overall_risk_of_bias_jbi_archivefactor, size=2)) +
  geom_point(aes (x=sampling_end_date, shape =  overall_risk_of_bias_jbi_archivefactor, size=2)) +
  geom_segment(aes(x = sampling_start_date, y = serum_positive_prevalence, xend = sampling_end_date, yend = serum_positive_prevalence),
               alpha = 0.5, size = 0.5) + # add a line segment between the points
  scale_color_manual(values=c("darkred", "darkgray"))+
  geom_smooth(method = "lm", se=FALSE, color="darkred")+
  labs(x = "Date", 
       y = "seroprevalence (%)",
       colour=" Region") +
  theme_Publication(18)+ 
  theme(legend.position = "none")

#plot the data interactively

p <- ggplot(datacompletecleanperu) +
  aes(x = sampling_start_date, y = serum_positive_prevalence, color= state_province ) +
  geom_point(aes (x=sampling_start_date, shape = "circle", size = overall_risk_of_bias_jbi_archive)) +
  geom_point(aes (x=sampling_end_date, shape = "circle", size = overall_risk_of_bias_jbi_archive)) +
  geom_segment(aes(x = sampling_start_date, y = serum_positive_prevalence, xend = sampling_end_date, yend = serum_positive_prevalence),
               alpha = 0.5, size = 0.5) + # add a line segment between the points
  scale_color_viridis_d(option = "cividis") +
  geom_smooth(method = "lm", se=FALSE, color="darkred")+
  labs(x = "Date", 
       y = "seroprevalence (%)",
       colour=" Region") +
  theme_Publication(18) 

ggplotly(p)
 


ggplot(datacompletecleanperu) +
 aes(x = publication_date, fill = overall_risk_of_bias_jbi_archive ) +
 geom_bar() +
 labs(x = "Country", 
 y = "Number of studies") +
 #coord_flip() +
theme_Publication(18) 


summarytools::ctable (datacompletecleanonlyafrica$country, datacompletecleanonlyafrica$serum_positive_prevalence, useNA= "no")


# only lima data

datacompletecleanperulima <- datacompletecleanperu %>% filter(str_detect(city, "Lima"))

ggplot(datacompletecleanperulima) +
  aes(x = publication_date, y = serum_positive_prevalence, color= overall_risk_of_bias_jbi_archive ) +
  geom_point(shape = "circle", size = 3) +
  scale_color_viridis_d(option = "cividis", direction = -1) +
  geom_smooth(method = "lm", se=TRUE, color="darkred")+
  scale_x_date(date_breaks = "1 month", date_labels = "%m-%y")+
  labs(x = "Date", 
       y = "seroprevalence (%)",
       colour=" Risk of bias") +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
  theme_Publication(18) 


datacompleteonlylima <- datacompletecleanperu %>% dplyr::filter(state_province== "Lima")
