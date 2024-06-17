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

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor)

# Read in the data from the Excel file
dataclean <- read_excel ("dataclean.xlsx")
dataclean$date <- as.Date (dataclean$date)

####################### Code for Figure 1 Panels #################################

# Create a table for seropositivity data
tableserop <- as.data.frame(ftable(xtabs(~datefactor + positivity + City, data=dataclean)))

# Reshape the data to have separate columns for positivity
tableserop2 <- tableserop %>% pivot_wider(names_from = positivity, values_from = Freq)

# Calculate total and percentage positivity
tableserop2$total <- tableserop2$n + tableserop2$p
tableserop2$percpos <- tableserop2$p * 100 / tableserop2$total

# Add mortality data to the table
datamort <- read.csv("C:/PATHTOMORTALITYDATA.csv", sep= ";")

# Clean the column names
datamortclean <- datamort %>% clean_names()

# Convert date columns to Date type
datamortclean$fecha_corte <- as.Date(as.character(datamortclean$fecha_corte), format="%Y%m%d")
datamortclean$fecha_fallecimiento <- as.Date(as.character(datamortclean$fecha_fallecimiento), format="%Y%m%d")

# Filter for specific department and create date factor
datacleanonlydep <- datamortclean %>% filter(departamento %in% c("LIMA"))
datacleanonlydep$datefactor <- cut.Date(datacleanonlydep$fecha_fallecimiento, breaks="month")

# Summarize mortality data
datacleanonlydep2 <- datacleanonlydep %>% 
  select(datefactor, departamento) %>% 
  group_by(datefactor, departamento) %>% 
  summarise(countmort = n())

# Combine seropositivity and mortality dataframes
dataserumandmort <- tableserop2 %>% 
  left_join(datacleanonlydep2, by = c("datefactor" = "datefactor", "City" = "departamento"))

# Filter for Lima
dataserumandmortlima <- dataserumandmort %>% filter(City == "LIMA")

# Correlation tests between countmort and percpos for different cities
cor.test(dataserumandmortlima$countmort, datacleanonlydep2$percpos, method="pearson", use="pairwise.complete.obs")

# Read additional mortality data
datamortdep <- read_xlsx("C:/PATHTOMORTALITYDATA.xlsx", sheet="mortexcessdep")

# Create date factor for the additional mortality data
datamortdep$datefactor <- cut.Date(as.Date(datamortdep$from), breaks="month")

# Select relevant columns from additional mortality data
datamortdepfilter <- datamortdep %>% select(department, observed, expected, cihigh, cilow, excess, additionalpercent, datefactor)

# Combine all data into one dataframe
dataallcomplete <- dataserumandmort %>% 
  left_join(datamortdep, by = c("datefactor" = "datefactor", "City" = "department"))

# Calculate percentage of additional deaths
dataallcomplete$additionalpercent2 <- dataallcomplete$additionalpercent * 100

# Plot percent positives serology
ggplot(data=dataallcomplete, mapping=aes(x=as.Date(datefactor))) +
  geom_smooth(aes(y=percpos), span=0.55, color="darkred") +
  geom_point(aes(y=percpos)) +
  scale_x_date(date_breaks="1 month", date_labels="%m", limits=as.Date(c("2020-01-01", "2021-06-01"))) +
  labs(x="Date", y="Seropositivity per month (%)") +
  scale_fill_viridis_d(option="cividis", direction=1) +
  facet_wrap(vars(City), ncol=1) +
  scale_y_continuous(name="Seropositivity per month (%)", breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))

# Plot excess mortality
ggplot(data=dataallcomplete, mapping=aes(x=as.Date(datefactor))) +
  geom_smooth(aes(y=additionalpercent2), span=0.1, color="darkred") +
  labs(x="Date", y="Excess mortality (%)") +
  scale_fill_viridis_d(option="cividis", direction=1) +
  facet_wrap(vars(City), ncol=1) +
  scale_x_date(date_breaks="1 month", date_labels="%m")


# plot data in time cpass versus maglumi -----------------------------------------------------------------------------------------------

datanona <- dataclean %>% drop_na (Inhibition) %>% mutate (Inhibitionzero= ifelse(Inhibition <0 , 0, Inhibition))

ggplot(datanona, aes(x = Inhibitionzero, y =log10(Concentrationmagcopy), color= positivity))  + 
  geom_point(position= "jitter") +
  geom_point(aes(shape=prntresult), size=4) +
  #geom_smooth(method = "lm")+
  geom_abline(intercept = 0.12829524, slope = 0.012794772)+
  xlab("Log CLIA oncentration")+ 
  ylab("sVNT percent inhibition")+
  scale_color_manual(values =c("grey", "darkred")) +
  geom_hline(yintercept= log(1), linetype="dashed")+
  geom_vline(xintercept=30, linetype="dashed")+
  #scale_y_continuous(breaks=seq(0,100,10))+
  scale_x_continuous(breaks=seq(0,100,10))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_line(color = NA, 
                                        linetype = "dashed", 
                                        size = 0), 
        panel.background = element_rect(fill = NA),  
        panel.border = element_rect(fill = NA))#+

# Calculate slope and intercept of line of best fit
coef(lm(log10(Concentrationmagcopy) ~ Inhibitionzero, data = datanona))

ggsave("maglumivscpassonly20202021.pdf")

# Data peru CLIA over time

ggplot(datanona, aes(x = date, y =Concentrationmagcopy, color=positivity))  + 
  geom_point(position= "jitter") +
  geom_point(aes(shape=prntresult), size=4) + 
  geom_smooth()+
  scale_color_manual(values =c("grey", "darkred")) +
  xlab("Sample reception time")+ 
  ylab("CLIA concentration")+
  geom_hline(yintercept=c(1), linetype="dashed")+
  scale_y_continuous(breaks=seq(0,100,10))+
  scale_x_date(date_breaks = "1 month", date_labels = "%m-%y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(panel.grid.major = element_line(color = NA, 
                                        linetype = "dashed", 
                                        size = 0), 
        panel.background = element_rect(fill = NA),  
        panel.border = element_rect(fill = NA))#+
