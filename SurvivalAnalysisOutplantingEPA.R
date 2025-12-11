library(dplyr)
library(lubridate)
library(tidyverse)
library(survival)
library(survminer)
library(ggsurvfit)
library(gtsummary)
library(car)
library(knitr)

setwd("C:/Users/kcoles/ABG Dropbox/Kelly Coles/EPA")

data <- read.csv("Outplanting/SurvivalPkgData.csv")
summary(data)
###Note that data have been transformed so that values are switched. Survival analysis
###considers a value of 1 to indicate that an event, i.e., death, occurred. So 1
###here indicates death, and a 0 indicates survival.

data <- subset(data, select = c(Tag.Number:Repr4))
data$T43 <- as.numeric(data$T43)
data$T92 <- as.numeric(data$T92)
data$T162 <- as.numeric(data$T162)
data$T238 <- as.numeric(data$T238)

#Adding a status column - whether death was observed or not
data <- data %>% mutate(status = case_when(
  T43 == 1 | 
  T92 == 1 |
  T162 == 1 |
  T238 == 1 ~ 1, TRUE ~ 0))

#Adding a time of death or censorship column
data <- data %>% mutate(Time = case_when(
  T43 == 1 ~ 43, 
  T92 == 1 ~ 92,
  T162 == 1 ~ 162,
  T238 == 1 ~ 238,
  TRUE ~ 238))

#Adding a column for whether plant was ever observed reproductive
data <- data %>% mutate(Flwr = case_when(
                        Repr1 == "Y" |
                        Repr2 == "Y" |
                        Repr3 == "Y" |
                        Repr4 == "Y" ~ "Y", T ~ "N"
))

#We need to drop tag number 695 because it was never planted.
data <- data %>% filter(Tag.Number!="695")

##########################################################################  
##We do not actually need to make the data long but I want to preserve this code

#data_long <- data %>% pivot_longer(cols=c(T43:T238),
                    #names_to='Days',
                    #values_to="Death")
#data_long <- data_long %>% pivot_longer(cols=c(Repr1:Repr4),
                                 #  names_to='Event',
                                #   values_to="ReproductiveStat")

#data_long$Days <- gsub('T','',data_long$Days)
#data_long$Days <- as.numeric(data_long$Days)
#data_long$Death <- as.numeric(data_long$Death)
##########################################################################



#creating the survival object, to be used as a response in models
Surv(data$Time, data$status)

#create survival curves
s1 <- survfit(Surv(Time, status) ~ 1, data = data)
str(s1)

#plotting survival curve
survfit2(Surv(Time, status) ~ 1, data = data) %>% 
  ggsurvfit() +
  labs(
    x = "Time in days",
    y = "Overall survival probability"
  )+ 
  add_confidence_interval()+
  add_risktable()

#estimating percentage that survived to 200 days
summary(survfit(Surv(Time, status) ~ 1, data = data), times = 200)

##Comparing Survival between Groups - Site
survdiff(Surv(Time, status) ~ Site, data = data) #no difference

##Comparing Survival between Groups - Species
survdiff(Surv(Time, status) ~ Species, data = data) #significant difference

##Comparing Survival between Groups - Reproductive status
survdiff(Surv(Time, status) ~ Flwr, data = data) #significant difference

#Fitting a regression model
fit1 <- coxph(Surv(Time, status) ~ Site, data = data)
print(fit1)
Anova(fit1)

fit2 <- coxph(Surv(Time, status) ~ Species, data = data) 
summary(fit2)
Anova(fit2)

fit3 <- coxph(Surv(Time, status) ~ Flwr, data = data) 
summary(fit3)
Anova(fit3)

##table of results
broom::tidy(
  coxph(Surv(Time, status) ~ Site, data = data), 
  exp = TRUE
) %>% 
  kable()

#is assumption of proportional hazards violated?
mv_fit <- coxph(Surv(Time, status) ~ Site, data = data)
cz <- cox.zph(mv_fit)
print(cz) #no
plot(cz)

mv_fit2 <- coxph(Surv(Time, status) ~ Species, data = data)
cz2 <- cox.zph(mv_fit2)
print(cz2) #it does violate assumption
plot(cz2)

mv_fit3 <- coxph(Surv(Time, status) ~ Flwr, data = data)
cz3 <- cox.zph(mv_fit3)
print(cz3) #no
plot(cz3)

survfit2(Surv(Time, status) ~ Site, data = data) |> 
  ggsurvfit(linetype_aes = F,linewidth = 1, show.legend=F) +
  scale_ggsurvfit()+
  add_confidence_interval()+
  labs(x = "Time (in days)")+
  theme(legend.position = "right",
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_text(face = "bold"))+
  labs(fill = "Site") +
  scale_fill_discrete(labels=c("Cleared + Burned", "Cleared + Mowed", "Control"))

LegendTitle = "Reproductive Status"
survfit2(Surv(Time, status) ~ Flwr, data = data) |> 
  ggsurvfit(linetype_aes = F,linewidth = 1, show.legend=F) +
  scale_ggsurvfit()+
  add_confidence_interval()+
  labs(x = "Time (in days)")+
  theme(legend.position = "right",
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_text(face = "bold"))+
  labs(fill = "Reproductive") +
  scale_fill_discrete(labels=c("No", "Yes"))



ggplot(data = data, aes(x = as.factor(Species),  fill = as.factor(status))) + 
  geom_bar(color="black")+
  labs(x = "Species", fill = "Status") +
  scale_fill_manual(values = c("ivory", "gray1"), labels = c("Alive", "Dead"))+
  theme_classic()

#looking at percent mortality overall
data %>% group_by(status) %>% summarise(count = n())
375/1620 #23%

#looking at percent mortality by reproductive status overall
data %>% group_by(Flwr, status) %>% summarise(count = n())
#  Flwr  status count
#<chr>  <dbl> <int>
#  1 N        0  1169
#2 N          1   373
#3 Y          0    75
#4 Y          1     2

373/1542 #24% mortality for nonreproductives vs
2/77 #3% for reproductive

#looking at percent mortality by species
data %>% group_by(Species, status) %>% summarise(count = n())
49/397 #12% total Andropogon and 49 died
33/139 #24% ; 139 total Antaenantia and 33 died
104/(104+489) #18% Balduina died
6/89 #7% Liatris died
140/245 #57% Oxypolis died
43/156 #28% Rhynchospora died