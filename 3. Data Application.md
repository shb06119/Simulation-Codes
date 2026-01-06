# 3. Data Application
We apply three modeling approahces (two-stage method, Bayesian hierarchical model, and spatial Bayesian hierarchical model) to empirical dataset. The dataset includes daily death counts, ozone concentration (8 hour moving average), and temperature (Kelvin) for 254 districts of South Korea from 2015 to 2021. We use 77 districts of Seoul, Incheon, and Gyeonggi among 254 districts. The area of Michuhol-gu was excluded from the analysis from 2015 to 2017 because no deaths were recorded during this period due to measurement errors. Likewise, the area of Bucheon-si was excluded for the year of 2015 for the same reason.
```{r, warning=FALSE, message=FALSE}
# Necessary packages
library(sf); library(gnm); library(data.table); library(lubridate); library(dplyr)
library(ggplot2); library(spdep); library(fields); library(leaflet)  
library(metafor); library(R2jags); library(parallel); library(devtools)
library(tidyverse); library(patchwork); library(coda); library(rjags)
library(spdep); library(spatialreg)
```
```{r}
# Data Uploading
death254_nostand = fread("/Users/nani13143553/Desktop/KAIST/Research/before CCO/death count before standardization 254.csv")[,-c(1)]

# Renaming each district from district code(e.g. 11011) to natural numbers(1-254)
death254_nostand$city_num <- as.numeric(gsub("_", "", death254_nostand$city_dist))
death254_nostand$region <- as.factor(death254_nostand$city_num)
death254_nostand$region <- as.numeric(death254_nostand$region)

# Eliminating non-existing areas contained in raw data
death1 = subset(death254_nostand, region < 52)
death2 = subset(death254_nostand, 52 < region)
death3 = subset(death2, region < 88)
death4 = subset(death254_nostand, region > 90)
death5 = subset(death4, region < 230)
death6 = subset(death254_nostand, region > 230)
death7 = subset(death6, region < 253)
death247 = rbind(death1, death3, death5, death7)
death247 = death247 %>% arrange(city_dist)
death247$region = as.factor(death247$city_dist)
death247$region = as.numeric(death247$region)

# Setting strata [area]-[year]-[month]-[day of the week]
year = year(death247$death_date)
month = str_pad(month(death247$death_date), 2, pad = "0")
dow = wday(death247$death_date)
md = paste(month, dow, sep = "_")
ymd = paste(year, month, dow, sep = "_")
death247$ymd = as.factor(ymd)
death247$zone_ymd = as.factor(paste(str_pad(death247$region, 2, pad = "0"), ymd, sep = "_"))


```
