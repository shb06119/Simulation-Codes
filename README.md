# Simulation-Code (No missing-data)
This is the code for the conditional poisson regression model based simulation which generates 20 datasets. Note that the original data includes two regions which have observational error. In particular, the region 57 and 86 have missing death counts for 2015-2017 and 2017, respectively. Therefore, we exclude year 2015-2017 from the 57th region and 2015 from the 86th region. Due to this treatment, simulating death count differs by region. This is why we need to carefully consider data-generating code. By the way, this code creates data for 2019-2022, thereby no need to consier missing data. However, we will consider the above missing data case after this.

```{r, warning=FALSE, message=FALSE}
library(sf); library(gnm); library(data.table); library(lubridate); library(dplyr)
library(ggplot2); library(spdep); library(fields); library(leaflet)  
library(metafor); library(R2jags); library(parallel); library(devtools)
library(tidyverse); library(patchwork); library(coda); library(rjags)

library(spdep)
library(spatialreg)
```

```{r}
a <- fread("/Users/nani13143553/Desktop/beta.csv")[,-1]
b <- fread("/Users/nani13143553/Desktop/u.csv")[,-1]

suicide254 = fread("/Users/nani13143553/Desktop/KAIST/Research/before CCO/suicide count before standardization 254.csv")[,-c(1)]
death254_nostand = fread("/Users/nani13143553/Desktop/KAIST/Research/before CCO/death count before standardization 254.csv")[,-c(1)]
death254 = fread("/Users/nani13143553/Desktop/KAIST/Research/before CCO/death count before cco 254.csv")[,-c(1)]
death254_nostand$city_num <- as.numeric(gsub("_", "", death254_nostand$city_dist))
death254_nostand$region <- as.factor(death254_nostand$city_num)
death254_nostand$region <- as.numeric(death254_nostand$region)

death1 = subset(death254_nostand, region < 52)
death2 = subset(death254_nostand, 52 < region)
death3 = subset(death2, region < 88)
death4 = subset(death254_nostand, region > 90)
death5 = subset(death4, region < 230)
death6 = subset(death254_nostand, region > 230)
death7 = subset(death6, region < 253)
death247 = rbind(death1, death3, death5, death7)
##### Setting Strata for New Data ##############################################
death247 = death247 %>% arrange(city_dist)
death247$region = as.factor(death247$city_dist)
death247$region = as.numeric(death247$region)

year = year(death247$death_date)
month = str_pad(month(death247$death_date), 2, pad = "0")
dow = wday(death247$death_date)
md = paste(month, dow, sep = "_")
ymd = paste(year, month, dow, sep = "_")
death247$ymd = as.factor(ymd)
death247$zone_ymd = as.factor(paste(str_pad(death247$region, 2, pad = "0"), ymd, sep = "_"))
death247$year = year
death247$temp_mv2 = (death247$temp_lag0 + death247$temp_lag1)/2
death15 = subset(death247, year == 2015)
death15_1 = subset(death15, region < 57)
death15_2 = subset(death15, region > 57)
death15_3 = subset(death15_2, region < 86)
death15_4 = subset(death15, region > 86)
death25415 = rbind(death15_1, death15_3, death15_4)

death16 = subset(death247, year == 2016)
death16_1 = subset(death16, region < 57)
death16_2 = subset(death16, region > 57)
death25416 = rbind(death16_1, death16_2)

death17 = subset(death247, year == 2017)
death17_1 = subset(death17, region < 57)
death17_2 = subset(death17, region > 57)
death25417 = rbind(death17_1, death17_2)

death_res = subset(death247, year > 2017)
death247_1 = rbind(death25415, death25416, death25417, death_res)
################################################################################
##### Finalizing the Data Processing ###########################################
death247_1 = death247_1 %>% arrange(city_dist)
death247_1$region = as.factor(death247_1$city_dist)
death247_1$region = as.numeric(death247_1$region)
death247_1 = death247_1 %>% mutate(temp_mv2 = ifelse(is.na(temp_mv2), temp_lag0, temp_mv2))

year = year(death247_1$death_date)
month = str_pad(month(death247_1$death_date), 2, pad = "0")
dow = wday(death247_1$death_date)
ymd = paste(year, month, dow, sep = "_")
death247_1$ymd = as.factor(ymd)
death247_1$zone_ymd = as.factor(paste(str_pad(death247_1$region, 2, pad = "0"), ymd, sep = "_"))
death247_1$year = year

death247_1$region = as.factor(death247_1$city_dist)
death247_1$region = as.numeric(death247_1$region)

death247_FINAL = death247_1[,-c(5:15, 18:46)]
death247_FINAL$temp_stand = (death247_FINAL$temp_mv2 - mean(death247_FINAL$temp_mv2))/sd(death247_FINAL$temp_mv2)
death247_FINAL$O3_8H_lag0 = death247_FINAL$O3_8H_lag0/10

data_2021 <- subset(death247_FINAL, year == 2021)
```

```{r, warning=FALSE}
beta1 = 0.5; beta2 = -0.01

adjmatrix <- data.frame(fread("/Users/nani13143553/Desktop/Final2 Simulation/adj mat bi (1).csv", header = TRUE)[,-1])
adjmatrix <- adjmatrix[c(1:25, 50:59, 76:117), c(1:25, 50:59, 76:117)]
W <- as.matrix(adjmatrix)
lw <- mat2listw(W, style = "W")
n <- nrow(W)
rho <- 0.8
sigma2 <- 0.1
set.seed(1)
eps <- rnorm(n, sd = sqrt(sigma2))
u1i_sample  <- as.numeric(invIrW(lw, rho) %*% eps)

set.seed(1)
u1i_sample1 = rnorm(247, mean = 0, sd = 2)
u1i_sample = fread("/Users/nani13143553/Desktop/u.csv")[,-1]
u1i_sample = as.vector(colMeans(u1i_sample))

u1i_sample = as.vector(colMeans(b))

death77_1921 <- subset(death77_FINAL, year > 2018)
death77_1921$temp_stand <- (death77_1921$temp_mv2 - mean(death77_1921$temp_mv2))/sd(death77_1921$temp_mv2)

simulation_data4_emp = death77_1921
simulation_data4_emp[,3] = NA
simulation_data4_emp = simulation_data4_emp %>% arrange(region)
simulation_data4_emp$num <- as.numeric(gsub("_", "", simulation_data4_emp$zone_ymd))
simulation_data4_emp$strata <- as.factor(simulation_data4_emp$num)
simulation_data4_emp$strata <- as.numeric(simulation_data4_emp$strata)
simulation_data4_emp = simulation_data4_emp %>% arrange(strata)

death77_1921$strata <- simulation_data4_emp$strata

kappa = rep(NA, length(unique(death77_1921$zone_ymd)))
for (i in 1:length(unique(death77_1921$strata))) {
  kappa[i] <- mean(subset(death77_1921, strata == i)$count)
}

for (z in 1:20) {
set.seed(z)
region_mat = data.frame(as.matrix(array(NA, dim = c(252, 5))))
matrix_list_77 <- vector("list", 77)
for (i in 1:77) {
  matrix_list_77[[i]] <- region_mat
}
for (k in 1:77) {
  regiondat <- subset(simulation_data4_emp, region == k)
  for (i in 1:length(unique(regiondat$strata))) {
  subdata <- subset(regiondat, strata == i + (k-1)*252)
  kappaij <- log(kappa[i + (k-1)*252])
  for (j in 1:nrow(subdata)) {
    matrix_list_77[[k]][i,j] <- exp(kappaij + (beta1 + u1i_sample[k])*subdata[j, 4] + beta2*subdata[j, 13])
  }
 }
}
region_mat = data.frame(as.matrix(array(NA, dim = c(252, 5))))
Y_list_77 <- vector("list", 77)
for (i in 1:77) {
  Y_list_77[[i]] <- region_mat
}
for (i in 1:77) {
  for (j in 1:252) {
    for (k in 1:5) {
     Y_list_77[[i]][j, k] <- rpois(1, matrix_list_77[[i]][j, k])
    }
  }
}
vector <- unlist(lapply(Y_list_77, function(mat) c(t(mat))))
vector_clean <- vector[!is.na(vector)]

simulation_data = simulation_data4_emp
simulation_data[,3] = NA
simulation_data = simulation_data %>% arrange(strata)
simulation_data$count <- vector_clean

write.csv(simulation_data, paste0("/Users/nani13143553/Desktop/Final 77/Spatial Scenario/Large/Data/Data", z, ".csv"))
}
```
