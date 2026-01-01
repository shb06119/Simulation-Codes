# Data Application
We apply three modeling approahces (two-stage method, Bayesian hierarchical model, and spatial Bayesian hierarchical model) to empirical dataset. The dataset includes daily death counts, ozone concentration (8 hour moving average), and temperature (Kelvin) in 2015-2021. 
```{r, warning=FALSE, message=FALSE}
library(sf); library(gnm); library(data.table); library(lubridate); library(dplyr)
library(ggplot2); library(spdep); library(fields); library(leaflet)  
library(metafor); library(R2jags); library(parallel); library(devtools)
library(tidyverse); library(patchwork); library(coda); library(rjags)

library(spdep); library(spatialreg)
```
