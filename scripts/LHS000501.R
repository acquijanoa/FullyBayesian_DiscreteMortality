library(dplyr)
library(ggplot2)
library(haven)

dir("data/DTA")

data <- haven::read_dta("data/DTA/CNPV2018_5PER_A2_11.DTA")

head(data) %>% View()

colnames(data)

table(data$p_nivel_anosr)
table(data$p_est_civil)
