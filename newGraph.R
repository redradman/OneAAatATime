library(tidyverse)
library(ggplot2)
library(ggthemes)
library(dplyr)

entireProtein <- read.csv("entireProtein.csv")
proteinshort <- entireProtein[, -c(7)]
proteinshort <- proteinshort[, (-c(11))]

proteinshort 

