library(cpr)
library(cprtesting)
library(qwraps2)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lme4)
library(knitr) 

knitr::opts_chunk$set(cache = TRUE, collapse = TRUE)
options(tibble.print_max = 5, tibble.print_min = 5) 

rmarkdown::render('slides.Rmd')
