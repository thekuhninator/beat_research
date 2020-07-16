#!/usr/bin/env Rscript

# packages for beat
install.packages("devtools", repos = "http://cran.us.r-project.org")
library(devtools)
install_github('theislab/kBET')

install.packages("bapred", repos = "http://cran.us.r-project.org")
install.packages("ggplot2", repos = "http://cran.us.r-project.org" )
install.packages('dplyr', repos = "http://cran.us.r-project.org")
install.packages('gplots', repos = "http://cran.us.r-project.org")
install.packages('getopt', repos = "http://cran.us.r-project.org")
install.packages('sva', repos = "http://cran.us.r-project.org")
install.packages('bapred', repos = "http://cran.us.r-project.org")
install.packages('readr', repos = "http://cran.us.r-project.org")
install.packages('matrixStats', repos = "http://cran.us.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("M3C")
# 
# # packages for multi_beat
install.packages("png", repos = "http://cran.us.r-project.org")
install.packages("raster", repos = "http://cran.us.r-project.org")
install.packages('hash', repos = "http://cran.us.r-project.org")
install.packages("grid", repos = "http://cran.us.r-project.org")
install.packages("gridExtra", repos = "http://cran.us.r-project.org")
install.packages("knitr", repos = "http://cran.us.r-project.org")
install.packages("cowplot", repos = "http://cran.us.r-project.org")