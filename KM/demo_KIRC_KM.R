rm(list = ls())
setwd("D:/Desktop/NCICv1.0/")
source("KM/MyKM.R")
#get the clinical data
clinical=read.csv("KM/KIRC_Clinical.csv",header = T)
# get the sample column
sample_col=read.csv("KM/KIRC_sample.csv",header = F)
# get the NCIC method res (sample label)
label=read.csv("res/KIRC_cnv_m_me_2.csv",header = F)
label=cbind(sample_col,label)
KM(label,clinical=clinical)
KM(label,clinical=clinical,fileName="KM/res/KIRC_cnv_m_me.pdf")