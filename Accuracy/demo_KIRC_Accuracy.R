rm(list = ls())
setwd("D:/Desktop/NCICv1.0/")
source("Accuracy/MyARI.R")
source("Accuracy/MyNMI.R")
source("Accuracy/MyRI.R")

# get the golden label and NCIC method result of sample label
Golden=read.csv("Accuracy/Golden/KIRC_Golden_label.csv",header = FALSE)
label=read.csv("res/KIRC_cnv_m_me_2.csv",header = FALSE)

print(paste("NMI:",MyNMI(Golden[,1],label[,1])))
write.table(MyNMI(Golden[,1],label[,1]),file="Accuracy/res/NMI.csv",col.names = FALSE)

print(paste("RI:",MyRI(Golden[,1],label[,1])))
write.table(MyRI(Golden[,1],label[,1]),file="Accuracy/res/RI.csv",col.names = FALSE)

print(paste("RI:",MyARI(Golden[,1],label[,1])))
write.table(MyARI(Golden[,1],label[,1]),file="Accuracy/res/ARI.csv",col.names = FALSE)
