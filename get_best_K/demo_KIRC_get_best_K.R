rm(list = ls())
setwd("D:/Desktop/NCICv1.0/")
source("get_best_K/GetSimi.R")
source("get_best_K/MySI.R")
cnv=read.csv("data/KIRC_cnv.csv",header = T,row.names=1)
m=read.csv("data/KIRC_m.csv",header = T,row.names=1)
me=read.csv("data/KIRC_me.csv",header = T,row.names=1)
cnv=as.matrix(cnv)
m=as.matrix(m)
me=as.matrix(me)

simi_cnv_m_me=GetSimi(cnv,m,me)
dis_data=1./simi_cnv_m_me
SI_matrix=matrix(nrow = 7,ncol = 2)
for(j in 2:8)
{
	label_path=paste("res/KIRC_cnv_m_me_",j,".csv",sep="")
	label=read.csv(label_path,header=FALSE)
	SI_matrix[j-1,1]=j;
	SI_matrix[j-1,2]=MySI(dis_data,label[,1])
}


#write out Sihouette value
SI_matrix=as.data.frame(SI_matrix)
write.csv(SI_matrix,"get_best_K/res/KIRC_cnv_m_me_SI.csv",row.names = FALSE)
#plot the SI picture
write_SI_pdf_path=paste("get_best_K/res/KIRC_cnv_m_me_SI.pdf",sep="")
pdf(write_SI_pdf_path)
plot(SI_matrix,type="o",lty=6,lwd=2,pch=17,xlab="K",ylab="SI")
dev.off()