# the input data (clusters and clinical) should not row.names=1
KM <- function(clusters, clinical= NULL, fileName = NULL) 
{
	#add -01 in the clinical data  
	add01 <- function(str) 
	{
    return(paste(str, '-01', sep=''))
	}
	
	clusters <- as.data.frame(clusters)
	colnames(clusters)[1] <- 'name'

	clinical <- as.data.frame(clinical)
	colnames(clinical)[1] <- 'name'
	
	clinical[1] <- apply(clinical[1], 2, add01)
	finData <- merge(clusters, clinical)
	colnames(finData)[2] <- 'clusters'
	colnames(finData)[3] <- 'OS_Days'
	colnames(finData)[4] <- 'OS_Status'
	
	library(survival)
	OS_DAYS<-finData$OS_Days
	OS_DAYS<-as.numeric(OS_DAYS)
  
	OS_STATUS<-finData$OS_Status
	OS_STATUS <- as.character(OS_STATUS)
  
	GROUP<-finData$clusters
	GROUP <- as.numeric(GROUP)
  
	maxGroup <- max(GROUP)
	strGroup <- 'Group1'
	if (maxGroup >= 2) for (i in 2:maxGroup) 
	{
    strGroup <-c(strGroup,paste('Group', i, sep='')) 
	}
	
	corlor=rainbow(maxGroup)
	kmsurviaval<-survfit(Surv(OS_DAYS,OS_STATUS=='Dead') ~GROUP)
	pdf(fileName)
	plot(kmsurviaval,col =corlor,lty=1,lwd=2)
	title(main="",xlab = "Time(Days)",ylab = "Survival Probability")
	legend(("topright"),strGroup,fill= corlor)
	cox = coxph(Surv(OS_DAYS,OS_STATUS=='Dead') ~GROUP)
	legend("top",legend = paste("cox p:", round(summary(cox)$sctest[3],digits = 5),"  ",sep=" "))
	p2=round(summary(cox)$sctest[3],digits = 5)
	dev.off()
}