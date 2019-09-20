GetSimi <- function(cnv,m,me)
{
	library(SNFtool) 
	library(flexclust)
	cnv=t(cnv)
	m=t(m)
	me=t(me)
	Dist_cnv=(dist2(as.matrix(cnv),as.matrix(cnv)))^(1/2)
	Dist_m=(dist2(as.matrix(m),as.matrix(m)))^(1/2)
	Dist_me=(dist2(as.matrix(me),as.matrix(me)))^(1/2)
	Wcnv=affinityMatrix(Dist_cnv)
	Wm=affinityMatrix(Dist_m)
	Wme=affinityMatrix(Dist_me)
	return (SNF(list(Wcnv,Wm,Wme)))
 }
#if you have one type data or two type data need calculate similarity
# one type data
#GetSimi<- function(cnv)
#{
#	library(SNFtool) 
#	library(flexclust)
#	cnv=t(cnv)
#	Dist_cnv=(dist2(as.matrix(cnv),as.matrix(cnv)))^(1/2)
#	return (affinityMatrix(Dist_cnv))
#}
#  
#  two type data
#GetSimi<- function(cnv,m)
#{
#	library(SNFtool) 
#	library(flexclust)
#	cnv=t(cnv)
#	m=t(m)
#	Dist_cnv=(dist2(as.matrix(cnv),as.matrix(cnv)))^(1/2)
#	Dist_m=(dist2(as.matrix(m),as.matrix(m)))^(1/2)
#	Wcnv=affinityMatrix(Dist_cnv)
#	Wm=affinityMatrix(Dist_m)
#	return (SNF(list(Wcnv,Wm)))	
#}