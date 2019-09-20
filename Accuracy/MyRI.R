MyRI <- function(a,b)
{
	library(flexclust)
  #a is a vector that contains the cluster label of each sample
  #b is a vector that contains the true label of each sample
  #return RI to evaluate the cluster
  return(flexclust::comPart(a,b,type=c('RI')))
}