MySI <- function(x,a)
{
  library(flexclust)
  library(cluster)
  # x is the data , one row is one sample, one columns is one feature
  #a is a vector that contains the cluster label of each sample
  #return silhouette to evaluate the cluster
  #return(mean(silhouette(a,dist(x))[,3]))
	#return(mean(silhouette(a,dist2(x,x))[,3]))
	return(mean(silhouette(a,x)[,3]))
  }