MyNMI <- function(a,b)
{
  #a is a vector that contains the cluster label of each sample
  #b is a vector that contains the true label of each sample
  #return NMI to evaluate the cluster
  len_a=length(a)
  xx=matrix(1,len_a,1)
  pa=aggregate(xx,list(a),sum)[,2]/len_a
  pb=aggregate(xx,list(b),sum)[,2]/len_a
  pab=aggregate(xx,list(a,b),sum)[3]/len_a
  
  ha=-sum(pa*log(pa))
  hb=-sum(pb*log(pb))
  hab=-sum(pab*log(pab))
  
  nmi=2*(ha+hb-hab)/(ha+hb)
  return(nmi)
}