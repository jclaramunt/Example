child <- function(node){
  return(list(node*2,(node*2)+1))
}

cChild <- function(nodenumber, maxd){
  if(nodenumber==0){nodenum <- seq(1:(2^(maxd+1)-1))} else {
    nodenum<-nodenumber
    leaves<-nodenumber
    while(sum(nodenum>=(2^maxd))==0){
      leaves <- sapply(leaves, child)
      nodenum <- unlist(c(nodenum, leaves))
      leaves <- unlist(leaves)
    }}
  return(as.vector(nodenum))
}
