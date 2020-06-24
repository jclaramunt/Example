makeGchmat<-function(gm,varx,splitpoint){
  #creates an indicator matrix of the child nodes
  #gm = indicator vector of persons in the parent node that is split
  #varx =  unsorted datavector of splitting predictor X
  #splitpoint is optimal threshold on varx used for splitting
  Gchmat<-matrix(nrow=length(gm),ncol=2)
  Gchmat<-cbind(ifelseC(gm==1&varx<=splitpoint,rep(1,length(gm)),numeric(length(gm))),ifelseC(gm==1&varx>splitpoint,rep(1,length(gm)),numeric(length(gm))))
  return(Gchmat)}

makeCatmat <- function(x,gm,z,splits){
  #create indicator matrix of childnodes for categorical variable
  if(length(z)>1){
    if(length(z)==2) { TrueFalse <- as.numeric(sapply(1:length(x), function(p) x[p]%in%splits))
    assigMatrix <- cbind(TrueFalse, (TrueFalse-1)^2)
    }
    if(length(z)==3) {assigMatrix <- sapply(1:ncol(splits), function(k) sapply(1:length(x), function(p) x[p]%in%splits[,k]))}
    if(length(z)>=4) {assigArray <- sapply(1:length(splits), function(w) sapply(1:ncol(splits[[w]]), function(k) sapply(1:length(x), function(p) x[p]%in%splits[[w]][,k])))
    assigMatrix <- do.call(cbind, assigArray)} #indicator column for every possible splitpoint
  } else {assigMatrix <- matrix(rep(c(1,0)), nrow=length(gm), ncol=2, byrow=T)}
  return(assigMatrix)}

makeGchmatcat<-function(gm,splitpoint,assigMatrix){
  #splits of categorical variable
  Gchmat<-matrix(nrow=length(gm),ncol=2)
  if(splitpoint!=0){Gchmat <- cbind(ifelseC(gm==1&assigMatrix[,splitpoint], rep(1,length(gm)),numeric(length(gm))), ifelseC(gm==1&!assigMatrix[,splitpoint],rep(1,length(gm)),numeric(length(gm))))}
  if(splitpoint==0){Gchmat <- matrix(rep(c(1,0)),nrow=length(gm),ncol=2,byrow=TRUE)}
  return(Gchmat)
}
