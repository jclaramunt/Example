predictInnerNett<-function(object,newdata){
  optmax <- object$control$optmax
  nsplit<-dim(object$si)[1]
  ytx <- newdata[,c(names(object$data)[1],names(object$data)[2],object$si[1:nsplit,3])]
  y<-as.matrix(ytx[,1])
  tr<-as.numeric(as.factor(ytx[,2]))-1
  root<-rep(1,dim(ytx)[1]);
  if(!is.factor(ytx[,object$si[1,3]])){
    Gmat <- makeGchmat(gm=root, varx=ytx[,object$si[1,3]], splitpoint=object$si[1,4])  }
  if(is.factor(ytx[,object$si[1,3]])){
    possibleSplits <- determineSplits(x=ytx[,object$si[1,3]], gm=root)
    assigMatrix <- makeCatmat(x=ytx[,object$si[1,3]], gm=root, z=possibleSplits[[1]], splits=possibleSplits[[2]])
    Gmat <- makeGchmatcat(root, object$si[1,4], assigMatrix)  }
  nnum<-c(2,3)
  if (nsplit>1){
    for (i in 2:nsplit){
      o<-which(nnum==object$si[i,1])
      if(!is.factor(ytx[,object$si[i,3]])){
        Gmatch <- makeGchmat(gm=Gmat[,o], varx=ytx[,object$si[i,3]], splitpoint=object$si[i,4])  }
      if(is.factor(ytx[,object$si[i,3]])){
        possibleSplits <- determineSplits(x=ytx[,object$si[i,3]], gm=Gmat[,o])
        assigMatrix <- makeCatmat(x=ytx[,object$si[i,3]], gm=Gmat[,o], z=possibleSplits[[1]], splits=possibleSplits[[2]])
        Gmatch <- makeGchmatcat(gm=Gmat[,o], splitpoint=object$si[i,4], assigMatrix=assigMatrix)  }

      nnum<-c(nnum[-o],object$si[i,1]*2,object$si[i,1]*2+1)
      Gmat<-cbind(Gmat[,-o],Gmatch)
    }}
  End<-cpmat(Gmat,y,tr,optmax=optmax)
  obj<-cbind(nnum,End)
  obj<-as.data.frame(obj)
  names(obj)=c("Leaf","#(T=0)","#(T=1)","d","T*","Yt*")

  predTrt <- (Gmat %*% object$li[order(object$li$node),"T*"]) #-1 #Gmat does not give leaf or node order so first order its columns according to nodes
  return(list(Gmat=Gmat, obj=obj, predTrt=predTrt))
}
