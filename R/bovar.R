bovar<-function(y,Xmat,tr,gm,parents,parvec,nsplit,optmax){
  #selects best observed splitting candidate
  #create matrix to keep the highest value of the criterion
  critmax<-matrix(0,nrow=dim(Xmat)[2],ncol=3)
  critmax[,1]<-1:ncol(Xmat)
  #split only if there are more than 1 persons in the node
  if(sum(gm)>1){
    critmax[,c(2:3)]<-t(sapply(1:dim(Xmat)[2],function(kk,y,Xmat,tr,gm,parents,parvec,nsplit,optmax){
      bos(y,x=Xmat[,kk],tr,gm,parents,parvec,nsplit,optmax=optmax)},
      Xmat=Xmat,gm=gm,y=y,tr=tr,parents=parents,nsplit=nsplit,parvec=parvec,optmax=optmax))
    #which predictor is the best splitting candidate for this split?
    if(optmax){ bestrow<-which(critmax[,3]==max(critmax[,3]))[1]
    } else {bestrow<-which(critmax[,3]==min(critmax[,3]))[1] } }
  else{bestrow<-1}
  return(critmax[bestrow,])
}


bovar23 <-function(y,Xmat,tr,gm,parents,parvec,nsplit,optmax){
  #selects best observed splitting candidate for nodes 2 and 3 (deciles)
  #create matrix to keep the highest value of the criterion
  critmax<-matrix(0,nrow=dim(Xmat)[2],ncol=3)
  critmax[,1]<-1:ncol(Xmat)
  #split only if there are more than 1 persons in the node
  if(sum(gm)>1){
    critmax[,c(2:3)]<-t(sapply(1:dim(Xmat)[2],function(kk,y,Xmat,tr,gm,parents,parvec,nsplit,optmax){
      bos23(y,x=Xmat[,kk],tr,gm,parents,parvec,nsplit,optmax=optmax)},
      Xmat=Xmat,gm=gm,y=y,tr=tr,parents=parents,nsplit=nsplit,parvec=parvec,optmax=optmax))
    #which predictor is the best splitting candidate for this split?
    if(optmax){ bestrow<- which.max(critmax[,3])
    } else {
      bestrow<- which.min(critmax[,3])
    }
  } else {
    bestrow<-1
  }
  return(critmax[bestrow,])
}
