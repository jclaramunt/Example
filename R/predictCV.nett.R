predictCV.nett<-function(object,newdat,cut){
  optmax <- object$control$optmax
  maxd <- object$control$maxd
  nsplit<-dim(object$si)[1]
  ytx.new <- newdat[,c(names(object$data)[1],names(object$data)[2],object$si[1:nsplit,3])]
  y.new<-as.matrix(ytx.new[,1])
  tr.new<-as.numeric(as.factor(ytx.new[,2]))-1
  ytx <- object$data[,c(names(object$data)[1],names(object$data)[2],object$si[1:nsplit,3])]
  y<-as.matrix(ytx[,1])
  tr<-as.numeric(as.factor(ytx[,2]))-1

  if(length(cut)>1) cut <- cut[cut!=0]

  if(cut==0){Gmat.new <- predictInnerNett(object, newdat)$Gmat
  MeanT.col <- predictInnerNett(object, newdat)$predTrt+1
  predY <- tryCatch(ctmat(Gmat.new, y.new, tr.new)[,c(2,5)][, MeanT.col][t(Gmat.new==T)], error=function(e) NA)     }
  if(1%in%cut){if(optmax) {t.max <- which.max(c(mean(y[tr==0]),mean(y[tr==1])))-1} else {tmax <- which.min(c(mean(y[tr==0]),mean(y[tr==1])))-1}
    predY <- rep(mean(y.new[tr.new==t.max]),nrow(newdat))}
  if(cut[1]!=0 & !1%in%cut){
    root<-rep(1,dim(ytx.new)[1])
    if(!is.factor(ytx.new[,object$si[1,3]])){
      Gmat.new <- makeGchmat(gm=root, varx=ytx.new[,object$si[1,3]], splitpoint=object$si[1,4])  }
    if(is.factor(ytx.new[,object$si[1,3]])){
      possibleSplits <- determineSplits(x=ytx.new[,object$si[1,3]], gm=root)
      assigMatrix <- makeCatmat(x=ytx.new[,object$si[1,3]], gm=root, z=possibleSplits[[1]], splits=possibleSplits[[2]])
      Gmat.new <- makeGchmatcat(root, object$si[1,4], assigMatrix)  }
    nnum<-c(2,3)

    cut.all <- unique(unlist(c(sapply(1:length(cut), function(v) cChild(cut[v],maxd)))))

    if(length(object$si[,1][!object$si[,1]%in%cut.all & 1!=object$si[,1]])!=0){
      for (i in which(!object$si[,1]%in%cut.all & 1!=object$si[,1])){
        o<-which(nnum==object$si[i,1])
        if(!is.factor(ytx.new[,object$si[i,3]])){
          Gmatch <- makeGchmat(gm=Gmat.new[,o], varx=ytx.new[,object$si[i,3]], splitpoint=object$si[i,4])  }
        if(is.factor(ytx.new[,object$si[i,3]])){
          possibleSplits <- determineSplits(x=ytx.new[,object$si[i,3]], gm=Gmat.new[,o])
          assigMatrix <- makeCatmat(x=ytx.new[,object$si[i,3]], gm=Gmat.new[,o], z=possibleSplits[[1]], splits=possibleSplits[[2]])
          Gmatch <- makeGchmatcat(gm=Gmat.new[,o], splitpoint=object$si[i,4], assigMatrix=assigMatrix)  }
        nnum<-c(nnum[-o],object$si[i,1]*2,object$si[i,1]*2+1)
        Gmat.new<-cbind(Gmat.new[,-o],Gmatch) } }

    best.t <- cpmat(object$all[,paste(nnum)], y, tr, optmax)[,4] #######
    MeanT.col <- Gmat.new%*%best.t+1
    predY <- ctmat(Gmat.new, y.new, tr.new)[,c(2,5)][, MeanT.col][t(Gmat.new==T)]
  }
  return(predY)
}
