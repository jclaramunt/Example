cSubtree <- function(object, internalnodes,nind){
  optmax <- object$control$optmax
  maxd <- object$control$maxd
  all <- as.vector(colnames(object$all))   #all nodes
  nsplit<-dim(object$si)[1]
  ytx <- object$data[,c(names(object$data)[1],names(object$data)[2],object$si[1:nsplit,3])]
  y<-as.matrix(ytx[,1])
  tr<-as.numeric(as.factor(ytx[,2]))-1
  if(optmax) {cmaxroot <- max(mean(y[tr==0]),mean(y[tr==1]))} else {cmaxroot <- min(mean(y[tr==0]),mean(y[tr==1]))} #C in root node
  Csub <- matrix(c(1, cmaxroot, 1),ncol=3)
  colnames(Csub) <- c("t","C(T-t)", "#L(T-t)")
  for (t in as.numeric(internalnodes)){
    childnodes <- cChild(t,maxd)
    TreeWithoutBranch <- object$all[,all %in% c(t, nind[!nind%in%childnodes])]
    Csub <- rbind(Csub, c(t,computeC(cpmat(TreeWithoutBranch, y, tr, optmax)), ncol(TreeWithoutBranch)))
  }
  return(Csub)
}

cSequence <- function(object){
  nsplit<-dim(object$si)[1]
  ytx <- object$data[,c(names(object$data)[1],names(object$data)[2],object$si[1:nsplit,3])]
  y<-as.matrix(ytx[,1])
  tr<-as.numeric(as.factor(ytx[,2]))-1

  nind <- as.vector(colnames(object$nind)) #leaves
  all <- as.vector(colnames(object$all))   #all nodes
  internalnodes <- as.vector(all[!all %in%nind]) #internal nodes
  optmax <- object$control$optmax
  maxd <- object$control$maxd

  CL.t <- cSubtree(object,internalnodes,nind)    #starting matrix, update later
  C.T<- computeC(cpmat(object$nind, y, tr, optmax))
  L.T <- ncol(object$nind)
  alpha <-0
  sequence <- 0
  critvalue <- C.T
  while(!1%in%sequence){
    G.t <- (C.T - CL.t[,2]) / (L.T-CL.t[,3]) #winst die een branch geeft / aantal extra leafs
    alpha <- c(alpha,min(na.omit(G.t)))      #smallest value of split-complexity measure
    cutbranch <- as.numeric(CL.t[,"t"][which.min(G.t)]) #branch that has smallest split-complexity
    sequence <- c(sequence, cutbranch)  #sequence of to be pruned branches (sequence of pruned subtrees)

    C.T <- CL.t[CL.t[,"t"]==cutbranch,2]  #C of new best tree (without branch)
    L.T <-  CL.t[CL.t[,"t"]==cutbranch,3] #L of new best tree (without branch)

    critvalue <- c(critvalue,C.T)
    internalnodes <- as.vector(internalnodes[!internalnodes %in% cChild(cutbranch,maxd)]) #new internalnodes
    nind <- c(cutbranch, nind[!nind%in%cChild(cutbranch,maxd)]) #new leaves

    CL.t <- cSubtree(object,internalnodes,nind) #CL matrix of new best tree (without branch)
  }
  alpha[alpha<0] <- 0
  return(list(sequence, alpha, critvalue))
}
