firsttwosplits <- function(ds, critmaxdyn, maxd, Xmat, root, cmaxroot,y,tr, optmax,parvec){
  #split left and right child node after first split, based on deciles
  critmax1 <- critmaxdyn[ds,]

  maxl <- 2^maxd
  allresults <- matrix(0, nrow=maxl, ncol=5)#4
  splitpoints <- matrix(0, nrow=maxl, ncol=1)
  #Make the first split
  if(!is.factor(Xmat[,critmax1[1]])){
    Gmat <- makeGchmat(root, Xmat[,critmax1[1]], critmax1[2])}
  if(is.factor(Xmat[,critmax1[1]])){
    possibleSplits <- determineSplits(x=Xmat[,critmax1[1]], gm=root)
    assigMatrix <- makeCatmat(x=Xmat[,critmax1[1]], gm=root, z=possibleSplits[[1]], splits=possibleSplits[[2]])
    Gmat <- makeGchmatcat(gm=root, splitpoint=critmax1[2], assigMatrix=assigMatrix)}

  ##Keep the child node numbers  nnum; #ncol(Gmat) is current number of leaves (=number of candidate parentnodes)=L ; #ncol(Gmat)+1 is total number of leaves after the split  (Lafter)
  nnum <- c(2,3)
  colnames(Gmat) <- c(2,3)
  L <- ncol(Gmat)

  # ##Keep the results (split information, fit information, end node information) after the first split
  # if(critmax1[3]>cmaxroot){
  #   allresults[1,] <- c(1,critmax1,critmax1[3])
  #   #Keep the splitpoints
  #   ifelse(!is.factor(Xmat[,critmax1[1]]), splitpoints[1] <- critmax1[2], splitpoints[1] <- paste(as.vector(unique(sort(Xmat[Gmat[,1]==1, critmax1[1]]))), collapse=", "))
  #   cmax <- allresults[1, 4]
  #   endinf <- ctmat(Gmat,y,tr)
  # } else { ##if there is no optimal triplet for the first split:
  #   Gmat <- Gmat*0
  #   endinf <- matrix(0, ncol=8, nrow=2)
  # }

  ##Keep the results (split information, fit information, end node information) after the first split
  if(critmax1[3]>cmaxroot){
    allresults[1,] <- c(1,critmax1,critmax1[3])
  } else { ##if there is no optimal triplet for the first split:
    allresults[1,] <- c(1,critmax1,cmaxroot)
  }
  ifelse(!is.factor(Xmat[,critmax1[1]]), splitpoints[1] <- critmax1[2], splitpoints[1] <- paste(as.vector(unique(sort(Xmat[Gmat[,1]==1, critmax1[1]]))), collapse=", "))
  cmax <- allresults[1, 5]
  endinf <- ctmat(Gmat,y,tr)

  stopc <- 0
  Gmatall <- Gmat
  queue <- c(2,3)

  #Repeat the tree growing procedure
  for(i in 2:3){
    Lafter <- ncol(Gmat)+1

    #make parentnode information matrix, select best observed parent node (with optimal triplet)
    parent <- cpmat(Gmat,y,tr,optmax=optmax)
    rownames(parent) <- nnum

    critmax <- bovar23(y, Xmat, tr, gm=Gmat[,paste(i)], parents=parent[-i,], parvec, nsplit=L, optmax=optmax)

    ##Perform the best split and keep results
    if(!is.factor(Xmat[,critmax[1]])){
      Gmatch <- makeGchmat(Gmat[,paste(i)], Xmat[,critmax[1]], critmax[2])  }
    if(is.factor(Xmat[,critmax[1]])){
      possibleSplits <- determineSplits(x=Xmat[,critmax[1]], gm=Gmat[,paste(i)])
      assigMatrix <- makeCatmat(x=Xmat[,critmax[1]], gm=Gmat[,paste(i)], z=possibleSplits[[1]], splits=possibleSplits[[2]])  ###paste(i) instead of "i"
      Gmatch <- makeGchmatcat(gm=Gmat[,paste(i)], splitpoint=critmax[2], assigMatrix=assigMatrix) }

    colnames(Gmatch) <- c(2*i, 2*i+1)

    Gmatnew <- cbind(Gmat, Gmatch)
    Gmatnew <- Gmatnew[,!colnames(Gmatnew) %in% i]   #delete the splitted node from Gmat

    End <- cpmat(Gmatnew,y,tr,optmax=optmax)

    check1<-ctc(pmat=End,parvec)
    if(check1==1){
      cdat<-computeC(End)
      check2<-as.numeric(cdat>cmax)
    } else {
      check2<-0
    }

    if(check1==1&check2==1){
      allresults[L,] <- c(i, critmax[1:3], cdat)
      ifelse(!is.factor(Xmat[,critmax[1]]), splitpoints[L] <- critmax[2], splitpoints[L] <- paste(as.vector(unique(sort(Xmat[Gmatch[,1]==1,critmax[1]]))), collapse=", "))
    } else {
      stopc <- 1
    }

    #update the parameters after the split:
    if(stopc==0) {
      Gmat <- Gmatnew
      Gmatall <- cbind(Gmatall, Gmatch)
      cmax <- allresults[L,5]
      L <- ncol(Gmat)
      nnum <- as.numeric(colnames(Gmat))
      queue <- c(queue[!queue %in% i])#,2*i, 2*i+1)
    } else {
      queue <- queue[!queue %in% i]
      stopc<-0
    }
    #end of for loop
  }

  if(max(allresults[,5]) > cmaxroot){
    result <- allresults[which.max(allresults[,5]),]
  } else {
    result <- c(0,0,0,0,cmaxroot)
  }
  return(result)
}

