ctc<-function(pmat,parvec){
  #check treatment cardinality condition; parvec contains the cardinality restrictions of resp. a1 and a2
  #a1 cardinality T=0 in a node ;   #a2 cardinality T=0 in a node
  cond <- ifelseC(rowSums(pmat[,1:2]>=parvec)==2,rep(1, nrow(pmat)),numeric(nrow(pmat)))
  condvec<- as.numeric(sum(cond)==nrow(pmat))
  #if condvec==1 then the cardinality conditions for each node after the split are met
  return(condvec)}

computeC<-function(pmat){
  #compute criterion
  #dmats3: designmatrix with admissible assignments of the nodes to the partition classes
  #compute value of partitioning criterion
  weight<-pmat[,1]+pmat[,2]
  crittot <- (pmat[,5]%*%weight)/sum(weight) #sumoverKleafs(max(y0,y1)*Nk)/N
  return(crittot=crittot)}

bos<-function(y,x,tr,gm,parents,parvec,nsplit,optmax){
  #best observed  split
  #computes best observed split for a splitting variable candidate
  #ouput is result (rowvector of length 5) with: optimal split point, max value of C
  #y = outcome variable
  #x = predictor to be used for the split
  #tr = treatment variable(1=treatment A; 2=treatment B)
  #gm = indicator vector of persons in the particular parent node that is split
  #parents contains information of rest of the parent nodes (cardinalities t0, cardinalities t1, and sum y|t=0, sum y|t=1)
  #parvec=parameter a1 and a2;
  #nsplit: number of split
  if(is.factor(x)){
    z <- unique(sort(x[gm==1]))
    n <- ((2^length(z))-2)/2

    #n is  the number of distinct values (split points) of predictor x
    crittot<-matrix(0,nrow=1,ncol=n)
    #crittot is a matrix with number of rows is number of possible assignments (rows of Ds)
    #and number of columns is n = number of splitpoints; the cells contain the value of the criterion C

    #Perform for each possible split point vik:
    if(n>0){
      possibleSplits <- determineSplits(x, gm)
      assignMatrix <- makeCatmat(x=x, gm=gm, z=possibleSplits[[1]], splits=possibleSplits[[2]])
      #For every possible split determine whether person goes to left node(=1) or right node(=0)
      for(i in 1:n){
        Gch<- makeGchmatcat(gm,i,assignMatrix)
        child<-cpmat(Gch,y,tr,optmax=optmax)
        if(nsplit==1){End<-child}
        if(nsplit>1){End<-as.matrix(rbind(parents,child) )}
        dimnames(End)<-NULL
        check1<-ctc(End,parvec)
        if(check1==1){
          cdat<-computeC(child)
          crittot[,i]<-cdat
        }
      }

      if(optmax) {
        ccmax<-apply(crittot,2,max)
        colstar<-which(ccmax==max(crittot))[1]
        rcmax<-apply(crittot,1,max)
        rowstar<- which(rcmax==max(crittot))[1]
      } else {
        ccmax<-apply(crittot,2,min)
        colstar<-which(ccmax==min(crittot))[1]
        rcmax<-apply(crittot,1,min)
        rowstar<- which(rcmax==min(crittot))[1]
      }

      splitpoint <- colstar
      result<-c(splitpoint,rcmax)}

    if (n==0){result<-numeric(2)}
  }

  #Continuous splitting variables
  if(!is.factor(x)){
    z <- sort(x[gm==1])
    if(length(z)>((4*min(parvec))+2)){
      z <- unique(z[(2*min(parvec)):(length(z)-((2*min(parvec))+1))])
    } else {
      z<- unique(z)
    }
    n<-length(z)

    #n is  the number of distinct values (split points) of predictor x
    crittot<-matrix(0,nrow=1,ncol=n)
    #crittot is a matrix with number of rows is number of possible assignments (rows of Ds)
    #and number of columns is n = number of splitpoints; the cells contain the value of the criterion C
    #Perform for each possible split point vik:
    if(n>1){
      for(i in 1:(n-1)){
        #create indicator matrix of child nodes (Gch) after split on z[i]:
        splitpoint<-(z[i]+z[i+1])/2
        Gch <- makeGchmat(gm,varx=x,splitpoint)

        #child will be filled with information of the two childnodes
        #(cardinalities t0, cardinalities t1, meandifT0-T1, best.t, Yt*)
        child<-cpmat(Gmat=Gch,y,tr,optmax=optmax)
        if(nsplit==1){End<-child}
        if(nsplit>1){End<-as.matrix(rbind(parents,child) )}
        dimnames(End)<-NULL
        #End is matrix E that contains the information of all the end nodes after a split (cardinalities t0, cardinalities t1, Yt0-Yt1, best tr, and mean y of best tr)
        ##check conditions and if these are met compute the value of C which will be collected in matrix crittot
        check1<-ctc(pmat=End,parvec)
        if(check1==1){
          cdat<-computeC(child)
          crittot[,i]<-cdat
        }
      }
    }

    if(optmax) {
      ccmax<-apply(crittot,2,max)
      #ccmax is vector with maximum value of C for each split point
      colstar<-which(ccmax==max(crittot))[1]
      #colstar is the columnnumber referring to the splitpoint on variable Xk that results in the highest value of C
      rcmax<-apply(crittot,1,max)
      #rcmax is vector with maximum value of C
      rowstar<- which(rcmax==max(crittot))[1]
    } else {
      ccmax<-apply(crittot,2,min)
      colstar<-which(ccmax==min(crittot))[1]
      rcmax<-apply(crittot,1,min)
      rowstar<- which(rcmax==min(crittot))[1]  }

    #rowstar is the rownumber referring to the row of dmats (Ds) that results in the highest value of C for this particular predictor variable
    splitpoint<- (z[colstar]+z[colstar+1] )/2
    result<-c(splitpoint,rcmax)

    if (n==1){result<-numeric(2)}
  }
  return(result)
}

bos1 <-function(y,x,tr,gm,parents,parvec,nsplit,optmax,kk){
  #possible  splits for root node
  #ouput is result (rowvector of length 5) with: splitting variable, split point, value of C
  #y = outcome variable
  #x = predictor to be used for the split
  #tr = treatment variable(1=treatment A; 2=treatment B)
  #gm = indicator vector of persons in the particular parent node that is split
  #parents contains information of rest of the parent nodes (cardinalities t0, cardinalities t1, and sum y|t=0, sum y|t=1)
  #parvec=parameter a1 and a2;
  #nsplit: number of split
  if(is.factor(x)){
    z <- unique(sort(x[gm==1]))
    n <- ((2^length(z))-2)/2

    #n is  the number of distinct values (split points) of predictor x
    crittot<-matrix(0,nrow=1,ncol=n)
    #crittot is a matrix with number of rows is number of possible assignments (rows of Ds)
    #and number of columns is n = number of splitpoints; the cells contain the value of the criterion C

    #Perform for each possible split point vik:
    if(n>0){
      possibleSplits <- determineSplits(x, gm)
      assignMatrix <- makeCatmat(x=x, gm=gm, z=possibleSplits[[1]], splits=possibleSplits[[2]])
      #For every possible split determine whether person goes to left node(=1) or right node(=0)
      for(i in 1:n){
        Gch<- makeGchmatcat(gm,i,assignMatrix)
        child<-cpmat(Gch,y,tr,optmax=optmax)
        if(nsplit==1){End<-child}
        if(nsplit>1){End<-as.matrix(rbind(parents,child) )}
        dimnames(End)<-NULL
        check1<-ctc(End,parvec)
        if(check1==1){
          cdat<-computeC(child)
          crittot[,i]<-cdat
        }
      }

      result <- cbind(kk,z,c(crittot))
    }
    if (n==0){result<-c(kk,numeric(2))}
  }

  #Continuous splitting variable
  if(!is.factor(x)){
    z<- sort(x[gm==1])

    if(length(z)>((4*min(parvec))+2)){
      z <- unique(z[(2*min(parvec)):(length(z)-((2*min(parvec))+1))])
    } else {
      z<- unique(z)
    }
    n<-length(z)

    #  x <- as.matrix(x, ncol=1)
    #  x[order(x[,1]),1]

    #n is  the number of distinct values (split points) of predictor x
    crittot<-matrix(0,nrow=1,ncol=n)
    #crittot is a matrix with number of rows is number of possible assignments (rows of Ds)
    #and number of columns is n = number of splitpoints; the cells contain the value of the criterion C
    #Perform for each possible split point vik:
    if(n>1){
      for(i in 1:(n-1)){
        #create indicator matrix of child nodes (Gch) after split on z[i]:
        splitpoint<-(z[i]+z[i+1])/2
        Gch <- makeGchmat(gm,varx=x,splitpoint)
        #child will be filled with information of the two childnodes
        #(cardinalities t0, cardinalities t1, meandifT0-T1, best.t, Yt*)
        child<- cpmat(Gmat=Gch,y,tr,optmax=optmax)

        if(nsplit==1){End<-child}
        if(nsplit>1){End<-as.matrix(rbind(parents,child) )}
        dimnames(End)<-NULL
        #End is matrix E that contains the information of all the end nodes after a split (cardinalities t0, cardinalities t1, Yt0-Yt1, best tr, and mean y of best tr)
        ##check conditions and if these are met compute the value of C which will be collected in matrix crittot
        check1<-ctc(pmat=End,parvec)
        if(check1==1){
          cdat<-computeC(child)
          crittot[,i]<-cdat
        }
      }
      result <- cbind(kk,z,c(crittot))
    }

    if (n==1){result<- matrix(c(kk,numeric(2)), ncol=3)}
  }
  result <- result[result[,3]!=0,]
  return(result)
}

bos23<-function(y,x,tr,gm,parents,parvec,nsplit,optmax){
  #best observed  split of nodes 2 and 3 in deciles
  #computes best observed split for a splitting variable candidate
  #ouput is result (rowvector of length 5) with: optimal split point, max value of C after 3 splits
  #y = outcome variable
  #x = predictor to be used for the split
  #tr = treatment variable(1=treatment A; 2=treatment B)
  #gm = indicator vector of persons in the particular parent node that is split
  #parents contains information of rest of the parent nodes (cardinalities t0, cardinalities t1, and sum y|t=0, sum y|t=1)
  #parvec=parameter a1 and a2;
  #nsplit: number of split
  if(is.factor(x)){
    z <- unique(sort(x[gm==1]))
    n <- ((2^length(z))-2)/2

    #n is  the number of distinct values (split points) of predictor x
    crittot<-matrix(0,nrow=1,ncol=n)
    #crittot is a matrix with number of rows is number of possible assignments (rows of Ds)
    #and number of columns is n = number of splitpoints; the cells contain the value of the criterion C

    #Perform for each possible split point vik:
    if(n>0){
      possibleSplits <- determineSplits(x, gm)
      assignMatrix <- makeCatmat(x=x, gm=gm, z=possibleSplits[[1]], splits=possibleSplits[[2]])
      #For every possible split determine whether person goes to left node(=1) or right node(=0)
      for(i in 1:n){
        Gch<- makeGchmatcat(gm,i,assignMatrix)
        child<-cpmat(Gch,y,tr,optmax=optmax)
        End<-as.matrix(rbind(parents,child) )
        dimnames(End)<-NULL
        check1<-ctc(End,parvec)
        if(check1==1){
          cdat<-computeC(child)
          crittot[,i]<-cdat
        }
      }

      if(optmax) {
        ccmax<-apply(crittot,2,max)
        colstar<-which(ccmax==max(crittot))[1]
        rcmax<-apply(crittot,1,max)
        rowstar<- which(rcmax==max(crittot))[1]
      } else {
        ccmax<-apply(crittot,2,min)
        colstar<-which(ccmax==min(crittot))[1]
        rcmax<-apply(crittot,1,min)
        rowstar<- which(rcmax==min(crittot))[1]
      }

      splitpoint <- colstar
      result<-c(splitpoint,rcmax)
    }
    if (n==0){result<-numeric(2)}
  }

  if(!is.factor(x)){
    z.orig <- x[gm==1]
    n.orig <- length(z.orig)
    if(n.orig> 4*min(parvec) & n.orig<= 20*min(parvec)){
      fl <- ceiling(n.orig/(2*min(parvec)))
      if(fl <= 2){z <- as.vector(quantile(z.orig, prob = .5, type = 5))} else{
        probb <- seq(.5-floor((fl/2))*.1, .5+floor((fl/2))*.1, by=.1)
        z <- as.vector(quantile(z.orig, prob = probb, type = 5)) }#deciles
    }
    if(n.orig> 20*min(parvec)) { z <- as.vector(quantile(z.orig, prob = seq(.1, .9, by=.1), type = 5)) }
    if(n.orig<= 4*min(parvec)) {z <- as.vector(quantile(z.orig, prob = .5, type = 5))}

    n<-length(z)

    #n is  the number of distinct values (split points) of predictor x
    crittot<-matrix(0,nrow=n,ncol=2)
    #crittot is a matrix with number of rows is number of possible assignments (rows of Ds)
    #and number of columns is n = number of splitpoints; the cells contain the value of the criterion C
    #Perform for each possible split point vik:
    if(n>1){
      for(i in 1:n){
        #create indicator matrix of child nodes (Gch) after split on z[i]:
        splitpoint<- z[i]
        Gch <- makeGchmat(gm,varx=x,splitpoint)

        #child will be filled with information of the two childnodes
        #(cardinalities t0, cardinalities t1, meandifT0-T1, best.t, Yt*)
        child<- cpmat(Gmat=Gch,y,tr,optmax=optmax)
        End<-as.matrix(rbind(parents,child) )

        dimnames(End)<-NULL
        #End is matrix E that contains the information of all the end nodes after a split (cardinalities t0, cardinalities t1, Yt0-Yt1, best tr, and mean y of best tr)
        ##check conditions and if these are met compute the value of C which will be collected in matrix crittot
        check1<-ctc(pmat=End,parvec)
        if(check1==1){
          cdat<-computeC(child)
          crittot[i,]<-c(splitpoint,cdat)
        }
      }

      if(optmax) {
        result <- crittot[which.max(crittot[,2]),]
      } else {
        result <- crittot[which.min(crittot[,2]),]
      }
    } else {result<-numeric(2)}
  }
  return(result)
}
