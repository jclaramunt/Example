#' Non-parametric Estimation of Tree-based Treatment regimes
#'
#' Tree-based approach for estimating optimal treatment regimes in randomized clinical trials that
#' directly maximizes an estimator of the overall expected outcome for the tree-based regimes under study
#' without making assumptions about models underlying the data, and with a look-ahead strategy to overcome greediness.
#'
#' @param formula a description of the model to be fit. The format is \code{Y ~ T | X1 + \dots + XJ},
#'   where the variable before the | represents the dichotomous treatment variable \eqn{T}
#'   and the variables after the | are the baseline characteristics used for partitioning.
#'   If the data are in the order \code{Y, T, X1,\dots, XJ}, no formula is needed.
#'   The lay-out of this formula is based on Zeileis & Croissant (2010).
#' @param data a dataframe containing the variables in the model.
#' @param control a list with control parameters as returned by \code{\link{nett.control}}.
#'
#' @references Doove, L. L., Dusseldorp, E., Van Deun, K., & Van Mechelen, I. (2015).
#'  A novel method for estimating optimal tree-based treatment regimes in randomized clinical trials.
#'   \emph{Manuscript submitted for publication.}
#'
#'   Zeileis A. and Croissant Y. (2010). Extended model formulas in R: Multiple parts and
#'   multiple responses. \emph{Journal of Statistical Software, 34(1)}, 1-13.
#'
#' @seealso \code{\link{summary.nett}}, \code{\link{nett.control}},
#'   \code{\link{prune.nett}}, \code{\link{predict.nett}}
#'
#' @examples
#'# Example with data from the Breast Cancer Recovery Project
#'library(quint)
#'data(bcrp)
#'set.seed(10)
#'head(bcrp)
#'
#'# Experimental condition: education vs. non-education
#'bcrp[,"education"] <- as.numeric(bcrp[,"cond"]==2)
#'
#'# Perform a NETT analysis
#'formula1 <- I(cesdt1 - cesdt3) ~ education | cesdt1 + negsoct1 + uncomt1 +
#'  disopt1 + comorbid + age + wcht1 + nationality + marital + trext
#'control  <- nett.control(a1=10, a2=10, optmax=TRUE, maxd=3)
#'fi <- nett(formula=formula1, data=bcrp, control=control)
#'
#'# Inspect the main results of the analysis
#'summary(fi)
#'
#' @keywords tree
#' @keywords cluster
#'
#' @import quint
#' @import Formula
#' @importFrom stats model.frame na.omit sd terms var quantile
#' @importFrom utils combn
#' @useDynLib NETT
#' @exportClass nett
#' @export
nett<-function(formula, data, control=NULL){
  #dat:data; first column in dataframe = the response variable
  #second column in dataframe = the dichotomous treatment vector
  #rest of the columns in dataframe are the predictors

  tcall <- match.call()
  dat <- as.data.frame(data)
  if (missing(formula)) {
    y <- dat[, 1]
    tr <- dat[, 2]
    Xmat <- dat[, -c(1, 2)]
    dat <- na.omit(dat)
    if (length(levels(as.factor(tr))) != 2) {
      stop("Cannot be performed. The number of treatment conditions does not equal 2.")
    }
  } else {
    F1 <- Formula(formula)
    mf1 <- model.frame(F1, data = dat)
    y <- as.matrix(mf1[, 1])
    origtr <- as.factor(mf1[, 2])
    tr <- as.numeric(origtr)-1
    if (length(levels(origtr)) != 2) {
      stop("Cannot be performed. The number of treatment conditions does not equal 2.")
    }
    Xmat <- mf1[, 3:dim(mf1)[2]]
    dat <- cbind(y, tr, Xmat)
    dat <- na.omit(dat)
    cat("Treatment variable (T) equals 0 corresponds to",
        attr(F1, "rhs")[[1]], "=", levels(origtr)[1], "\n")
    cat("Treatment variable (T) equals 1 corresponds to",
        attr(F1, "rhs")[[1]], "=", levels(origtr)[2], "\n")
    names(dat)[1:2] <- names(mf1)[1:2]
  }
  cat("The sample size in the analysis is", dim(dat)[1], "\n")

  N<-length(y)
  if(is.null(control)) {
    control <- nett.control()  #Use default control parameters and criterion
  }

  #specify parameters a and b  (parvec), maximum number of leaves, maximiziation:
  parvec <- control$parvec
  maxd <- control$maxd
  optmax <- control$optmax

  #if no control argument was specified ,use default parameter values for treatment cardinaltiy condition:
  if(is.null(parvec)){
    a1 <- round(sum(tr==0)/10)
    a2 <- round(sum(tr==1)/10)
    parvec <- c(a1, a2)
    control$parvec <- parvec
  }

  #Start of the tree growing: all persons are in the rootnode, L=1; Criterion value (cmax)=0
  root <- rep(1, length(y))
  meany1 <- mean(y[tr==1])
  meany0 <- mean(y[tr==0])
  if(optmax) {cmaxroot <- max(meany0,meany1)} else {cmaxroot <- min(meany0,meany1)}

  #Step 1
  #Select the optimal triplet for the first split: the split resulting in the maximum value of the criterion after 3 splits
  #use the rootnode information: cardinality t=0, cardinality t=1, meant0-meant1
  rootvec <- c(sum(tr==0), sum(tr==1), meany0-meany1)

  #List of all possible splits for every predictor
  critmax<-lapply(1:dim(Xmat)[2],function(kk,y,Xmat,tr,gm,parents,parvec,nsplit,optmax){
    bos1(y,x=Xmat[,kk],tr,gm=root,parents=rootvec,parvec,nsplit=1,optmax=optmax,kk=kk)},
    Xmat=Xmat,gm=root,y=y,tr=tr,parents=rootvec,nsplit=1,parvec=parvec,optmax=optmax)
  critmaxdyn <- do.call(rbind, critmax)
  #critmaxdyn, matrix of all possible first splits: splitting variable, splitpoint, value of C

  #For every possible first split, split the left and right child node according to deciles
  allresultsdyn <- t(sapply(1:nrow(critmaxdyn), function(ds) firsttwosplits(ds,critmaxdyn, maxd, Xmat, root, cmaxroot,y,tr,optmax, parvec)))

  #allresultsdyn: splitnumber, splitting variable, splitpoint, value of C of specific split, total value of C
  critmax1 <- critmaxdyn[which.max(allresultsdyn[,5]),]
  #critmax1, first split that results in highest C-value after 3 splits: splitting variable, splitpoint, value of C after first split

  #Step 2
  #Create matrix for results
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
  #   stopc <-0
  # } else { ##if there is no optimal triplet for the first split:
  #   Gmat <- Gmat*0
  #   Gmatall <- Gmat
  #   endinf <- matrix(0, ncol=8, nrow=2)
  #   stopc <- 1
  # }

  if(max(allresultsdyn[,5])>cmaxroot){
    if(critmax1[3]>cmaxroot){allresults[1,]<-c(1,critmax1,critmax1[3])}else{allresults[1,] <- c(1,critmax1,cmaxroot)}
    ifelse(!is.factor(Xmat[,critmax1[1]]), splitpoints[1] <- critmax1[2], splitpoints[1] <- paste(as.vector(unique(sort(Xmat[Gmat[,1]==1, critmax1[1]]))), collapse=", "))
    cmax <- allresults[1, 5]
    endinf <- ctmat(Gmat,y,tr)
    stopc <-0
  } else {##if there is no optimal after 2nd layer splits:
    Gmat <- Gmat*0
    Gmatall <- Gmat
    endinf <- matrix(0, ncol=8, nrow=2)
    stopc <- 1
  }

  #Step 3
  #Repeat the tree growing procedure
  if(stopc==0){
    Gmatall <- Gmat
    queue <- c(2,3)

    while(length(queue)!=0&queue[1]<=((2^maxd)-1)){
      i <- queue[1]
      Lafter <- ncol(Gmat)+1

      #make parentnode information matrix, select best observed parent node (with optimal triplet)
      parent <- cpmat(Gmat,y,tr,optmax=optmax)
      rownames(parent) <- nnum

      critmax <- bovar(y, Xmat, tr, gm=Gmat[,paste(i)], parents=parent[-i,], parvec, nsplit=L, optmax=optmax)
      #critmax: splitting variable, splitpoint, value of C after first split

      if(is.na(critmax[2])==F){
        ##Perform the best split and keep results
        if(!is.factor(Xmat[,critmax[1]])){
          Gmatch <- makeGchmat(Gmat[,paste(i)], Xmat[,critmax[1]], critmax[2]) }
        if(is.factor(Xmat[,critmax[1]])){
          possibleSplits <- determineSplits(x=Xmat[,critmax[1]], gm=Gmat[,paste(i)])
          assigMatrix <- makeCatmat(x=Xmat[,critmax[1]], gm=Gmat[,paste(i)], z=possibleSplits[[1]], splits=possibleSplits[[2]])
          Gmatch <- makeGchmatcat(gm=Gmat[,paste(i)], splitpoint=critmax[2], assigMatrix=assigMatrix) }

        colnames(Gmatch) <- c(2*i, 2*i+1)

        Gmatnew <- cbind(Gmat, Gmatch)
        Gmatnew <- Gmatnew[,!colnames(Gmatnew) %in% i]   #delete the splitted node from Gmat

        End <- cpmat(Gmatnew,y,tr,optmax=optmax)

        #check cardinalities (check 1) and whether new C value is higher than previous C value (check 2)
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
          queue <- c(queue[!queue %in% i],2*i, 2*i+1)
        } else {
          queue <- queue[!queue %in% i]
          stopc<-0
        }
      }else{
        queue <- queue[!queue %in% i]
        stopc<-0
      }
      #end of while loop
    }}

  Lfinal <- ncol(Gmat)  #Lfinal=final number of leaves of the tree

  #create endnode information of the tree
  endinf <- matrix(0,nrow=length(nnum),ncol=10)
  if(sum(Gmat[,1]!=0)){
    endinf[,c(2:9)] <- ctmat(Gmat,y,tr)}
  endinf <- data.frame(endinf)
  if(optmax) {endinf[,10] <- as.numeric(endinf[,3]<endinf[,6])} else {endinf[,10] <- as.numeric(endinf[,3]>=endinf[,6])}

  endinf[,1] <- nnum
  index <- leafnum(nnum)
  endinf <- endinf[index,]
  rownames(endinf) <- paste("Leaf ",1:Lfinal,sep="")
  colnames(endinf) <- c("node","#(T=0)", "meanY|T=0", "SD|T=0","#(T=1)", "meanY|T=1","SD|T=1","d", "se", "T*")

  if(Lfinal==2){allresults <- data.frame(t(c(2,allresults[1,])))}
  if(Lfinal>2){allresults <- data.frame(cbind(2:Lfinal, allresults[1:(Lfinal-1),]))}

  allresults[,3] <- colnames(Xmat)[allresults[,3]]
  splitnr <- 1:(Lfinal-1)
  allresults <- cbind(splitnr, allresults)
  colnames(allresults) <- c("split", "#leaves", "parentnode", "splittingvar", "splitpoint", "split statistic","C")

  si <- allresults[,3:6]
  cn <- paste(si[,1]*2, si[,1]*2+1, sep=",")
  si <- cbind(parentnode=si[,1], childnodes=cn, si[,2:3])
  rownames(si) <- paste("Split ", 1:(Lfinal-1), sep="")
  object <- list(call=tcall, control=control, data=dat, si=si, fi=allresults[,c(1:2,6:7)], li=endinf, nind=Gmat[,index], all=Gmatall)
  class(object) <- "nett"

  return(object)
}
