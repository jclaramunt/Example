#' Pruning of a NETT Tree
#'
#' Determines the optimally pruned size of the tree.
#'
#' @param tree fitted tree of the class \code{nett}.
#' @param pp pruning parameter, the constant (\eqn{c}) to be used in the \eqn{c*}standard
#'   error rule. The default value is 1.
#' @param B the number of bootstrap samples.
#' @param \dots optional additional arguments.
#'
#' @return Returns an object of class \code{nett}. The number of leaves of this object is
#'   equal to the optimally pruned size of the tree.
#'
#' @references Breiman L., Friedman J.H., Olshen R.A. and Stone C.J. (1984).
#'   \emph{Classification and Regression Trees}. Chapman & Hall/CRC: Boca Raton.
#'
#'   Doove, L. L., Dusseldorp, E., Van Deun, K., & Van Mechelen, I. (2015).
#'   A novel method for estimating optimal tree-based treatment regimes in randomized clinical trials.
#'   \emph{Manuscript submitted for publication.}
#'
#'   LeBlanc M. and Crowley J. (1993). Survival trees by goodness of split.
#'   \emph{Journal of the American Statistical Association, 88,} 457-467.
#'
#' @seealso \code{\link{nett.control}}, \code{\link{nett}}
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
#'fitnett <- nett(formula=formula1, data=bcrp, control=control)
#'
#'# Prune the tree
#'fit <- prune(tree=fitnett, B=2,pp=1)
#'
#'# Inspect the main results of the analysis
#'summary(fit)
#'
#' @keywords tree
#'
#' @import quint
#' @importFrom rpart prune
#' @importFrom utils combn
#' @useDynLib NETT
#' @exportClass nett
#' @export
prune.nett <- function(tree,pp=1,B=20,...){
  object<-tree
  complParam <- cSequence(object) #create sequence of pruned subtrees
  alpha <- vector()
  for(m in 1:(length(complParam[[2]])-1)){alpha[m] <- sqrt(complParam[[2]][m] * complParam[[2]][m+1])}
  alpha <- c(alpha, 10^4)
  con<-object$control
  optmax <- object$control$optmax
  maxd <- object$control$maxd
  predictedY <- list()
  Cboot <- list()

  emptyleaf<-TRUE
  while(emptyleaf==TRUE){
    indexboot<-Bootstrap(y=object$data[,1],B,tr=object$data[,2])

    for (b in 1:B){
      besttree <- nett(data = object$data[indexboot[,b],], control = con)
      nsplit<-dim(besttree$si)[1]
      if(nsplit>0&!is.na(besttree$si[1,3])){
        complParamB <- cSequence(besttree)

        #In welk interval vallen de originele alpha waardes?
        Binterval <- findInterval(alpha, complParamB[[2]])
        intervalinfo <- cbind(orig.branch=complParam[[1]],orig.alpha=alpha, intervalBtree=Binterval, B.branch=complParamB[[1]][Binterval], B.crit=complParamB[[3]][Binterval])
        uniqueBbranch <- unique(intervalinfo[,"B.branch"])
        predictedY[[b]] <- sapply(1:length(uniqueBbranch), function(p) predictCV.nett(object=besttree, newdat=object$data, cut=uniqueBbranch[1:p]))
        colnames(predictedY[[b]]) <- uniqueBbranch
        predictedY[[b]]  <- predictedY[[b]][,paste(intervalinfo[,"B.branch"])]
        Cboot[[b]] <- complParamB[[3]][Binterval]
      }else{
        yt <- besttree$data[,c(names(besttree$data)[1],names(besttree$data)[2])]
        y<-as.matrix(yt[,1])
        tr<-as.numeric(as.factor(yt[,2]))-1
        y.orig <-as.matrix(object$data[,paste(names(object$data)[1])])
        tr.orig <-as.numeric(as.factor(object$data[,paste(names(object$data)[2])]))-1
        if(optmax) {t.max <- which.max(c(mean(y[tr==0]),mean(y[tr==1])))-1} else {tmax <- which.min(c(mean(y[tr==0]),mean(y[tr==1])))-1}
        predictedY[[b]] <- matrix(mean(y.orig[tr.orig==t.max]),nrow=nrow(object$data),ncol=length(alpha))
        Cboot[[b]] <- mean(y[tr==t.max])
      }
    }
    if(sum(is.na(unlist(predictedY)))==0) emptyleaf<-FALSE
  }

  C.orig <- do.call('cbind',lapply(1:B, function(x) apply(predictedY[[x]],2,mean)))
  C.boot <- do.call('cbind',lapply(1:B, function(x) Cboot[[x]]))

  y.app<-as.matrix(object$data[,names(object$data)[1]])
  tr.app<-as.numeric(as.factor(object$data[,names(object$data)[2]]))-1
  if(optmax) {cmaxroot <- max(mean(y.app[tr.app==0]),mean(y.app[tr.app==1]))} else {cmaxroot <- min(mean(y.app[tr.app==0]),mean(y.app[tr.app==1]))}

  Optimism.bL <- C.boot-C.orig
  Optimism.L <- apply(Optimism.bL,1,mean)
  SE.l <- sqrt(apply((Optimism.bL-Optimism.L)^2, 1, sum)/(B-1))/sqrt(B)

  C.bc <- complParam[[3]]-Optimism.L   #complParam[[3]] == c.app
  boot.info <- round(cbind(complParam[[1]], C.bc, SE.l),3)
  row.names(boot.info) <- NULL

  cat("Bias-corrected criterion values and SE are \n")
  print(boot.info)

  if(optmax){
    maxrow <- which(C.bc==max(na.omit(C.bc)))[1]
    prunebranch <- max(which(C.bc>= (C.bc[maxrow]-pp*SE.l[maxrow])))
  }else{
    maxrow <- which(C.bc==min(na.omit(C.bc)))[1]
    prunebranch <- max( which(C.bc>= (C.bc[maxrow]-pp*SE.l[maxrow])))
  }

  branches <- complParam[[1]][prunebranch]
  cat("Weakest branch is", branches, "\n")
  bestbranch <- branches

  if(bestbranch==0){prunedtree <- object}
  if(bestbranch==1){
    warning("There is no clear qualitative interaction present in the data.","\n")

    Gmat<-as.matrix(rep(1,dim(object$data)[1]))
    colnames(Gmat)<-c("1")

    prunedtree <- object

    prunedtree$si<-NULL
    prunedtree$fi<-NULL


    #create endnode information of the tree
    endinf <- matrix(0,nrow=1,ncol=10)
    endinf[,c(2:9)] <- ctmat(Gmat,y=object$data[,1],tr=object$data[,2])
    endinf <- data.frame(endinf)
    if(optmax) {endinf[,10] <- as.numeric(endinf[,3]<endinf[,6])} else {endinf[,10] <- as.numeric(endinf[,3]>=endinf[,6])}

    endinf[,1] <- 1

    rownames(endinf) <- paste("Leaf ",1,sep="")
    colnames(endinf) <- c("node","#(T=0)", "meanY|T=0", "SD|T=0","#(T=1)", "meanY|T=1","SD|T=1","d", "se", "T*")
    prunedtree$li<-endinf
    prunedtree$nind<-Gmat
    prunedtree$all<-Gmat

    warning("Best tree is the root node.")

    }
  if(bestbranch>1){
    select.prune <- which(complParam[[1]]==bestbranch)[1]
    bestbranches <- complParam[[1]][2:select.prune]
    prune.b <- vector()

    while(length(bestbranches)!=0){
      bestbranches <- sort(bestbranches)
      cut.c <- cChild(bestbranches[1],maxd)
      prune.b <- c(prune.b,bestbranches[1])
      bestbranches <- bestbranches[!bestbranches%in%cut.c]
    }

    cut.all <- unique(unlist(c(sapply(1:length(prune.b), function(v) cChild(prune.b[v],maxd)))))

    nsplit<-dim(object$si)[1]
    ytx <- object$data[,c(names(object$data)[1],names(object$data)[2],object$si[1:nsplit,3])]
    y<-as.matrix(ytx[,1])
    tr<-as.numeric(as.factor(ytx[,2]))-1

    nind <- as.numeric(colnames(object$nind)) #leaves initial tree
    nnum <- c(nind[!nind %in% cut.all], prune.b) #leaves pruned tree
    Gmat <- object$all[,paste(nnum)]

    for(i in 1:length(prune.b)){
      newleave <- object$all[,paste(nind[nind %in% cChild(prune.b[i],maxd)])]    #leaves that become the new leaf (in place of the branch)
      Gmat[,paste(prune.b[i])] <- apply(newleave,1,sum)
    }
    Lfinal <- ncol(Gmat)

    #create endnode information of the tree
    endinf <- matrix(0,nrow=length(nnum),ncol=10)
    endinf[,c(2:9)] <- ctmat(Gmat,y,tr)
    endinf <- data.frame(endinf)
    if(optmax) {endinf[,10] <- as.numeric(endinf[,3]<endinf[,6])} else {endinf[,10] <- as.numeric(endinf[,3]>=endinf[,6])}

    endinf[,1] <- nnum
    index <- leafnum(nnum)
    endinf <- endinf[index,]
    rownames(endinf) <- paste("Leaf ",1:Lfinal,sep="")
    colnames(endinf) <- c("node","#(T=0)", "meanY|T=0", "SD|T=0","#(T=1)", "meanY|T=1","SD|T=1","d", "se", "T*")

    prunedtree <- object
    prunedtree$nind <- Gmat[,order(colnames(Gmat))]
    prunedtree$fi <- object$fi
    prunedtree$li <- endinf
    prunedtree$si <- object$si[!c(object$si[,1]) %in% c(cut.all), ]
  }

  return(prunedtree)
}
