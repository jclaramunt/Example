#' Predictions for new data with a NETT object
#'
#' Predicts for (new) subjects the treatment subgroups (P1, P2 or P3) based on a fitted
#'   \code{nett} object. The meaning of the subgroups are based on the two treatment categories
#'   used to fit the \code{nett} object.
#'
#' @param object an object of the class \dQuote{nett}.
#' @param newdata a data frame with data on new subjects for whom predictions should be made.
#'   The data frame should contain at least the variables used in the splits of the fitted tree.
#'   It is not necessary to include the treatment variable.
#' @param type character string denoting the type of predicted object to be returned. The default is
#'   set to \code{type="class"}: a vector with predicted treatment subgroup classes per subject
#'   is returned. If set to \code{"matrix"}, a matrix is returned with the leaf and
#'   corresponding node of the tree to which a subject is assigned.
#' @param \dots optional additional arguments.
#'
#' @seealso \code{\link{nett}}, \code{\link{prune.nett}}
#'
#' @examples
#'# Example with data from the Breast Cancer Recovery Project
#'library(quint)
#'data(bcrp)
#'set.seed(10)
#'
#'# Experimental condition: education vs. non-education
#'bcrp[,"education"] <- as.numeric(bcrp[,"cond"]==2)
#'
#'# Perform a NETT analysis
#'formula1 <- I(cesdt1 - cesdt3) ~ education | cesdt1 + negsoct1 + uncomt1 +
#'  disopt1 + comorbid + age + wcht1 + nationality + marital + trext
#'control  <- nett.control(a1=20, a2=20, optmax=TRUE, maxd=3)
#'fi <- nett(formula=formula1, data=bcrp[1:145,], control=control)
#'
#'# Prune the tree
#'fit <- prune(tree=fi, B=2,pp=1)
#'
#'#Predict the corresponding leaf of the units of a new dataset
#'prednett1<-predict(fit, newdata=subset(bcrp,cond<3))
#'prednett1
#'
#' @import quint
#' @importFrom stats as.formula model.frame
#' @importFrom utils combn
#' @useDynLib NETT
#' @exportClass nett
#' @export
predict.nett<-function(object,newdata,type='class',...){
  optmax <- object$control$optmax
  nsplit<-dim(object$si)[1]
  form<-as.formula(paste(names(object$data)[1],"~",paste(c(names(object$data)[2],object$si[1:nsplit,3]),collapse="+" )))
  ytxna <- model.frame(form, data=newdata, na.action=NULL)
  #ytx <- newdata[,c(names(object$data)[1],names(object$data)[2],object$si[1:nsplit,3])]
  ytx<-model.frame(form,data=newdata) # NA omitted data frame required for procedure
  y<-as.matrix(ytx[,1])
  tr<-as.numeric(as.factor(ytx[,2]))-1
  root<-rep(1,dim(ytx)[1]);
  if(!is.factor(ytx[,object$si[1,3]])){
    Gmat <- makeGchmat(gm=root, varx=ytx[,object$si[1,3]], splitpoint=object$si[1,4])  }
  if(is.factor(ytx[,object$si[1,3]])){
    possibleSplits <- determineSplits(x=ytx[,object$si[1,3]], gm=root)
    assigMatrix <- makeCatmat(x=ytx[,object$si[1,3]], gm=root, z=possibleSplits[[1]], splits=possibleSplits[[2]])
    Gmat <- makeGchmatcat(gm=root, splitpoint=object$si[1,4], assigMatrix=assigMatrix)  }
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
    }
  }

  #Here the differences

  End<-cpmat(Gmat,y,tr,optmax=optmax)
  obj<-cbind(nnum,End)
  obj<-as.data.frame(obj)
  names(obj)=c("Leaf","#(T=0)","#(T=1)","d","T*","Yt*")


  #check number of missings
  nmis<-sum(is.na(ytxna))
  index <- c(1:dim(ytxna)[1])
  if(nmis!= 0){
    naindex <- which(is.na(ytxna)) # which subjects have NA on required variables
    index <- c(1:dim(ytxna)[1])[-naindex] # which subjects have no NA on required variables
  }

  if(type=="matrix") {
    nodemat <- matrix(0, nrow=dim(newdata)[1], ncol=2)
    for(i in 1:length(nnum)){
      nodemat[index[which(Gmat[,i]==1)],1] <- which(object$li[,1]==nnum[i])
      nodemat[index[which(Gmat[,i]==1)],2] <- nnum[i]
    }
    if(nmis!= 0){
      nodemat[naindex,c(1:2)] <- NA
    }
    colnames(nodemat) <- c("Leaf", "Node")
    #rownames(nodemat) <- as.numeric(rownames(ytxna)) # give subjects numbers of original dataset
    return(list(PredValues=nodemat, obj=obj))
  }

  if(type=="class"){
    classmat <- numeric(dim(newdata)[1])
    for(i in 1:length(nnum)) {
      classmat[index[which(Gmat[,i]==1)]] <- object$li[which(nnum[i]==object$li[,1]),10]
    }
    if(nmis!= 0){
      classmat[naindex] <- NA
    }
    names(classmat) <- 1:dim(newdata)[1]
    return(list(PredValues=classmat, obj=obj))
  }
  #predTrt <- (Gmat %*% object$li[order(object$li$node),"T*"]) #-1 #Gmat does not give leaf or node order so first order its columns according to nodes


}
