# All subfunctions relevant for displaying node information:

csigmap<-function(n1,sd1,n2,sd2, weight = TRUE){
  #computes estimate of pooled sigma , see Cohen (1988, p. 66)
  # sd1 , sd 2 and n1, n2: standard deviation and sample size for group 1 and 2 respectively
  sigmap<-sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))      #denominator: -2 !
  if (weight == FALSE) {sigmap <- sqrt((sd1^2 + sd2^2)/2)}
  return(sigmap)}

computeD<-function(n1, mean1,sd1,n2,mean2,sd2){
  #computes Cohen's D see Cohen (1988, p. 66)
  sigmap <- csigmap(n1, sd1, n2, sd2)
  dval <- (mean1 - mean2)/sigmap
  mi <- (n1 + n2 - 2)
  ni <- (n1 * n2)/(n1 + n2)
  fac <- mi/((mi - 2) * ni)
  options(warn=-1)
  cm <- (gamma(mi/2))/(sqrt(mi/2) * gamma((mi - 1)/2))
  var <- (fac * (1 + ni * dval^2)) - dval^2/cm^2
  sedval <- sqrt(var)
  obj <- list(dval = dval, se = sedval)
  return(obj)}

ctmat<-function(Gmat,y,tr){
  #creates a matrix (tmat) with final information of each terminal node of the tree (each column of Gmat)
  #tmat= I*7 matrix
  #Gmat = nodeindicator matrix
  #y = outcome variable, column vector
  #tr = treatmentvariable with two values (T=0; T=1)
  ## cardinalities t0, cardinalities t1, and mean  and var y|t=0, mean and var y|t=1
  #each row of tmat gives this information for each column of Gmat
  #thus number of rows of pmat corresponds to number of columns of Gmat
  rownum <- ncol(Gmat)
  t2mat <- matrix(0, ncol = 2, nrow = rownum)
  t1mat <- matrix(0, ncol = 2, nrow = rownum)
  tmat <- matrix(0, ncol = 2, nrow = rownum)
  dat <- data.frame(cbind(y = y, tr = tr))

  tr0 <- tr==0
  tr1 <- tr==1
  Gmat.1 <- Gmat==1

  t1mat[, 1] <- sapply(1:rownum, function(kk, Gmat, y, tr) {
    ifelse(sum(Gmat.1[, kk] & tr0) == 0, NA, mean(y[Gmat.1[,kk] & tr0]))}, Gmat = Gmat, y = y, tr = tr)
  t1mat[, 2] <- sapply(1:rownum, function(kk, Gmat, y, tr) {
    ifelse(sum(Gmat.1[, kk] & tr0) == 0, NA, sqrt(var(y[Gmat.1[, kk]& tr0])))}, Gmat = Gmat, y = y, tr = tr)
  t2mat[, 1] <- sapply(1:rownum, function(kk, Gmat, y, tr) {
    ifelse(sum(Gmat.1[, kk] & tr1) == 0, NA, mean(y[Gmat.1[,kk] & tr1]))}, Gmat = Gmat, y = y, tr = tr)
  t2mat[, 2] <- sapply(1:rownum, function(kk, Gmat, y, tr) {
    ifelse(sum(Gmat.1[, kk] & tr1) == 0, NA, sqrt(var(y[Gmat.1[,kk] & tr1])))}, Gmat = Gmat, y = y, tr = tr)
  tmat <- cbind(apply(matrix(Gmat[tr0, ],ncol = dim(Gmat)[2]), 2, sum), t1mat, apply(matrix(Gmat[tr1, ],ncol = dim(Gmat)[2]), 2, sum), t2mat)
  es <- sapply(1:rownum, function(kk, tmat) {
    ifelse(is.na(sum(tmat[kk, c(2:3, 5:6)])), NA, computeD(n1=tmat[kk,1], mean1=tmat[kk, 2], sd1=tmat[kk, 3], n2=tmat[kk, 4], mean2=tmat[kk, 5], sd2=tmat[kk, 6])$dval) }, tmat = tmat)
  se <- sapply(1:rownum, function(kk, tmat) {
    ifelse(is.na(sum(tmat[kk, c(2:3, 5:6)])), NA, computeD(tmat[kk,1], tmat[kk, 2], tmat[kk, 3], tmat[kk, 4], tmat[kk,5], tmat[kk, 6])$se)}, tmat = tmat)
  tmat <- cbind(tmat, es, se)
  colnames(tmat) <- c("nt0", "meant0", "sdt0", "nt1", "meant1","sdt1", "d", "se")
  return(as.matrix(tmat))}


cpmat<-function(Gmat,y,tr,optmax){
  #creates pmat = candidate parent nodes information matrix
  #pmat= I * 3 matrix ; I = number of candidate parent nodes
  ##what is 3?
  #columns of pmat: per node: cardinality t1, cardinality t2, cohen's d or difference in absolute means
  #Gmat = N*I matrix=indicator matrix of all candidate parent nodes
  #y = outcome variable , column vector
  #tr = treatment-variable with values 0 and 1, column vector
  rownum<-ncol(Gmat)
  pmat<-matrix(0,ncol=3,nrow=rownum)
  nt0<-apply(Gmat[tr==0,],2,sum)
  nt1<-apply(Gmat[tr==1,],2,sum)
  meant0<-numeric(rownum)
  meant1<-numeric(rownum)

  tr0 <- tr==0
  tr1 <- tr==1
  Gmat.1 <- Gmat==1

  #compute mean value of y for T=0
  meant0<-sapply(1:rownum,function(kk,Gmat,y,tr){mean(y[Gmat.1[,kk]&tr0]) },Gmat=Gmat,y=y,tr=tr)
  #compute mean value of y for T=1
  meant1<-sapply(1:rownum,function(kk,Gmat,y,tr){mean(y[Gmat.1[,kk]&tr1])},Gmat=Gmat,y=y,tr=tr)
  # col 4 van pmat geeft beste behandeling in leaf, col 3 geeft het verschil in means tussen t=1 en t=0 in een leaf

  if(optmax){
    pmat<-as.matrix(cbind(nt0, nt1, meant0-meant1, as.numeric(apply(cbind(meant0,meant1),1, which.max))-1,
                          apply(cbind(meant0,meant1) ,1, max)))
  } else {
    pmat<-as.matrix(cbind(nt0, nt1, meant0-meant1, as.numeric(apply(cbind(meant0,meant1),1, which.min))-1,
                          apply(cbind(meant0,meant1) ,1, min))) }

  colnames(pmat)[1:3] <- c("nt0","nt1","mean0-mean1")#,"best.t","Yt*")
  #pmat, for every leaf: cardinality t=0, cardinality t=1, difference in means t=0 and t=1, best treatment (0 or 1), mean outcome of best treatment in leaf

  pmat[,3][is.na(pmat[,3])]<-0
  dimnames(pmat)<-NULL
  return(pmat)}
