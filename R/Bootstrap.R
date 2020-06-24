Bootstrap<-function(y,B,tr){
  #stratification on y
  n<-length(y)
  vec<-1:n
  n1<-length(y[tr==0])
  n2<-length(y[tr==1])
  Bmat<-matrix(vec,nrow=n,ncol=B)
  Bmatr<-matrix(0,nrow=n,ncol=B)
  ##stratification on tr ensures that in each bootstrap sample: n1=78 and n2=70
  Bmatr[tr==0,]<-sapply(1:B,function(kk,n,Bmat,tr){sample(Bmat[tr==0,kk],n,replace=TRUE) },Bmat=Bmat,n=n1,tr=tr)
  Bmatr[tr==1,]<-sapply(1:B,function(kk,n,Bmat,tr){sample(Bmat[tr==1,kk],n,replace=TRUE) },Bmat=Bmat,n=n2,tr=tr)
  return(Bmatr)}
