#' Control Parameters for NETT Algorithm
#'
#' Various parameters that control aspects of the \dQuote{nett} algorithm.
#'
#' @param maxd maximum depth of the tree. Default value is 4.
#' @param a1 the minimal sample size of Treatment A (\eqn{T=1}) in a leaf.
#' @param a2 the minimal sample size of Treatment B (\eqn{T=2}) in a leaf.
#' @param optmax should expected potential outcome be maximized?
#'
#' @return A list containing the options.
#'
#' @seealso \code{\link{nett}}
#' @examples
#'# Experimental condition: education vs. non-education
#'library(quint)
#'data(bcrp)
#'bcrp[,"education"] <- as.numeric(bcrp[,"cond"]==2)
#'
#'# Perform a NETT analysis
#'formula1 <- I(cesdt1 - cesdt3) ~ education | cesdt1 + negsoct1 + uncomt1 +
#'  disopt1 + comorbid + age + wcht1 + nationality + marital + trext
#'control  <- nett.control( maxd=4,a1=10, a2=10, optmax=TRUE)
#'
#' #Set the maximum depth of the tree at 6
#' control2<-nett.control(maxd=6)
#'
#' #Set minimal sample size in each treatment group at 5
#' control3<-nett.control(a1=5,a2=5)
#'
#' @import quint
#' @export
nett.control<-function(maxd = 4, a1=NULL, a2=NULL, optmax=TRUE){
  #a1=min bucket T=0, a2=minbucket T=1
  #optmax: is criterion to be maximized or minimized, default is TRUE (maximized)

  if(!is.null(a1)){
    parvec=round(c(a1, a2))
  }else{
    parvec=NULL
  }
  list(maxd = maxd, parvec = parvec, optmax=optmax)
}
