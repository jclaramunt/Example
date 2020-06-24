#' Summarizing NETT Tree Information
#'
#' Summary method for an object of class \code{nett}.
#'
#' @param object a \code{nett} object. This can be the output of \code{\link{nett}} or \code{\link{prune.nett}}.
#' @param digits specified number of decimal places (default is 2).
#' @param \dots optional additional arguments.
#'
#' @return prints a summarized version of the \code{nett} output.
#'
#' @details This function is a method for the generic function summary for class
#'   \code{nett}. It extracts the following essential components from a \code{nett}
#'   object: 1) Fit information;
#'   2) Split information, and 3) Leaf information.
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
#'control  <- nett.control(a1=20, a2=20, optmax=TRUE, maxd=3)
#'#fi <- nett(formula=formula1, data=subset(bcrp[1:100,],cond<3), control=control)
#'fi <- nett(formula=formula1, data=bcrp[1:145,], control=control)
#'
#'# Inspect the main results of the analysis
#'summary(fi)
#'
#' @keywords summary
#'
#' @import quint
#' @export
summary.nett<-function(object,digits=2,...){

  if(is.null(dim(object$si)[2])){
    cat("Only","root", "node.", "\n")

    cat("\n")
    cat("Leaf information:","\n")
    print(round(object$li[,c(2:10)],digits=digits))

  }else{
    #digits=number of digits at decimal points
    cat("Fit","information:", "\n")
    cat(c(rep("",15),"Criterion") ,"\n" )
    cat(c(rep("",16),paste(rep("-",4),sep="")),"\n")
    print( round(object$fi[,c(1:4)],digits=digits), row.names=FALSE)
    cat("\n")

    roundedSplitpoints <- round(suppressWarnings(as.numeric(object$si[,4])),digits)
    for(i in 1:length(roundedSplitpoints)){if(!is.na(roundedSplitpoints[i])) object$si[i,4] <- roundedSplitpoints[i]}
    object$si<-cbind(object$si[,1:3], splitpoint=object$si[,4])

    cat("Split information:","\n")
    print(object$si,row.names=TRUE)
    cat("\n")
    cat("Leaf information:","\n")
    print(round(object$li[,c(2:10)],digits=digits))
  }
}
