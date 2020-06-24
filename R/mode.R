mode <- function(x) {
  ux <- sort(unique(x))
  if(0%in%ux) {ux <- c(0,rev(ux[-1]))}else{ux <- rev(ux)}
  ux[which.max(tabulate(match(x, ux)))]
}
