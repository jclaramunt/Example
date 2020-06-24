ifelseC = function(test, yes, no){
  requireNamespace("Rcpp")
  if(typeof(test) == "list"){
    test = unlist(test)
  }
  if(typeof(test) != "logical"){
    test = try(as.logical(test))
    stopifnot(typeof(test)== "logical")
  }

  l_test = length(test)
  l_yes = length(yes)
  l_no = length(no)
  l_max = max(l_test, l_yes, l_no)

  if(l_yes == 1){
    yes = rep(yes,l_max)
  }

  if(l_no == 1){
    no = rep(no, l_max)
  }

  type_yes = typeof(yes)
  type_no = typeof(no)
  stopifnot(type_yes == type_no)
  if(type_yes == "double"){
    out = try(ifelseCNum(test, yes, no))
  } else if(type_yes == "character"){
    out = try(ifelseCChar(test, yes, no))
  } else if(type_yes == "integer"){
    out = try(ifelseCInt(test, yes, no))
  } else if(type_yes == "logical"){
    out = try(ifelseCLogic(test, yes, no))
  } else{
    stop(paste("type", type_yes, "not supported."))
  }
  out
}
