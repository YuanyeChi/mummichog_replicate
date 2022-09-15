is.blank <- function(x, false.triggers=FALSE){
  if(is.function(x)) return(FALSE) # Some of the tests below trigger
  # warnings when used on functions
  return(
    is.null(x) ||                # Actually this line is unnecessary since
      length(x) == 0 ||            # length(NULL) = 0, but I like to be clear
      all(is.na(x)) ||
      all(x=="") ||
      (false.triggers && all(!x))
  )
}

removeListElem <- function(inlist,elem_remove){
  outlist = lapply(inlist,setdiff,elem_remove)
  outlist1  = lapply(elem_remove,setdiff,inlist)
  outlist[lengths(outlist) > 0 | lengths(outlist1) > 0]
}