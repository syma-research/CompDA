tss_withzero <- function(x) {
  if(any(x < 0))
    stop("Data shouldn't be TSS'ed!")
  if(all(x == 0))
    return(x)
  return(x / sum(x))
}
