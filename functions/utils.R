normalise <- function(x) {
  total <- sum(x, na.rm = TRUE)
  new <- x/total
  return(new)
}