#' Compute using Fortran
#'
#' Calls a Fortran subroutine to perform computation.
#'
#' @param x numeric vector input
#' @return processed numeric vector
#' @export
#' @examples
#' result <- call_fortran(c(1, 2, 3))
call_fortran <- function(x) {
  .Fortran("twcomputation", as.double(x), as.integer(length(x)))[[1]]
}   