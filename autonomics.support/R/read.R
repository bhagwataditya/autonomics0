#' Convenient fread
#' 
#' Convenient \code{\link[data.table]{fread}} wrapper which sets 
#' integer64 to 'numeric' to prevent issues with 
#' 
#' Convenience wrapper to \code{\link[data.table]{fread}} with preset arguments, 
#' making the fread operation:
#' \itemize{
#'    \item reliable:   integer64 are read as 'numeric'
#'    \item silent:     verbose set to FALSE and warnings are suppressed
#' }
#' 
#' The integer64 = 'numeric' specification is required to read integers as 
#' regular integers rather than 64 bit integers. The problem with 64 bit 
#' integers is that the numbers display incorrectly when the package 
#' bit64 is not loaded.
#' 
#' @param file        passed to \code{\link[data.table]{fread}(input = file)}
#' @param verbose     FALSE (default) or TRUE. 
#'                    Passed to \code{\link[data.table]{fread}(input = file)}
#' @param integer64  'numeric' (default) or 'integer64': passed to data.table::fread
#' @param data.table \code{TRUE} or \code{FALSE}: whether to return as
#' \code{data.table} (rather than data.frame)
#' @param ... Passed to \code{\link[data.table]{fread}}.
#' @export
cfread <- function(file, verbose = FALSE, integer64  = 'numeric', data.table = TRUE, ...){

    # Read
    dt <- suppressWarnings(
            data.table::fread(
              file, 
              data.table = data.table, 
              integer64  = integer64, 
              verbose    = verbose,
              ...))
  
    # The integer64 = 'numeric' sometimes fails: https://github.com/Rdatatable/data.table/issues/2607
    # Enforce it.
    integer64_columns <- names(which(lapply(dt, class) == 'integer64'))
    for (col in integer64_columns){
        data.table::set(dt, j=col, value=as.numeric(dt[[col]]))
    }
  
    # Return
    dt
  
}
 