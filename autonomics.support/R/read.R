#' convenient fread
#' 
#' Convenience wrapper to \code{\link[data.table]{fread}} with preset arguments, making the fread operation:
#' \itemize{
#'    \item reliable:   integer64 are read as 'numeric'
#'    \item convenient: data frame returned rather than a data.table
#'    \item silent:     verbose set to FALSE and warnings are suppressed
#' }
#' 
#' The integer64 = 'numeric' specification is required to read integers as regular integers rather than
#' 64 bit integers. The problem with 64 it integers is that the numbers display incorrectly when the package 
#' bit64 is not loaded. And the package bit64 itself had some issues, at least in earlier versions.
#' 
#' @importFrom data.table fread
#' @param file Either the file name to read (containing no \\n character), a 
#' shell command that preprocesses the file (e.g. fread("grep blah filename")) 
#' or the input itself as a string (containing at least one \\n). Passed to 
#' \code{\link[data.table]{fread}}.
#' @param verbose Be chatty and report timings? Passed to 
#' \code{\link[data.table]{fread}}.
#' @param integer64 What to read columns detected as containing integers larger 
#' than 2^31 as. Passed to \code{\link[data.table]{fread}}.
#' @param data.table \code{TRUE} returns a \code{data.table}. \code{FALSE} 
#' returns a \code{data.frame}. Passed to \code{\link[data.table]{fread}}.
#' @param ... Passed to \code{\link[data.table]{fread}}.
#' @export
cfread <- function(file, verbose = FALSE, integer64  = 'numeric', data.table = FALSE, ...){
  
  dt <- suppressWarnings(
        data.table::fread(
          file, 
          data.table = data.table, 
          integer64  = integer64, 
          verbose    = verbose,
          ...
        )
      )
  
  # The integer64 = 'numeric' sometimes fails: https://github.com/Rdatatable/data.table/issues/2607
  # Enforce it.
  integer64_columns <- names(which(lapply(dt, class) == 'integer64'))
  for (col in integer64_columns) data.table::set(dt, j=col, value=as.numeric(dt[[col]]))
  
  dt
  
}
 