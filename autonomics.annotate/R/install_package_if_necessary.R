utils::globalVariables('biocLite')

#' Install package if necessary
#' @param pkg package
#' @export
install_package_if_necessary <- function(pkg){
   if(!requireNamespace(pkg, quietly = TRUE)){
      message("The package ", pkg, " is not available; trying to install it.")
      oldOps <- options(warn = 2)
      on.exit(options(oldOps))
      source("https://bioconductor.org/biocLite.R")
      biocLite(pkg)
   }
}


