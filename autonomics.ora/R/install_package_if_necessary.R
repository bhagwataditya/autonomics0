install_package_if_necessary <- function(pkg){
  if(!requireNamespace(pkg, quietly = TRUE)){
    message("The package ", pkg, " is not available; trying to install it.")
    oldOps <- options(warn = 2)
    on.exit(options(oldOps))
    BiocManager::install(pkg, update = FALSE)
  }
}


