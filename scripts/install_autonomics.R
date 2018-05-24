# Installer functions
source("https://bioconductor.org/biocLite.R")
is_package_installed <- function(x){
   x %in% rownames(installed.packages())
}
install_if_not_available <- function(x){
   lapply(x, 
          function(x){
             if (!is_package_installed(x)) biocLite(x, suppressUpdates = TRUE)
   })
}

check_version_compatibility <- function(){
   
   install.packages('assertive.reflection')
 
   # 13 May 2018:   fixed  
   # 29 April 2018: data.table is not installing (smoothly) in R3.5.0 on windows
   #   - The binary is not yet available on CRAN
   #   - Installing from source doesn't work either, probably due to an incompatibility of Rtools 3.5 (latest) and R 3.5.0, 
   #     as can be detected with devtools::setup_rtools(). 
   #if (assertive.reflection::is_windows() & version$major=='3' & version$minor=='5.0'){
   #   stop('Roll back to R 3.4.4 before installing autonomics in windows.', 
   #        '\nThe R package data.table is not (yet) installing properly on R 3.5.0 in windows.', 
   #        '\nAnd autonomics needs data.table.')
   #}
}

install_autonomics <- function(){
   
   # Check version compatibility
   check_version_compatibility()
   
   # devtools
   install_if_not_available('devtools')

   # autonomics.data & autonomics.support
   devtools::install_github('bhagwataditya/autonomics/autonomics.data',       repos = biocinstallRepos(), upgrade_dependencies = FALSE)
   devtools::install_github('bhagwataditya/autonomics/autonomics.support',    repos = biocinstallRepos(), upgrade_dependencies = FALSE)

   # autonomics.preprocess
   install_if_not_available('imputeLCMD')
   devtools::install_github('bhagwataditya/autonomics/autonomics.preprocess', repos = biocinstallRepos(), upgrade_dependencies = FALSE)
   
   # autonomics.annotate & autonomics.import
   install_if_not_available(c('SummarizedExperiment', 'GenomeInfoDbData'))
   devtools::install_github('bhagwataditya/autonomics/autonomics.annotate',   repos = biocinstallRepos(), upgrade_dependencies = FALSE)
   devtools::install_github('bhagwataditya/autonomics/autonomics.import',     repos = biocinstallRepos(), upgrade_dependencies = FALSE)

   # autonomics.plot
   devtools::install_github('bhagwataditya/autonomics/autonomics.plot',       repos = biocinstallRepos(), upgrade_dependencies = FALSE)
   devtools::install_github('bhagwataditya/autonomics/autonomics.explore',    repos = biocinstallRepos(), upgrade_dependencies = FALSE)
   devtools::install_github('bhagwataditya/autonomics/autonomics.find',       repos = biocinstallRepos(), upgrade_dependencies = FALSE)

   # autonomics.ora
   install_if_not_available(c('GO.db', 'PANTHER.db'))
   devtools::install_github('bhagwataditya/autonomics/autonomics.ora',        repos = biocinstallRepos(), upgrade_dependencies = FALSE)
   devtools::install_github('bhagwataditya/autonomics/autonomics.integrate',  repos = biocinstallRepos(), upgrade_dependencies = FALSE)

   # autonomics
   devtools::install_github('bhagwataditya/autonomics/autonomics',            repos = biocinstallRepos(), upgrade_dependencies = FALSE)

}

install_autonomics()