# Installer functions
install.packages(setdiff("BiocManager", installed.packages()))
is_package_installed <- function(x){
   x %in% rownames(installed.packages())
}
install_if_not_available <- function(x){
   lapply(x, 
          function(x){
             if (!is_package_installed(x)) BiocManager::install(x, update = FALSE)
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
   devtools::install_github('bhagwataditya/autonomics0/autonomics.data',       repos = BiocManager::repositories(), upgrade = FALSE)
   devtools::install_github('bhagwataditya/autonomics0/autonomics.support',    repos = BiocManager::repositories(), upgrade = FALSE)

   # autonomics.annotate & autonomics.import
   install_if_not_available(c('SummarizedExperiment', 'GenomeInfoDbData'))
   devtools::install_github('bhagwataditya/autonomics0/autonomics.annotate',   repos = BiocManager::repositories(), upgrade = FALSE)
   devtools::install_github('bhagwataditya/autonomics0/autonomics.import',     repos = BiocManager::repositories(), upgrade = FALSE)

   # autonomics.preprocess
   install_if_not_available('imputeLCMD')
   devtools::install_github('bhagwataditya/autonomics0/autonomics.preprocess', repos = BiocManager::repositories(), upgrade = FALSE)
   
   # autonomics.plot
   devtools::install_github('bhagwataditya/autonomics0/autonomics.plot',       repos = BiocManager::repositories(), upgrade = FALSE)
   devtools::install_github('bhagwataditya/autonomics0/autonomics.explore',    repos = BiocManager::repositories(), upgrade = FALSE)
   devtools::install_github('bhagwataditya/autonomics0/autonomics.find',       repos = BiocManager::repositories(), upgrade = FALSE)

   # autonomics.ora
   install_if_not_available(c('GO.db', 'PANTHER.db'))
   devtools::install_github('bhagwataditya/autonomics0/autonomics.ora',        repos = BiocManager::repositories(), upgrade = FALSE)
   devtools::install_github('bhagwataditya/autonomics0/autonomics.integrate',  repos = BiocManager::repositories(), upgrade = FALSE)

   # autonomics
   devtools::install_github('bhagwataditya/autonomics0/autonomics',            repos = BiocManager::repositories(), upgrade = FALSE)

}

install_autonomics()
