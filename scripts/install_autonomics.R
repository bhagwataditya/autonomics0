# Installer functions
is_package_installed <- function(x){
   x %in% rownames(installed.packages())
}

if(!is_package_installed('BiocManager')) install.packages('BiocManager')

install_if_not_available <- function(x){
   lapply(x, 
          function(x){
             if (!is_package_installed(x)) BiocManager::install(x, update = FALSE, ask = FALSE)
   })
}

check_version_compatibility <- function(){
  if(
    'ggplot2' %in% utils::installed.packages() &&
    utils::packageVersion('ggplot2') > package_version('2.2.1'))
  {
    warning("'autonomics' is currently incompatible with 'ggplot2' > v2.2.1. Downgrading.")
    remotes::install_version('ggplot2', version = '2.2.1')
  }
  
  if(
    'ggstance' %in% utils::installed.packages() &&
    utils::packageVersion('ggstance') > package_version('0.3'))
  {
    warning("'autonomics' is currently incompatible with 'ggstance' > v0.3 (which depends on 'ggplot2' v3.0). Downgrading.")
    remotes::install_version('ggstance', version = '0.3')
  } 
  
   #install.packages('assertive.reflection')
   #  
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
   
   # remotes
   install_if_not_available('remotes')

   # Check version compatibility
   check_version_compatibility()

   # autonomics.data & autonomics.support
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.data',       repos = BiocManager::repositories(), upgrade = FALSE)
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.support',    repos = BiocManager::repositories(), upgrade = FALSE)

   # autonomics.preprocess
   install_if_not_available('imputeLCMD')
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.preprocess', repos = BiocManager::repositories(), upgrade = FALSE)
   
   # autonomics.annotate & autonomics.import
   install_if_not_available(c('SummarizedExperiment', 'GenomeInfoDbData'))
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.annotate',   repos = BiocManager::repositories(), upgrade = FALSE)
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.import',     repos = BiocManager::repositories(), upgrade = FALSE)

   # autonomics.plot
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.plot',       repos = BiocManager::repositories(), upgrade = FALSE)
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.explore',    repos = BiocManager::repositories(), upgrade = FALSE)
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.find',       repos = BiocManager::repositories(), upgrade = FALSE)

   # autonomics.ora
   install_if_not_available(c('GO.db', 'PANTHER.db'))
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.ora',        repos = BiocManager::repositories(), upgrade = FALSE)
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics.integrate',  repos = BiocManager::repositories(), upgrade = FALSE)

   # autonomics
   remotes::install_github('bhagwataditya/autonomics', subdir='autonomics',            repos = BiocManager::repositories(), upgrade = FALSE)

   check_version_compatibility()
}

install_autonomics()
