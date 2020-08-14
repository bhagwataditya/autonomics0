# Installer functions
is_package_installed <- function(x){
   x %in% rownames(installed.packages())
}

if (!is_package_installed('BiocManager')) install.packages('BiocManager')

install_if_not_available <- function(x){
   lapply(x, 
          function(x){
             if (!is_package_installed(x)) BiocManager::install(x,
               update = FALSE, ask = FALSE)
   })
}

check_version_compatibility <- function(){
  if (
    'ggplot2' %in% utils::installed.packages() &&
    utils::packageVersion('ggplot2') > package_version('2.2.1')) {
    warning("'autonomics' is currently incompatible with 'ggplot2' > v2.2.1. Downgrading.")
    remotes::install_version(
       'ggplot2',
       version = '2.2.1',
       repos = union(getOption("repos"), c(MRAN = "https://mran.microsoft.com")))
  }
  
  if (
    'ggstance' %in% utils::installed.packages() &&
    utils::packageVersion('ggstance') > package_version('0.3')) {
    warning("'autonomics' is currently incompatible with 'ggstance' > v0.3 (which depends on 'ggplot2' v3.0). Downgrading.")
    remotes::install_version(
       'ggstance',
       version = '0.3',
       repos = union(getOption("repos"), c(MRAN = "https://mran.microsoft.com")))
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

install_git_multiple_subdirs <- function(
  repo, subdirs, host = NULL, ref = "HEAD",
  remote_type = c("bitbucket", "git", "github", "gitlab"),
  repos = BiocManager::repositories()) {
  remote_type <- match.arg(remote_type)
  ## Download bundle (bitbucket, github, gitlab) or shallow clone (git) repo
  ## Heavily relying on (unexported) functionality from `remotes`
  if (is.null(host)) host <- switch(
    remote_type,
    bitbucket = "api.bitbucket.org/2.0",
    git = stop("`host` definition required."),
    github = "api.github.com",
    gitlab = "gitlab.com")
  remote <- switch(remote_type,
                   bitbucket = remotes:::bitbucket_remote(repo = repo, host = host, ref = ref),
                   git = remotes:::git_remote(
                     url = utils::URLencode(paste(host, repo, sep = "/")), git = "external", ref = ref),
                   github = remotes:::github_remote(repo, host = host, ref = ref),
                   gitlab = remotes:::gitlab_remote(repo, host = host, ref = ref))
  bundle <- remotes:::remote_download(remote)
  bundle_is_tarball <- grepl(pattern = "\\.tar\\.gz$", x = bundle)
  on.exit(unlink(bundle, recursive = TRUE), add = TRUE)
  ## Assure presence of requested subdirs
  if (bundle_is_tarball) {
    subdirs_present <- untar(tarfile = bundle, list = TRUE,
                             compressed = TRUE)
    subdirs_present <- unique(sapply(
      subdirs_present,
      function(x){
        unlist(strsplit(x, split = "[/\\]"))[2]
      }))
  } else {
    subdirs_present <-  list.dirs(path = bundle, recursive = FALSE)
  }
  names(subdirs_present) <- basename(subdirs_present)
  missing_subdirs <- setdiff(subdirs, names(subdirs_present))
  if (length(missing_subdirs) != 0) stop("Requested subdirs not present: ",
                                         paste(subdirs, collapse = ", "))
  ## Iteratively install the subdir-contained packages
  for (sd in subdirs) {
    remotes::install_local(path = bundle, subdir = sd, dependencies = TRUE,
                           upgrade = "never", repos = repos)
  }
}

install_autonomics <- function(ref ="HEAD"){
   
   # remotes
   install_if_not_available('remotes')

   # Check version compatibility
   check_version_compatibility()

   # Install prerequisite packages
   ## autonomics.preprocess
   install_if_not_available('imputeLCMD')
   ## autonomics.annotate & autonomics.import
   install_if_not_available(c('SummarizedExperiment', 'GenomeInfoDbData'))
   ## autonomics.ora
   install_if_not_available(c('GO.db', 'PANTHER.db'))

   # Install autonomics (THE ORDER OF SUBDIRS MATTERS)
   install_git_multiple_subdirs(repo = "bhagwataditya/autonomics",
     subdirs = c("autonomics.data", "autonomics.support",
       "autonomics.annotate", "autonomics.import", "autonomics.preprocess",
       "autonomics.plot", "autonomics.explore", "autonomics.find",
       "autonomics.ora", "autonomics.integrate", "autonomics"),
     remote_type = "github", ref = ref)
   
   check_version_compatibility()
}

install_autonomics(ref = "prod")
