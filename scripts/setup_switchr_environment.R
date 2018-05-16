# Setup 'switchr' infrastructure ------------------------------------------
switchr::switchTo("autonomics")

# Install/update gitted infrastructure ------------------------------------
# Deal with prerequisites
packages_prerequisite <- c(
   'magrittr',
   'stringi')
packages_to_install <-  setdiff(
   packages_prerequisite,
   installed.packages())
if(length(packages_to_install) > 0){
   switchr::install_packages(
      packages_to_install)
}
## Don't use magrittr functionality prior to this!
library(magrittr)

# Discover dirname
split_wd <- getwd() %>%
   stringi::stri_split_fixed(.Platform$file.sep) %>%
   unlist()

base_i <- split_wd %>%
   stringi::stri_detect_regex('^autonomics$') %>%
   which() %>%
   max(na.rm = TRUE)

autonomics_base_dir <- split_wd %>%
   magrittr::extract(seq(base_i)) %>%
   paste(collapse = .Platform$file.sep)

# Install autonomics modules
## Order matters for interdepence!
autonomics_module_packages <- c(
   "autonomics.data",
   "autonomics.ora",
   "autonomics.integrate",
   "autonomics.find",
   "autonomics.preprocess",
   "autonomics.plot",
   "autonomics.import",
   "autonomics.explore",
   "autonomics.support",
   "autonomics.annotate",
   "autonomics")

## This takes LONG when no dependencies are present yet ...
switchr::makeManifest(
   name   = autonomics_module_packages,
   url    = autonomics_base_dir,
   subdir = autonomics_module_packages,
   type   = 'local') %>% #,
   # branch = dplyr::case_when(
   #    . == 'autonomics.preprocess' ~ 'feature_sd_mean',
   #    . != 'autonomics.preprocess' ~ 'master')) %>%
   switchr::install_packages(
      slot(., 'manifest') %>%
         magrittr::extract2('name'),
      .)

# Update everything (also already present) --------------------------------
BiocInstaller::biocLite(checkBuilt = TRUE, ask = FALSE)
