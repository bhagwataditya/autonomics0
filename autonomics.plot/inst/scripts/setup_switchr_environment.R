library(magrittr)
# Setup 'switchr' infrastructure ------------------------------------------
getwd() %>%
   basename() %>%
   switchr::switchTo()
library(magrittr)

# Install/update gitted infrastructure ------------------------------------
'https://bitbucket.org/graumannlabtools' %>%
   file.path(
      c(
         'autonomics.preprocess',
         'autonomics.annotate',
         'autonomics.import',
         'autonomics.support',
         'autonomics.plot'),
      fsep = '/') %>%
   switchr::makeManifest(
      name = basename(.),
      url  = .,
      type = 'git') %>% #,
      # branch = dplyr::case_when(
      #    basename(.) == 'autonomics.preprocess' ~ 'feature_sd_mean',
      #    basename(.) != 'autonomics.preprocess' ~ 'master')) %>%
   switchr::install_packages(
      slot(., 'manifest') %>%
         magrittr::extract2('name'),
      .)

# Update everything (also already present) --------------------------------
BiocInstaller::biocLite(checkBuilt = TRUE, ask = FALSE)
