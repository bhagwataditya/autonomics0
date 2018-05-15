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
         'autonomics.support',
         'autonomics.explore',
         'autonomics.plot',
         'autonomics.import',
         'autonomics.annotate',
         'autonomics.data',
         'autonomics.preprocess'),
      fsep = '/') %>%
   switchr::makeManifest(
      name = basename(.),
      url  = .,
      type = 'git') %>%
   switchr::install_packages(
      slot(., 'manifest') %>%
         magrittr::extract2('name'),
      .)

# Update everything (also already present) --------------------------------
BiocInstaller::biocLite(checkBuilt = TRUE, ask = FALSE)
