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
         'autonomics.plot',
         'autonomics.annotate',
         'autonomics.import',
         'autonomics.support',
         'autonomics.find'),
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
update.packages(checkBuilt = TRUE, ask = FALSE)
