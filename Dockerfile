FROM bioconductor/bioconductor_docker
RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.data",       repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.support",    repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.annotate",   repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.import",     repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.preprocess", repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.plot",       repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.find",       repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.ora",        repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics.integrate",  repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'
RUN R -e 'remotes::install_github("bhagwataditya/autonomics/autonomics",            repos = BiocManager::repositories(), dependencies = TRUE, upgrade = FALSE)'

