# Intro

autonomics is an R suite for making omics data analysis easier, by **automating** the automatable, and **facilitating** the interactive.


# Installation

The **development** version is up-to-date, but not yet stable:

    # First install the R package remotes
    install.packages('remotes')
    remotes::install_github('bhagwataditya/autonomics/autonomics.data',       ref = 'dev')
    remotes::install_github('bhagwataditya/autonomics/autonomics.support',    ref = 'dev')
    remotes::install_github('bhagwataditya/autonomics/autonomics.import',     ref = 'dev')
    remotes::install_github('bhagwataditya/autonomics/autonomics.annotate',   ref = 'dev')
    remotes::install_github('bhagwataditya/autonomics/autonomics.preprocess', ref = 'dev')
    remotes::install_github('bhagwataditya/autonomics/autonomics.plot',       ref = 'dev')
    remotes::install_github('bhagwataditya/autonomics/autonomics.find',       ref = 'dev')
    remotes::install_github('bhagwataditya/autonomics/autonomics.ora',        ref = 'dev')
    remotes::install_github('bhagwataditya/autonomics/autonomics',            ref = 'dev')


The **stable** branch is error-free, but now outdated:

    # remotes::install_github('bhagwataditya/autonomics/autonomics.data'      )
    # remotes::install_github('bhagwataditya/autonomics/autonomics.support'   )
    # remotes::install_github('bhagwataditya/autonomics/autonomics.import'    )
    # remotes::install_github('bhagwataditya/autonomics/autonomics.annotate'  )
    # remotes::install_github('bhagwataditya/autonomics/autonomics.preprocess')
    # remotes::install_github('bhagwataditya/autonomics/autonomics.plot'      )
    # remotes::install_github('bhagwataditya/autonomics/autonomics.explore'   )
    # remotes::install_github('bhagwataditya/autonomics/autonomics.find'      )
    # remotes::install_github('bhagwataditya/autonomics/autonomics.ora'       )
    # remotes::install_github('bhagwataditya/autonomics/autonomics'           )

## Read omics data and prepare for analysis (dev)

    # METABOLON
          require(magrittr)
          object <- 'extdata/glutaminase/glutaminase.xlsx'    %>% 
                     system.file(package = 'autonomics.data') %>% 
                     autonomics::read_metabolon()
          object %>% autonomics::prepare_metabolon()
    
    # SOMASCAN
          object <- 'extdata/stemcomp/soma/stemcomp.adat'     %>% 
                     system.file(package = 'autonomics.data') %>% 
                     autonomics::read_somascan()
          object %>% autonomics::prepare_somascan()
    
    # RNASEQ
          object <- 'extdata/stemdiff/rnaseq/gene_counts.txt' %>% 
                     system.file(package = 'autonomics.data') %>% 
                     autonomics::read_somascan()
          object %>% autonomics::prepare_rnaseq()
    
    # PROTEINGROUPS
          object <- 'extdata/stemcomp/maxquant/proteinGroups.txt' %>% 
                     system.file(package = 'autonomics.data') %>% 
                     autonomics::read_proteingroups()
          object %>% autonomics::prepare_proteingroups()
    
    # EXIQON
          object <-  autonomics::read_exiqon(myfile)
          object %>% autonomics::prepare_exiqon()

   
    # ANY OMICS DATASET
          object <- 'extdata/glutaminase/glutaminase.xlsx'   %>% 
                     system.file(package = 'autonomics.data') %>% 
                     autonomics::read_omics(
                        sheet      = 2,
                        fid_rows   = 11:401,    fid_cols   = 5,
                        sid_rows   = 3,         sid_cols   = 15:86,
                        expr_rows  = 11:401,    expr_cols  = 15:86,
                        fvar_rows  = 10,        fvar_cols  = 1:14,
                        svar_rows  = 1:10,      svar_cols  = 14,
                        fdata_rows = 11:401,    fdata_cols = 1:14,
                        sdata_rows = 1:10,      sdata_cols = 15:86,
                        transpose  = FALSE)
                        
## Explore omics data

    # Principal Component Analysis
        object <- autonomics.data::glutaminase
        object %>% autonomics::plot_pca_samples()
        object %>% autonomics::plot_pca_features()
        object %>% autonomics::plot_pca_samples_and_features()
        
    #  Linear Discriminant Analysis
        object %>% autonomics::plot_lda_samples()
        object %>% autonomics::plot_lda_features()
        object %>% autonomics::plot_lda_samples_and_features()
        
    # Partial Least Squares Analysis
        object %>% autonomics::plot_pls_samples()
        object %>% autonomics::plot_pls_features()
        object %>% autonomics::plot_pls_samples_and_features()
   
    # Spectral Map Analysis
        object %>% autonomics::plot_sma_samples()          
        object %>% autonomics::plot_sma_features() 
        object %>% autonomics::plot_sma_samples_and_features()


## Analyze contrasts

    object <- autonomics.data::glutaminase
    table(object$subgroup)
    glutcontrasts <- c(  UT.h72   = 'UT_h72   - UT_h10',
                         uM10.h72 = 'uM10_h72 - uM10_h10')
    autonomics::contrastdefs(object) <- glutcontrasts
    object %<>% autonomics::add_limma()
    object %>% autonomics::plot_contrast_features(n=4)
    object %>% autonomics::plot_contrast_features(contrast = glutcontrasts[2], n=4)
    


