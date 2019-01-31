#' Add subgroup values to object
#'
#' Three methods: (1) take values from existing svar
#'                (2) guess from sampleids
#'                (3) merge designfile (with one column sample_id) into sdata
#'
#' @param object SummarizedExperiment
#' @param take_from_svar        character(1)
#' @param guess_from_sampleids  logical(1)
#' @param merge_designfile      character(1): designfile path
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- system.file('extdata/glutaminase/glutaminase.xlsx', package = 'autonomics.data')
#'    object <- file %>% autonomics.import::read_metabolon_asis()
#' }
#' @importFrom magrittr %>%
#' @export
add_subgroup <- function(
   object,
   take_from_svar        = NULL,
   guess_from_sampleids  = FALSE,
   merge_designfile      = NULL
){

   # Assert
   if (!is.null(take_from_svar)){ assertive.types::assert_is_character(take_from_svar);
      assertive.types::assert_is_subset(take_from_svar, autonomics.import::svars(object)) }
   assertive.types::assert_is_logical(guess_from_sampleids)
   if (!is.null(merge_designfile)) assertive.files::assert_all_are_existing_files(merge_design_file)

   # Either take from svar
   if (!is.null(take_from_svar)){
      autonomics.support::cmessage("Add subgroup: take from svar '%s'", take_from_svar)
      autonomics.import::sdata(object) %<>% cbind(subgroup = .[[svar]], .)

      # Or guess from sampleids
   } else if (guess_from_sampleids){
      autonomics.support::cmessage('Add subgroup: guess from sampleids')
      autonomics.import::sdata(object) %<>% cbind(subgroup = autonomics.import::guess_subgroup_values(.$sample_id), .)

      # Or merge in designfile
   } else if (merge_designfile){
      autonomics.support::cmessage("Add subgroup: merge in designfile '%s'", merge_designfile)
      design_df <- autonomics.support::cfread(merge_designfile, data.table = FALSE)
      assertive.sets::assert_is_subset('sample_id', names(design_df))
      autonomics.import::sdata(object) %<>% merge(design_df, by = 'sample_id', sort = FALSE, all.x =TRUE)
   }

   # Return
   return(object)
}

#======================
# METABOLON
#======================

#' Read metabolon
#' @param file                          metabolon xlsx file
#' @param sheet                         xls sheet name  or number
#' @param log2transform                 logical: whether to log2 transform
#' @param qirlc_consistent_nondetects   logical(1)
#' @param add_kegg_pathways             logical(1): whether to add KEGG pathways to fdata
#' @param add_smiles                    logical(1): whether to add SMILES to fdata
#' @return SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/glutaminase/glutaminase.xlsx' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::read_metabolon()
#' }
#' if (require(subramanian.2016)){
#'    file <- 'extdata/metabolon/subramanian.2016.metabolon.xlsx' %>%
#'            system.file(package='subramanian.2016')
#'    object <- file %>% autonomics.import::load_metabolon(sheet = 5)
#' }
#' @importFrom magrittr %>%
#' @export
read_metabolon <- function(
   file,
   sheet                       = 2,
   fid_var                     = 'COMP_ID',
   sid_var                     = 'CLIENT_IDENTIFIER',
   log2_transform              = TRUE,
   qirlc_consistent_nondetects = FALSE,
   add_kegg_pathways           = FALSE,
   add_smiles                  = FALSE
){
   # Satisfy CHECK
   . <- NULL
   
   # Read
   all_sheets <- readxl::excel_sheets(file)
   cur_sheet <- all_sheets %>% (function(x){ names(x) <- x; x}) %>% magrittr::extract2(sheet)
   autonomics.support::cmessage('\t\tRead  %s  %s', basename(file), cur_sheet)
   object <- file %>% autonomics.import::read_metabolon_asis(sheet = sheet, fid_var = fid_var, sid_var = sid_var)
   
   # exprs
   if (log2transform) object %<>% autonomics.import::log2transform(verbose = TRUE)
   
   # Impute consistent nondetects (in intensity data)
   if (qrilc_consistent_nondetects) object %<>% autonomics.preprocess::qrilc_consistent_nondetects()
   
   # fdata
   if (add_kegg_pathways){
      autonomics.support::cmessage('\tAdd KEGG pathways to fdata')
      object %<>% autonomics.import::add_kegg_pathways_to_fdata()
   }
   if (add_smiles){
      autonomics.support::cmessage('\tAdd SMILES to fdata')
      object %<>% autonomics.import::add_smiles_to_fdata()
   }
   
   # Return
   object
}




#==========================================================
# RNASEQ
#==========================================================

#' Load RNAseq cpm
#' @param file    rnaseq counts file
#' @param fid_var feature id variable
#' @param filter_exprs_replicated_in_some_subgroup logical(1)
#' @param plot logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemdiff/rnaseq/gene_counts.txt' %>% 
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics::read_rnaseq()
#' }
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/rnaseq/gene_counts.txt', package = 'subramanian.2016')
#'    file %>% autonomics::read_rnaseq()
#' }
#' @export
read_rnaseq <- function(
   file,
   fid_var = 'gene_id', 
   filter_exprs_replicated_in_some_subgroup = FALSE, 
   plot = TRUE
){
   
   # Read
   autonomics.support::cmessage('\t\tRead counts')
   object <- autonomics.import::read_rnaseq_asis(file, fid_var = fid_var)
   
   # Sdata
   autonomics.import::sdata(object)$subgroup <- object %>% autonomics.import::guess_subgroup_values(verbose = TRUE)

   # Exprs
     # Filter nonzero in some sample
       object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample(verbose = TRUE)
     # Store counts
       autonomics.support::cmessage('\t\tStore counts in autonomics.import::counts(object)')
       autonomics.import::counts(object) <- autonomics.import::exprs(object)
     # TMM-normalize
       autonomics.import::exprs(object) %<>% autonomics.import::counts_to_cpm()
     # Filter for replications
       if (filter_exprs_replicated_in_some_subgroup) object %<>% autonomics.import::filter_exprs_replicated_in_some_subgroup('>', 1)
     # Log2 transform
       message('\t\tLog2 transform exprs: cpm -> log2cpm')
       object %<>% autonomics.import::log2transform(verbose = FALSE)
     # Add precision weights
       SummarizedExperiment::assays(object)$weights <- object %>% autonomics.import::compute_precision_weights(plot = plot, verbose = TRUE)
   
   # Return
   object
   
   # Note: Quantile normalize?
   # Gordon:  prefer TMM over quantile normalization for most cases
   #          Use quantile normalization for extreme cases
   #          Don't use both
   # object %<>% limma::normalizeBetweenArrays(method = normalize.method)
   # Don't quantile normalize when using TMM
   # https://support.bioconductor.org/p/77664/
}



#===========================================================
# EXIQON
#===========================================================

#' Read exiqon and prepare for analysis
#' @param file                         character(1)
#' @param filter_features              character(1): fvar expression on which to filter features
#' @param align_sample_means           logical(1)
#' @param lod                          numeric(1): Ct value beyond which to consider signal as NA
#' @param qrilc_consistent_nondetects  logical(1)
#' @param filter_conserved_in_human    logical(1)
#' @param plot                         logical(1)
#' @param verbose                      logical(1)
#' @examples
#' require(magrittr)
#' if (require(subramanian.2016)){
#'    file <- system.file('extdata/exiqon/subramanian.2016.exiqon.xlsx', package = 'subramanian.2016')
#'    file %>% autonomics::read_exiqon(lod = 36)
#' }
#' @importFrom magrittr %>%
#' @export
read_exiqon <- function(
   file,
   filter_features             = '`#RefGenes`==0 & `#Spike`   ==0',
   align_sample_means          = TRUE,
   lod                         = 37,
   qrilc_consistent_nondetects = TRUE,
   filter_conserved_in_human   = FALSE,
   plot                        = TRUE, 
   verbose                     = TRUE
){
   
   # Read 
   #-----
   object <- autonomics.import::read_exiqon_asis(file)
   
   # Plot function
   #--------------
   plotfun <- function(object, xlab, title){
               object %>% 
               autonomics.plot::plot_overlayed_sample_distributions() +
               ggplot2::xlab(xlab) +
               ggplot2::ggtitle(title)}
   
   # sdata
   #------
   object$subgroup <- object$sample_id %>% autonomics.import::guess_subgroup_values(verbose = verbose)
   
   # exprs
   #------
     # Original
       if (plot) object %>% plotfun('Ct', 'Load Ct values') %>% print()
     # Filter features
       object %<>% autonomics.import::filter_features_(filter_features, verbose = TRUE)
     # Align sample means
       if (align_sample_means){
         autonomics.support::cmessage('\t\tAlign sample means')
         sample_means <- object %>% autonomics.import::exprs() %>% (function(y){y[y>32] <- NA; y}) %>% colMeans(na.rm = TRUE)
         sample_diffs <- sample_means - median(sample_means)
         autonomics.import::exprs(object) %<>% sweep(2, sample_diffs)
         if (plot) object %>% plotfun('Ct', 'Align sample means') %>% print()
       }
     # NA beyond lod
       object %<>% autonomics.preprocess::na_exprs_weakly_gt(lod)
       if (plot) object %>% plotfun('Ct', 'NA Ct >= lod') %>% print()
     # Invert scale
       autonomics.support::cmessage('\t\tInvert exprs = %d - Ct', lod)
       autonomics.import::exprs(object) %<>% magrittr::subtract(lod, .)
       newmetric <- sprintf('%d - Ct', lod)
       if (plot) object %>% plotfun(newmetric, newmetric) %>% print()
     # Zero consistent NAs
       if (qrilc_consistent_nondetects){
         object %<>% autonomics.preprocess::qrilc_consistent_nondetects()
         if (plot) object %>% plotfun(newmetric, 'Zero consistent NA values') %>% print()
       }
     # Filter for mirs conserved in human
       if (filter_conserved_in_human){
         autonomics.import::fdata(object) %<>% (function(x){x$number <- x$feature_id %>% substr(9, nchar(.)); x})
         conserved_human_mirs <- autonomics.annotate::load_targetscan('H.sapiens') %>% magrittr::extract2('number')
         object %<>% autonomics.import::filter_features(number %in% conserved_human_mirs, verbose = TRUE)
       }
   
   # Return
   #-------
   object
   
}
   


#============================================================
# SOMASCAN
#============================================================

#' Read somascan and prepare for analysis
#' @param file character(1)
#' @param fid_var  character(1): feature id variable
#' @param sid_var  character(1): sample id variable
#' @param filter_sample_type         character(.): sample  types to be filtered for.     Probably a subset of c('Sample', 'QC', 'Buffer', 'Calibrator').
#' @param filter_feature_type        character(.): feature types to be filtered for.     Probably a subset of c('Protein', 'Hybridization Control Elution', 'Rat Protein').
#' @param filter_sample_quality      character(.): sample  qualities to be filtered for. Probably a subset of c('PASS', 'FLAG', 'FAIL')
#' @param filter_feature_quality     character(.): feature qualities to be filtered for. Probably a subset of c('PASS', 'FLAG', 'FAIL')
#' @param filter_na_svars            logical(1)  : whether to rm NA svars
#' @param filter_single_value_svars  logical(1)  : whether to rm single value svars
#' @param infer_design               logical(1)  : whether to infer design from sampleids
#' @param log2_transform             logical(1)  : whether to log2 transform
#' @return Summarizedexperiment
#' @examples
#' require(magrittr)
#' if (require(atkin.2014)){
#'    file <- system.file('extdata/soma/WCQ-18-007.hybNorm.plateScale.medNorm.calibrate.20181004.adat',
#'                         package = 'atkin.2014')
#'    file %>% autonomics.import::read_somascan_asis()
#'    file %>% autonomics.import::read_somascan()
#' }
#' @importFrom magrittr %>%
#' @export
read_somascan <- function(
   file,
   fid_var                 = 'SeqId',
   sid_var                 = 'SampleId',
   filter_sample_type      = 'Sample',
   filter_feature_type     = 'Protein',
   filter_sample_quality   = c('FLAG', 'PASS'),
   filter_feature_quality  = c('FLAG', 'PASS'),
   rm_na_svars             = TRUE,
   rm_single_value_svars   = TRUE,
   infer_design            = FALSE,
   log2_transform          = TRUE
){
   # Assert
   #-------
   assertive.types::assert_is_character(c(filter_sample_type, filter_feature_type, filter_sample_quality, filter_feature_quality))
   assertive.types::assert_is_logical(c(rm_na_svars, rm_single_value_svars, infer_design, log2_transform))

   # Read
   #-----
   autonomics.support::cmessage('\t\tRead %s', file)
   object <- file %>% autonomics.import::read_somascan_asis(fid_var = fid_var, sid_var = sid_var)

   # Filter
   #-------
   if ('SampleType' %in% autonomics.import::svars(object)){ # sample type - older versions don't have it
      message('\t\t=========================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Sample types` = autonomics.import::sdata(object)$SampleType))
      object %<>% autonomics.import::filter_samples(SampleType %in% filter_sample_type, verbose = TRUE)
   }
   if ('RowCheck'   %in% autonomics.import::svars(object)){ # sample quality
      message('\t\t=========================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Sample qualities ("RowCheck")` = autonomics.import::sdata(object)$RowCheck))
      object %<>% autonomics.import::filter_samples(RowCheck %in% filter_sample_quality, verbose = TRUE)
   }
   if ('Type'       %in% autonomics.import::fvars(object)){ # feature type
      message('\t\t=========================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Type` = autonomics.import::fdata(object)$Type))
      object %<>% autonomics.import::filter_features(Type %in% filter_feature_type, verbose = TRUE)
   }
   if ('ColCheck'   %in% autonomics.import::fvars(object)){ # feature quality
      message('\t\t=========================================================================')
      autonomics.support::cmessage_df('\t\t%s', table(`Feature qualities ("ColCheck")` = autonomics.import::fdata(object)$ColCheck))
      object %<>% autonomics.import::filter_features(ColCheck %in% filter_feature_quality, verbose = TRUE)
      message('\t\t=========================================================================')
   }

   # Design
   #-------
   if (autonomics.support::all_svalues_available(object, 'SampleGroup')){
      autonomics.import::sdata(object) %<>% cbind(subgroup = .$SampleGroup, .)
   }

   # Select
   #------
   if (rm_na_svars)            autonomics.import::sdata(object) %<>% autonomics.support::rm_na_columns()
   if (rm_single_value_svars)  autonomics.import::sdata(object) %<>% autonomics.support::rm_single_value_columns()

   # Log2 transform
   #---------------
   if (log2_transform) object %<>% autonomics.import::log2transform(verbose = TRUE)


   # Return
   #-------
   object
}


#===============================================
# LCMS PROTEINGROUPS and PHOSPHOSITES
#===============================================

#' Clean maxquant snames
#'
#' For charactervector or SummarizedExperiment
#'
#' Drop "Ratio normalized", "LFQ intensity" etc from maxquant snames & sample_id values
#'
#' @param x        character(.) or SummarizedExperiment
#' @param verbose  logical(1)
#' @examples
#' require(magrittr)
#'
#' # character vector
#'     autonomics.import::clean_maxquant_snames("Ratio M/L normalized STD(L)_EM00(M)_EM01(H)_R1")
#'     autonomics.import::clean_maxquant_snames("Ratio M/L STD(L)_EM00(M)_EM01(H)_R1")
#'     autonomics.import::clean_maxquant_snames('LFQ intensity STD_R1')
#'     autonomics.import::clean_maxquant_snames('LFQ intensity L STD(L)_EM00(M)_EM01(H)_R1')
#'     autonomics.import::clean_maxquant_snames('Reporter intensity 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1')
#'     autonomics.import::clean_maxquant_snames('Reporter intensity corrected 0 A(0)_B(1)_C(2)_D(3)_E(4)_F(5)_R1')
#'
#' # SummarizedExperiment
#' if (require(autonomics.data)){
#'      x <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'            system.file(package = 'autonomics.data')     %>%
#'            autonomics.import::read_proteingroups_asis()
#'      x %>% autonomics.import::clean_maxquant_snames(verbose = TRUE)
#' }
#' @importFrom magrittr %>%
#' @export
clean_maxquant_snames <- function (x, ...) {
   UseMethod("clean_maxquant_snames", x)
}


#' @importFrom magrittr %>%
#' @export
#' @rdname clean_maxquant_snames
clean_maxquant_snames.character <- function(
   x,
   quantity = autonomics.import::guess_maxquant_quantity(x),
   verbose  = FALSE
){
   # x = mix + channel. Return mix if single channel.
   pattern <- autonomics.import::maxquant_patterns %>% magrittr::extract2(quantity)
   channel <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
   mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$2')
   if (all(channel=='')){ cleanx <- mix
   } else               { cleanx <- sprintf('%s{%s}', mix, channel)
   }
   autonomics.support::cmessage('\t\tClean snames: %s  ->  %s', x[1], cleanx[1])
   return(cleanx)
}

#' @importFrom magrittr %>%
#' @export
#' @rdname clean_maxquant_snames
clean_maxquant_snames.SummarizedExperiment <- function(
   x,
   quantity = autonomics.import::guess_maxquant_quantity(x),
   verbose  = FALSE
){
   newsnames <- autonomics.import::snames(x) %>% clean_maxquant_snames(quantity = quantity, verbose=verbose)
   autonomics.import::snames(x) <- autonomics.import::sdata(x)$sample_id <- newsnames
   x
}


#' Uniplex snames
#'
#' For charactervector or SummarizedExperiment
#'
#' @param x        character vector or SummarizedExperiment
#' @param verbose  logical(1)
#' @examples
#' require(magrittr)
#'
#' # character vector
#'     # Alternate multiplexing forms supported
#'      autonomics.import::uniplex_snames("STD(L)_EM00(M)_EM01(H)_R1{M/L}")     # Label Ratio
#'      autonomics.import::uniplex_snames('A(0)_B(1)_C(2)_D(3)_R1{0}'     )     # Reporter intensity
#'      autonomics.import::uniplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')       # Label Intensity
#'
#'    # Alternate separators supported
#'      autonomics.import::uniplex_snames('STD(L)_EM00(M)_EM01(H)_R1{L}')       # underscore
#'      autonomics.import::uniplex_snames('STD(L).EM00(M).EM01(H).R1{L}')       # dot
#'      autonomics.import::uniplex_snames('STD(L)EM00(M)EM01(H).R1{L}')         # no separator
#'
#'    # Composite snames supported
#'      autonomics.import::uniplex_snames("WT.t0(L)_WT.t1(M)_WT.t2(H)_R1{H/L}") # composite snames
#'
#'    # Uniqueness ensured by appending labels when necessary
#'      autonomics.import::uniplex_snames(c("STD(L).BM00(M).BM00(H).R10{M/L}",  # implicit uniquification
#'                                          "STD(L).BM00(M).BM00(H).R10{H/L}"))
#'    # Uniplexed snames are returned unchanged
#'      autonomics.import::uniplex_snames(c('STD_R1', 'EM0_R1'))
#'
#' # SummarizedExperiment
#'   if (require(autonomics.data)){
#'      x <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'            system.file(package = 'autonomics.data')     %>%
#'            autonomics.import::read_proteingroups_asis() %>%
#'            autonomics.import::clean_maxquant_snames()
#'      x %>% autonomics.import::uniplex_snames(verbose = TRUE)
#'   }
#'
#' @export
uniplex_snames <- function (x, ...) {
   UseMethod("uniplex_snames", x)
}

#' @rdname uniplex_snames
#' @export
uniplex_snames.character <- function(x, verbose = FALSE){
   
   # Return unchanged if not multiplexed
   # KD(H)WT(L){H/L}
   pattern <- '(.+)\\{(.+)\\}'
   n_open   <- x %>% stringi::stri_count_fixed('(')
   n_closed <- x %>% stringi::stri_count_fixed(')')
   is_multiplexed <- all(stringi::stri_detect_regex(x, pattern) & (n_open==n_closed) & (n_open>0))
   if (!is_multiplexed) return(x)
   
   # Separate mix and channel
   mix     <- x %>% stringi::stri_replace_first_regex(pattern, '$1')
   channel <- x %>% stringi::stri_replace_first_regex(pattern, '$2')
   
   # Separate labels and samples
   pattern <- '\\(.+?\\)'
   labels  <- mix %>% stringi::stri_extract_all_regex(pattern) %>% lapply(stringi::stri_replace_first_fixed, '(', '') %>%
      lapply(stringi::stri_replace_first_fixed, ')', '')
   samples <- mix %>% stringi::stri_split_regex(pattern) %>%
      # rm sep from samples (but not from replicate - needed to glue back later!)
      lapply(function(y){y[1:length(labels[[1]])] %<>% stringi::stri_replace_first_regex('^[_. ]', ''); y})
   
   # Return unchanged if mixes differ in no of labels or samples
   are_all_identical <- function(y) if (length(y)==1) TRUE else all(y[-1] == y[1])
   n_samples <- vapply(samples,  length, integer(1))
   n_labels  <- vapply(labels, length, integer(1))
   if (!are_all_identical(n_samples) | !are_all_identical(n_labels)){
      autonomics.support::cmessage('\t\tCannot demultiplexing snames: mixes differ in number of samples or labels')
      return(x)
   }
   
   # Extract replicate
   n_samples %<>% unique()
   n_labels  %<>% unique()
   if (n_samples > n_labels){ replicate <- mix %>% stringi::stri_split_regex(pattern) %>% vapply((function(y) y %>% magrittr::extract(length(y))), character(1))
   samples %<>% lapply(extract, 1:(n_samples-1))
   } else {                   replicate <- rep('', length(samples))
   }
   
   # Extract channel samples from mix
   is_ratio <- channel %>% stringi::stri_detect_fixed('/') %>% all()
   samples %<>% mapply(set_names, ., labels, SIMPLIFY = FALSE)
   if (is_ratio){
      num_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(extract, character(1), 1)
      den_label <- channel %>% stringi::stri_split_fixed('/') %>% vapply(extract, character(1), 2)
      den_samples <- mapply(extract, samples, den_label)
      num_samples <- mapply(extract, samples, num_label)
      xdemultiplex <- sprintf('%s_%s%s', num_samples, den_samples, replicate)
   } else {
      samples %<>% mapply(extract, ., channel)
      xdemultiplex <- sprintf('%s%s', samples, replicate)
   }
   if (verbose) autonomics.support::cmessage('\t\tDemultiplex snames: %s  ->  %s', x[1], xdemultiplex[1])
   
   # Ensure uniqueness. Add labels if required.
   idx <- autonomics.support::cduplicated(xdemultiplex) %>% which()
   if (length(idx)>0){
      label_tags <- channel[idx] %>% stringi::stri_replace_first_fixed('/', '')
      if (verbose)   autonomics.support::cmessage('\t\tUniquify snames: %s -> %s%s (for %d/%d snames)',
                                                  xdemultiplex[idx][1], xdemultiplex[idx][1], label_tags[1],
                                                  length(idx), length(xdemultiplex))
      xdemultiplex[idx] %<>% paste0(label_tags)
   }
   
   # Return
   return(xdemultiplex)
}

#' @rdname uniplex_snames
#' @importFrom magrittr %>%
#' @export
uniplex_snames.SummarizedExperiment <- function(
   x,
   verbose  = FALSE
){
   newsnames <- autonomics.import::snames(x) %>% autonomics.import::uniplex_snames(verbose = verbose)
   autonomics.import::snames(x) <- autonomics.import::sdata(x)$sample_id <- newsnames
   x
}


#' Deconvolute proteingroups
#' @param object             SummerizedExperiment with proteinGroups data
#' @param fastafile          path to fastafile
#' @param fastafields        character vector: fields to load from fastafile
#' @param drop_isoform_info  logical: whether to drop isoform info
#' @return deconvoluted and annotated SummarizedExperiment
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    fastafile <- '../data/uniprot_hsa_20140515.fasta'
#'    if (file.exists(fastafile)){
#'       object <- 'extdata/stemcomp/maxquant/proteinGroups.txt'            %>%
#'                  system.file(package='autonomics.data')                  %>%
#'                  load_proteingroups(infer_design_from_sampleids = TRUE)
#'       autonomics.import::fdata(object) %>% head()
#'
#'       object %<>% magrittr::extract(1:100, )                             %>%
#'                   deconvolute_proteingroups(fastafile = fastafile)
#'       autonomics.import::fdata(object) %>% head()
#'    }
#' }
#' @importFrom magrittr %>% %<>%
#' @export
deconvolute_proteingroups <- function(
   object,
   fastafile,
   fastafields = c('GENES', 'PROTEIN-NAMES', 'EXISTENCE', 'REVIEWED'),
   drop_isoform_info = FALSE
){
   # Satisfy CHECK
   EXISTENCE <- GENES <- IS.FRAGMENT <- ISOFORM <- N <- NGENE <- NISOFORMS <- NPERACCESSION <- NULL
   `PROTEIN-NAMES` <- REVIEWED <- `Uniprot accessions` <- ngene <- nprotein <- nseq <- .SD <- NULL

   # Load fasta annotations
   autonomics.support::cmessage('\t\tLoad fasta file')
   fasta_annotations <- fastafile %>% autonomics.annotate::load_uniprot_fasta_annotations(fastafields)

   # Uncollapse
   fdata1 <- object %>%
      autonomics.import::fdata() %>%
      magrittr::extract(, c('feature_id', 'Uniprot accessions'), drop = FALSE) %>%
      tidyr::separate_rows(`Uniprot accessions`, sep = ';') %>%
      data.table::data.table()

   # Split into CANONICAL and isoform
   fdata1 %<>% magrittr::extract(, ISOFORM     := `Uniprot accessions`) %>%
      magrittr::extract(, `Uniprot accessions` := `Uniprot accessions` %>% stringi::stri_replace_first_regex('[-][0-9]+',     '')) # %>%
   #magrittr::extract(, ISOFORM := ISOFORM %>% sort() %>% unique() %>% paste0(collapse=';'), by = c('feature_id', 'Uniprot accessions')) %>%
   #unique()

   # Merge in uniprot fasta annotations
   nunmapped <- fdata1 %>%
      magrittr::extract(ISOFORM %in% setdiff(ISOFORM, fasta_annotations$UNIPROTKB)) %>%
      magrittr::extract(, .SD[1], by = 'feature_id') %>%
      nrow()
   autonomics.support::cmessage('\t\tDeconvolute %d/%d proteingroups with sequences from fastafile',
                                nrow(object) - nunmapped, nrow(object))
   fdata1 %<>% merge(fasta_annotations, by.x = 'Uniprot accessions', by.y = 'UNIPROTKB', sort = FALSE)
   report_n <- function(dt, prefix='', suffix=''){
      n <- dt %>% magrittr::extract(, .SD[, list(nseq = .N,
                                                 ngene = length(unique(GENES)),
                                                 nprotein = length(unique(`PROTEIN-NAMES`)))],
                                    by = 'feature_id') %>%
         magrittr::extract(, list(ngroups = .N,
                                  nsinglegene    = sum(ngene==1),
                                  nsingleprotein = sum(nprotein==1),
                                  nsingleseq     = sum(nseq==1)))
      autonomics.support::cmessage('\t\t\t%s%d proteingroups -> %d singlegene -> %d singleprotein -> %d singleseq%s',
                                   stringi::stri_pad_right(prefix, width = 60), n$ngroups,          n$nsinglegene,   n$nsingleprotein,   n$nsingleseq, suffix)
   }
   message('')
   fdata1 %>% report_n(prefix = 'All proteingroups')

   # Prefer best existence
   fdata1 %<>% magrittr::extract(, .SD[EXISTENCE == min(EXISTENCE)], by = 'feature_id')
   fdata1 %>% report_n(prefix = 'Per group: drop inferior existences')

   # Drop trembl entries from swissprot groups
   fdata1 %<>% magrittr::extract(, .SD[REVIEWED == max(REVIEWED)], by = 'feature_id')
   swissprot <- fdata1[REVIEWED==1]
   trembl    <- fdata1[REVIEWED==0]
   fdata1 %>% report_n(prefix = 'Per group: drop trembl when swissprot available')

   # trembl groups
   #--------------
   if (nrow(trembl)>0){
      message('')
      trembl  %>% report_n(prefix = 'Trembl groups')

      # Drop fragments when full sequences available
      trembl  %>% magrittr::extract(, IS.FRAGMENT := 0)
      trembl  %>% magrittr::extract(, IS.FRAGMENT:= `PROTEIN-NAMES` %>% stringi::stri_detect_fixed('(Fragment)') %>% as.numeric())
      trembl %<>% magrittr::extract(, .SD[IS.FRAGMENT == min(IS.FRAGMENT)], by = 'feature_id')
      trembl[, IS.FRAGMENT:=NULL]
      trembl  %>% report_n(prefix = 'Per group: drop fragments when full available')

      # Use first sequence per gene
      trembl  %>% magrittr::extract(, N     := .N,                    by = 'feature_id')
      trembl  %>% magrittr::extract(, NGENE := length(unique(GENES)), by = 'feature_id')
      trembl %<>% magrittr::extract(, .SD[1],                         by = c('feature_id', 'GENES'))
      trembl  %>% report_n(prefix = 'Per group/gene: use first accession')
   }

   # swissprot groups
   #-----------------
   if (nrow(swissprot)>0){
      message('')
      swissprot  %>% report_n(prefix = 'Swissprot groups')

      # swissprot groups: collapse similar isoforms (shared accession)
      swissprot  %>% magrittr::extract(, ISOFORM := ISOFORM %>% paste0(collapse = ';'), by = c('feature_id', 'Uniprot accessions'))
      swissprot %<>% unique()
      swissprot  %>% report_n(prefix = 'Per group/accession: collapse spliceforms')

      # swissprot groups: collapse dissimilar isoforms: retain accession with maximum isoforms
      swissprot  %>% magrittr::extract(, NPERACCESSION := ISOFORM %>% stringi::stri_count_fixed(';'), by = c('feature_id', 'GENES'))
      swissprot  %>% magrittr::extract(, ISOFORM       := ISOFORM %>% paste0(collapse = ';'),         by = c('feature_id', 'GENES'))
      swissprot  %>% magrittr::extract(, `PROTEIN-NAMES` %>% autonomics.support::commonify_strings(), by = c('feature_id', 'GENES'))
      swissprot %<>% magrittr::extract(, .SD[NPERACCESSION==max(NPERACCESSION)],                      by = c('feature_id', 'GENES'))
      swissprot  %>% report_n(prefix = 'Per group/gene: use spliceform with most accessions')
      swissprot %<>% magrittr::extract(, .SD[1],                                                      by = c('feature_id', 'GENES'))
      swissprot  %>% report_n(prefix = 'Per group/gene: use first spliceform')
   }

   # paralogs
   #---------
   # Collapse paralogs: choose gene with most isoforms
   message('')
   swissprot[, NPERACCESSION:=NULL]
   trembl[, N:=NULL]
   trembl[, NGENE:=NULL]
   fdata1 <- rbind(swissprot, trembl)
   monologs <- fdata1[, .SD[length(unique(GENES))==1], by = c('feature_id')]
   paralogs <- fdata1[, .SD[length(unique(GENES))>1], by = c('feature_id')]
   if (nrow(paralogs)>0){
      paralogs  %>% report_n(prefix = 'Paralog groups')
      paralogs  %>% magrittr::extract(,  NISOFORMS := 1+ISOFORM %>% stringi::stri_count_fixed(';'))
      paralogs  %>% magrittr::extract(,  GENES          := GENES           %>% paste0(collapse = ';'),                  by = 'feature_id')
      paralogs  %>% magrittr::extract(,  ISOFORM        := ISOFORM         %>% paste0(collapse = ';'),                  by = 'feature_id')
      paralogs  %>% magrittr::extract(, `PROTEIN-NAMES` := `PROTEIN-NAMES` %>% autonomics.support::commonify_strings(), by = 'feature_id')
      paralogs %<>% magrittr::extract(, .SD[NISOFORMS == max(NISOFORMS)], by = 'feature_id')
      paralogs  %>% report_n(prefix = 'Per group: use paralog with most spliceforms')
      paralogs %<>% magrittr::extract(, .SD[1], by = 'feature_id')
      paralogs  %>% report_n(prefix = 'Per group: use first paralog')
      paralogs  %>% magrittr::extract(, NISOFORMS := NULL)
   }
   fdata1 <- rbind(monologs, paralogs)

   message('')
   fdata1 %>% report_n(prefix = 'All groups (deconvoluted)')

   # Merge back
   nullify_fvars <- function(object, fvars){
      for (curfvar in fvars)   autonomics.import::fdata(object)[[curfvar]] <- NULL
      return(object)
   }
   object %<>% nullify_fvars(fvars = c('Uniprot accessions', 'Protein names', 'Gene names'))
   autonomics.import::fdata(object) %<>% merge(fdata1, by = 'feature_id', sort = FALSE, all.x = TRUE)

   # Rename (MaxQuant style)
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('GENES',         'Gene names')
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('PROTEIN-NAMES', 'Protein names')
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('ISOFORM',       'Isoforms')

   # Remove unimportant fvars
   autonomics.import::fdata(object)$REVIEWED  <- NULL
   autonomics.import::fdata(object)$EXISTENCE <- NULL
   if (drop_isoform_info) autonomics.import::fdata(object)$Isoforms <- NULL

   # Return
   object

}


#' Read proteingroups and prepare for analysis
#' @param file                    character(1)
#' @param quantity                character(1): any value in names(maxquant_patterns)
#' @param fvars                   character(n): names of fvar columns
#' @param rm_reverse              logical(1)
#' @param rm_contaminants         logical(1)
#' @param rm_complete_nondetects  logical(1)
#' @param clean_snames            logical(1)
#' @param uniplex_snames          logical(1)
#' @param log2transform           logical(1)
#' @param fastafile               NULL or character(1): if provided, heuristic proteingroups deconvolution is attempted
#' @param verbose                 logical(1)
#' @examples
#' require(magrittr)
#' if (require(autonomics.data)){
#'    file <- 'extdata/stemdiff/maxquant/proteinGroups.txt' %>%
#'             system.file(package = 'autonomics.data')
#'    file %>% autonomics.import::read_proteingroups_asis()
#' }
#' if (require(graumann.lfq)){
#'    file <- system.file('extdata/proteinGroups.txt', package = 'graumann.lfq')
#'    file %>% read_proteingroups()
#' }
#' @importFrom magrittr %>% %<>%
#' @export
read_proteingroups <- function(
   file,
   quantity               = NULL,
   fvars                  = autonomics.import::proteingroups_fvars,
   rm_reverse             = TRUE,
   rm_contaminants        = TRUE,
   rm_complete_nondetects = TRUE,
   clean_snames           = TRUE,
   uniplex_snames         = TRUE,
   log2transform          = TRUE,
   fastafile              = NULL,
   verbose                = TRUE
){

   # Assert
   assertive.types::assert_is_logical(c(rm_reverse, rm_contaminants, rm_complete_nondetects, 
                                        clean_snames, uniplex_snames, log2transform, verbose))
   if (!is.null(fastafile)) assertive.files::assert_all_are_existing_files(fastafile)
   
   # Read
   if (verbose) autonomics.support::cmessage('\tRead proteinGroups.txt')
   object <- file %>% autonomics.import::read_proteingroups_asis(quantity       = quantity,
                                                                 fvars          = fvars,
                                                                 verbose        = verbose)
   # Process exprs
   if (verbose) autonomics.support::cmessage('\texprs')
   if (rm_reverse)              object %<>% autonomics.import::filter_features_("Reverse != '+'", verbose = verbose)
   if (rm_contaminants){        contaminant_var <- c('Potential contaminant', 'Contaminant') %>% intersect(autonomics.import::fvars(object))
   object %<>% autonomics.import::filter_features_(sprintf("`%s` != '+'", contaminant_var), verbose = verbose)
   }
   if (rm_complete_nondetects)  object %<>% autonomics.preprocess::filter_features_nonzero_in_some_sample(verbose = verbose)
   #object %>% autonomics.preprocess::qrilc_consistent_nondetects()
   if (log2transform){
      if (verbose) autonomics.support::cmessage('\t\tLog2 transform')
      autonomics.import::exprs(object) %<>% log2()
   }

   # Process sdata
   if (verbose) autonomics.support::cmessage('\tsdata')
   if (clean_snames)    object %<>% autonomics.import::clean_maxquant_snames(verbose = verbose)
   if (uniplex_snames)  object %<>% autonomics.import::uniplex_snames(verbose = verbose)
   object$subgroup <- object$sample_id %>% autonomics.import::guess_subgroup_values(verbose = verbose)
   #object$block    <- object$sample_id %>% autonomics.import::guess_subject_values( verbose = TRUE)

   # Process fdata
   if (verbose) autonomics.support::cmessage('\tfdata')
   if (verbose) autonomics.support::cmessage("\t\tRename: 'Majority protein IDs' -> 'Uniprot accessions'")
   autonomics.import::fvars(object) %<>% stringi::stri_replace_first_fixed('Majority protein IDs', 'Uniprot accessions')
   if (!is.null(fastafile)) object %>% autonomics.import::deconvolute_proteingroups(fastafile)

   # Return
   object
}
