################################################################################
#                                                                              #
#                     evidence  ->  data.table                                 #
#                                   matrix                                     #
#                                                                              #
################################################################################

#' Compute geometric median
#' @param x numeric vector for which to compute geometric median
#' @export
geomedian <- function(x){
   exp(stats::median(log(x), na.rm = TRUE))
}

#   1) rm peptides with variable modifications (phospho and deamidation)
#      as these peptides should not be used for quantification purposes
#      See MQ google group for more details.
#   2) rm unquantified 'MSMS' peptides
#      These are the peptides for which no MS1 label pair was detected, so no quantification available
#      Source: https://groups.google.com/d/msg/maxquant-list/ay4hItVFGMc/b8jxCoKGCAAJ
#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %<>%
filter_evidences <- function(evidence){

   # Satisfy CHECK
   Reverse <- Contaminant <- SeqNeverDeamidated <- `Deamidation (N)` <- NULL
   Sequence <- SeqNeverPhosphorylated <- `Phospho (STY)` <- NULL

   autonomics.support::cmessage('\t%d features in total', nrow(evidence))

   n <- nrow(evidence)
   evidence %<>% magrittr::extract(Reverse !=  '+'   )
   autonomics.support::cmessage('\t%d after removing reverse proteins', nrow(evidence), n)

   n <- nrow(evidence)
   evidence %<>% magrittr::extract(Contaminant !=  '+'   )
   autonomics.support::cmessage('\t%d after removing contaminants', nrow(evidence), n)

   if ('Deamidation (N)' %in% names(evidence)){
      n <- nrow(evidence)
      evidence[, SeqNeverDeamidated := all(`Deamidation (N)` == 0), by = list(Sequence)]#, Experiment, Labels)]
      evidence <- evidence[SeqNeverDeamidated==TRUE]
      autonomics.support::cmessage('\t%d from sequences that are not deamidated in the same sample', nrow(evidence), n)
   }
   if ('Phospho (STY)' %in% names(evidence)){
      n <- nrow(evidence)
      evidence[, SeqNeverPhosphorylated := all(`Phospho (STY)`   == 0), by = list(Sequence)]#, Experiment, Labels)]
      # by=Sequence gave better results than by=list(Sequence, Experiment, Labels) for Billing differentiation
      evidence <- evidence[SeqNeverPhosphorylated==TRUE]
      autonomics.support::cmessage('\t%d from sequences that are not phosphorylated in the same sample', nrow(evidence), n)
   }
   evidence

   return(evidence)
}

#' Read evidence file
#'
#' Read (relevant variables of) evidence file into a data table
#' @param evidence_file path to evidence file
#' @examples
#' require(magrittr)
#' if (require(billing.differentiation.data)){
#'    system.file('extdata/maxquant/evidence.txt',
#'                 package = 'billing.differentiation.data') %>%
#'    read_evidences()
#'
#' }
#' @importFrom magrittr   %<>%
#' @export
read_evidences <- function(evidence_file){

   header <- data.table::fread(evidence_file, nrows = 0, header = TRUE)
   sel_cols <- c('id', 'Protein group IDs', 'Experiment', 'Type', 'Peptide ID', 'Score',
                 'Intensity L', 'Intensity M', 'Intensity H',
                 'Ratio H/L normalized', 'Ratio M/L normalized', 'Ratio H/M normalized',
                 'Contaminant', 'Potential contaminant',  'Reverse', 'Phospho (STY)', 'Deamidation (N)', 'Sequence',
                 'Phospho (STY) site IDs')
   sel_cols %<>% magrittr::extract(sel_cols %in% names(header))

   suppressWarnings(
      evidence <- data.table::fread(
         evidence_file,
         select = sel_cols
      ))
   if ('Potential contaminant' %in% names(evidence)){
      data.table::setnames(evidence, 'Potential contaminant', 'Contaminant')
   }
   evidence
}

#' @importFrom data.table   data.table   :=
#' @importFrom magrittr     %>%   %<>%
reformat_evidence_intensities_n_ratios <- function(evidence){

   # Satisfy CHECK
   HH <- `Intensity H` <- MM <- `Intensity M` <- LL <- `Intensity L` <- NULL
   Intensity.HL <- Intensity.ML <- Intensity.HM <- NULL

   # Compute summed intensities
   # Using workaround for "evidence[, HL:= H + L]", which gives errors when melting (data.table bug?)
   if ('Intensity H' %in% names(evidence)){   evidence[, HH:=`Intensity H`]   }
   if ('Intensity M' %in% names(evidence)){   evidence[, MM:=`Intensity M`]   }
   if ('Intensity L' %in% names(evidence)){   evidence[, LL:=`Intensity L`]   }
   if ('Ratio H/L normalized' %in% names(evidence)){
      names(evidence) %<>% stringi::stri_replace_first_fixed('Intensity H', 'Intensity.HL')
      evidence[, Intensity.HL := HH + LL]
   }
   if ('Ratio M/L normalized' %in% names(evidence)){
      names(evidence) %<>% stringi::stri_replace_first_fixed('Intensity M', 'Intensity.ML')
      evidence[, Intensity.ML := MM + LL]
   }
   if ('Ratio H/M normalized' %in% names(evidence)){
      names(evidence) %<>% stringi::stri_replace_first_fixed('Intensity L', 'Intensity.HM')
      evidence[, Intensity.HM := HH + MM]
   }
   suppressWarnings(evidence[, HH:=NULL]) # Suppress warning if HH does not exist
   suppressWarnings(evidence[, MM:=NULL])
   suppressWarnings(evidence[, LL:=NULL])

   # Reformat ratios
   names(evidence) %<>% stringi::stri_replace_first_regex('Ratio ([HM])/([LM]) normalized', 'Ratio.$1$2')

   # Select relevant columns
   sel_cols <- c('id', 'Protein group IDs', 'Experiment', 'Type', 'Peptide ID', 'razor',
                 'Ratio.HL', 'Ratio.ML', 'Ratio.HM',
                 'Intensity.HL', 'Intensity.ML', 'Intensity.HM',
                 'Phospho (STY)', 'Deamidation (N)', 'Sequence', 'Score', 'Reverse', 'Contaminant')
   sel_cols <- sel_cols[sel_cols %in% names(evidence)]
   evidence <- evidence[, sel_cols, with = FALSE]

   # Return
   return(evidence)
}

#' @importFrom data.table   data.table
#' @importFrom magrittr     %<>%         %>%
melt_evidences <- function(evidence){

   # Satisfy CHECK
   . <- NULL

   evidence %<>% reformat_evidence_intensities_n_ratios()

   labels <- c('HL', 'ML', 'HM') %>% magrittr::extract(sprintf('Ratio.%s', .) %in% names(evidence))
   ratio_names <- sprintf('Ratio.%s', labels)
   intensity_names <- sprintf('Intensity.%s', labels)

   evidence %<>% data.table::melt.data.table(measure.vars = list(ratio_names, intensity_names), value.name = c('Ratio', 'Intensity'))
   levels(evidence$variable) <- labels
   names(evidence) %<>% stringi::stri_replace_first_fixed('variable', 'Labels')
   evidence
}

#' @importFrom data.table   data.table   .N   :=
#' @importFrom magrittr     %>%
summarize_evidences_per_proteingroup <- function(evidence){

   # Satisfy CHECK
   . <- Ratio <- nobs <- `Protein group IDs` <- Experiment <- NULL
   ndirectobs <- Type <- ndirectpeptides <- `Peptide ID` <- NULL
   pcor <- Intensity <- Labels <- Unique <- NULL

   # Rm NA ratios and record no of obs
   evidence <- evidence[!is.na(Ratio)]
   evidence[, nobs           := .N,                                                              by = list(`Protein group IDs`, Experiment, Labels)]
   evidence[, ndirectobs     := sum(!Type %in% c('MSMS', 'ISO-MSMS')),                           by = list(`Protein group IDs`, Experiment, Labels)]
   evidence[, ndirectpeptides:= length(unique(`Peptide ID`[!Type %in% c('MSMS', 'ISO-MSMS')])),  by = list(`Protein group IDs`, Experiment, Labels)]

   # 1-3 obs: geomedian of all obs
   autonomics.support::cmessage('\t1-3 obs: median of all obs')
   median_of_max_three_obs <- evidence[nobs <= 3,
                                       list(Ratio = geomedian(Ratio), Type = 'median_of_max_three_obs',
                                            Unique = all(!Unique)),
                                       by = list(`Protein group IDs`, Experiment, Labels)];
   evidence <- evidence[nobs>3]

   # 3+ direct peptides: geomedian of direct obs
   autonomics.support::cmessage('\t3+ direct peptides: median of direct obs')
   median_of_direct_obs <- evidence[ndirectpeptides >=3,
                                    list(Ratio = geomedian(Ratio[!Type %in% c('MSMS', 'ISO-MSMS')]),
                                         Type  = 'median_of_direct_obs',
                                         Unique = all(Unique)),
                                    by = list(`Protein group IDs`, Experiment, Labels)]
   evidence <- evidence[ndirectpeptides<3]

   # cor(LR,LI)<0.005: max(lm(LR,LI))  # suppress unimportant warning by cor.test
   autonomics.support::cmessage('\tp.cor(LR,LI)<0.005: take max of lin fit')
   suppressWarnings(evidence[, pcor:= stats::cor.test(log(Ratio), log(Intensity), method = 'spearman')$p.value,
                             by = list(`Protein group IDs`, Experiment, Labels)])
   max_of_lin_fit <- evidence[pcor < 0.005, list(Ratio = stats::predict(stats::lm(log(Ratio) ~ log(Intensity))) %>% magrittr::extract2(which.max(abs(.))) %>% exp(),
                                                 Type  = 'max_of_lin_fit',
                                                 Unique = all(Unique)),
                              by = list(`Protein group IDs`, Experiment, Labels)]
   evidence <- evidence[pcor>=0.005]

   # Remaining: take geomedian of all observations
   autonomics.support::cmessage('\tremaining: geomedian of all obs')
   median_of_all_obs <- evidence[, list(Ratio = geomedian(Ratio), Type = 'median_of_all_obs', Unique = all(Unique)),
                                 by = list(`Protein group IDs`, Experiment, Labels)]
   # rbind and return
   data.table::rbindlist(list(median_of_max_three_obs, median_of_direct_obs, median_of_all_obs, max_of_lin_fit))
}

#' @importFrom data.table   data.table   :=
melt_razor_evidences <- function(evidence){

   # Satisfy CHECK
   `Protein group IDs` <- id <- Labels <- Unique <- NULL

   unique_evidences <- evidence[!stringi::stri_detect_fixed(`Protein group IDs`, ';')]
   razor_evidences  <- evidence[stringi::stri_detect_fixed(`Protein group IDs`, ';')]
   deconvolution <- razor_evidences[, list(`Protein group IDs` = unlist(strsplit(`Protein group IDs`, ';'))), by = list(id, Labels)]
   razor_evidences[, `Protein group IDs`:=NULL]
   data.table::setkey(razor_evidences, id, Labels)
   data.table::setkey(deconvolution, id, Labels)
   razor_evidences <- merge(razor_evidences, deconvolution)
   razor_evidences[, Unique:=FALSE]
   unique_evidences[, Unique:=TRUE]
   data.table::rbindlist(list(unique_evidences, razor_evidences), use.names = TRUE)
}

#' @importFrom data.table   data.table   :=   .N
#' @importFrom magrittr     %<>%
assign_razor_evidences <- function(evidence){

   # Satisfy CHECK
   NUniquePeptides <- `Peptide ID` <- `Protein group IDs` <- NPeptides <- NULL
   NUniqueEvidences <- NEvidences <- MedianUniqueScore <- Score <- MedianScore <- NULL
   Unique <- id <- Labels <- NULL

   # Count statistics per PG
   evidence[, NUniquePeptides   := length(unique(`Peptide ID`[Unique==TRUE])),        by = list(`Protein group IDs`)]#, Experiment, Labels)]
   evidence[, NPeptides         := length(unique(`Peptide ID`)),                      by = list(`Protein group IDs`)]#, Experiment, Labels)]
   evidence[, NUniqueEvidences  := sum(Unique==TRUE),                                 by = list(`Protein group IDs`)]#, Experiment, Labels)]
   evidence[, NEvidences        := .N,                                                by = list(`Protein group IDs`)]#, Experiment, Labels)]
   evidence[, MedianUniqueScore := stats::median(Score[Unique==TRUE], na.rm = TRUE),  by = list(`Protein group IDs`)]#, Experiment, Labels)]
   evidence[, MedianScore       := stats::median(Score, na.rm = TRUE),                by = list(`Protein group IDs`)]#, Experiment, Labels)]

   # Select
   unique_evidences <- evidence[Unique==TRUE]
   razor_evidences  <- evidence[Unique==FALSE]
   # razor_evidences %<>% magrittr::extract(, Filter := NUniquePeptides   >  0                     , by = list(id, Labels)) %>% magrittr::extract(Filter==TRUE)
   razor_evidences %<>% magrittr::extract(, Filter := NUniquePeptides   == max(NUniquePeptides  ), by = list(id, Labels)) %>% magrittr::extract(Filter==TRUE)
   razor_evidences %<>% magrittr::extract(, Filter := NUniqueEvidences  == max(NUniqueEvidences ), by = list(id, Labels)) %>% magrittr::extract(Filter==TRUE)
   razor_evidences %<>% magrittr::extract(, Filter := MedianUniqueScore == max(MedianUniqueScore), by = list(id, Labels)) %>% magrittr::extract(Filter==TRUE)
   razor_evidences %<>% magrittr::extract(, Filter := NPeptides         == max(NPeptides        ), by = list(id, Labels)) %>% magrittr::extract(Filter==TRUE)
   razor_evidences %<>% magrittr::extract(, Filter := NEvidences        == max(NEvidences       ), by = list(id, Labels)) %>% magrittr::extract(Filter==TRUE)

   # Merge back
   razor_evidences <- razor_evidences[, names(unique_evidences), with = FALSE]
   data.table::rbindlist(list(unique_evidences, razor_evidences))
}


#' Aggregate evidences into protein groups
#'
#' Reads evidences, filters, deconvolutes razor sequences, melts labels, and summarizes per protein group.
#'
#' @details
#' Reads the evidence.txt file. \cr \cr
#' Then filters features by removing:
#' \itemize{
#'    \item sequences with variable modifications
#'    \item contaminant and reverse sequences
#'    \item features with no quantifications
#' }
#' Deconvolutes razor sequences by assigning them to the protein group with most other sequences. \cr
#' Melts the (labels of the) evidence table. \cr
#' Summarizes eviences per protein group:
#' \itemize{
#'    \item 1-3 obs: geomedian of all obs
#'    \item 3+ direct peptides: geomedian of direct obs
#'    \item cor(LR,LI)<0.005: max(lm(LR,LI))
#'    \item remaining obs: geomedian of all obs
#' }
#'
#' @param evidence_file path to evidence file
#' @examples
#' require(magrittr)
#' \dontrun{
#' if (require(billing.differentiation.data)){
#'    evidence_file <- system.file('extdata/maxquant/evidence.txt',
#'                                  package = 'billing.differentiation.data')
#'    evidences_to_proteingroups(evidence_file)
#' }
#' }
#' @importFrom data.table   data.table
#' @importFrom magrittr     %<>%   %>%
#' @export
evidences_to_proteingroups <- function(evidence_file){

   # Satisfy CHECK
   Ratio <- NULL

   autonomics.support::cmessage('Read evidence file')
   evidence <- read_evidences(evidence_file)

   autonomics.support::cmessage('Melt labels')
   evidence %<>% melt_evidences()
   evidence <- evidence[!is.na(Ratio)]

   autonomics.support::cmessage('Filter features')
   evidence %<>% filter_evidences()

   autonomics.support::cmessage('Deconvolute razor peptides')
   evidence %<>% melt_razor_evidences()
   evidence %<>% assign_razor_evidences()

   autonomics.support::cmessage('Summarize per protein group')
   protein_groups <- evidence %>% summarize_evidences_per_proteingroup()
   protein_groups
   #autonomics.support:cmessage('Cast value matrix')
   #ratios      <- protein_groups %>% cast_value_matrix('Ratio',     sample_file)
   #types       <- protein_groups %>% cast_value_matrix('Type',      sample_file)

   # autonomics.support::cmessage('Summarize ratios (per protein group and experiment)')
   # normratios <- evidence %>% summarize_evidences_per_proteingroup()

   # Return
   #return(ratios)

}

#' @importFrom magrittr %>%
robustify_protein_group_ids <- function(pg_ids){

   # Satisfy CHECK
   . <- NULL

   pg_ids %>%
      as.numeric() %>%
      formatC(width = ceiling(max(log10(.))), flag = '0', format = 'd') %>%
      paste0('PG', .)
}


prettify_sample_names <- function(value_matrix, sample_file){
   sample_design <- load_maxquant_design(sample_file)
   idx <- match(colnames(value_matrix),
                paste0(sample_design$experiment, '_',
                       stringi::stri_replace_first_fixed(sample_design$labels, '_', ''))
   )
   colnames(value_matrix) <- sample_design$sample_id[idx]
   value_matrix
}


#' Matrixify summarized evidences
#'
#' Convert summarized evidences from long data.table into
#' feature x sample matrix
#'
#' @param summarized_evidences data.table with summarized evidences
#' @param sample_file          file with sample annotation
#' @param var                  variable to extract and matrixify
#' @examples
#' require(magrittr)
#' \dontrun{
#' if (require(billing.differentiation.data)){
#'    evidence_file <- system.file(
#'       'extdata/maxquant/evidence.txt',
#'       package = 'billing.differentiation.data')
#'    sample_file <- system.file(
#'       'extdata/maxquant/sample_design.txt',
#'       package = 'billing.differentiation.data')
#'    summarized_evidences <- evidences_to_proteingroups(evidence_file)
#'    summarized_evidences %>% matrixify_summarized_evidences(sample_file)
#' }
#' }
#' @importFrom data.table   data.table
#' @importFrom magrittr     %>%
#' @export
matrixify_summarized_evidences <- function(summarized_evidences, sample_file, var = 'Ratio'){

   # Satisfy CHECK
   `Protein group IDs` <- NULL

   # Cast into (feature x sample) format
   value_df <- summarized_evidences %>% data.table::dcast(`Protein group IDs` ~ Experiment + Labels, value.var = c(var))
   protein_group_ids <- value_df[, `Protein group IDs`] %>% robustify_protein_group_ids()

   # Matrixify
   value_df %<>% magrittr::extract(, -1, with = FALSE)
   if        (var == 'Ratio'){
      value_mat <- value_df %>% data.matrix()
   } else if (var == 'Type'){
      value_mat <- value_df %>% as.matrix()
   }

   # Add sample and feature names
   value_mat %>% magrittr::set_rownames(protein_group_ids) %>%
      prettify_sample_names(sample_file)
}

