############################################################################
# min_occupancy_binning
############################################################################
#' @title min_occupancy_binning
#' @description A method to split a numeric verctor into bins with a minimum 
#' number of members.
#' @details This function is written as a helper function to 
#' \code{\link{calc_sig_b}}, where it serves to bin protein 
#' IDs into equal occupancy bins according to their summed peptide signal 
#' intensity.
#' 
#' \code{BinType} indicates how values should be distributed if \code{length(x)}
#' is not a multiple of \code{BinSize}. When using the MaxQuant/Perseus-analogous
#' \code{BinType="Hard"}, all but the last bin contain \code{BinSize} members, 
#' while \code{BinType="Optimal"} calculates the maximum number of bins with at 
#' least \code{BinSize} members and evenly distributes the members to those bins 
#' (see example below).
#' @param x A vector of numerical values tobe split into bins of ordered values.
#' @param BinSize An integer indicating the minimum number of members each bin
#' should have.
#' @param BinType Mode of bin positioning and size determination (see "Details";
#' defaulting to "Hard").
#' @return Returns a vector of integers that indicate how "NumVector" needs to be
#' split to achieve bins of consecutive values with the minimum size of "BinSize"
#' (see example).
#' @author Johannes Graumann
#' @keywords methods manip
#' @examples
#' # Example with length(x) as a multiple of BinSize
#' myVector <- c(1,9,2,8,3,7,4,6,5)
#' (splitMembership <- min_occupancy_binning(myVector,3))
#' split(myVector,splitMembership)
#' 
#' # Minimum occupancy not achievable with equal bin size
#' myVector <- c(1,9,2,8,3,7,4,6,5,10)
#' (splitMembership <- min_occupancy_binning(myVector,3))
#' split(myVector,splitMembership)
#' 
#' # "Hard" vs. "Optimal" BinType
#' table(min_occupancy_binning(seq(100),22,BinType="Hard"))
#' table(min_occupancy_binning(seq(100),22,BinType="Optimal"))
#' # Clean up
#' rm(myVector,splitMembership)
#' @noRd
min_occupancy_binning <- function(
  x,
  BinSize,
  BinType=c("Hard","Optimal")# Defaults to HARD (match.arg)
){
  # Input checks
  assertive.types::assert_is_numeric(x)
  assertive.types::assert_is_a_number(BinSize)
  assertive.numbers::assert_all_are_whole_numbers(BinSize)
  assertive.numbers::assert_all_are_positive(BinSize)
  BinType <- match.arg(BinType,c("Hard","Optimal"))
  # What is the number of bins fulfilling the minimum occupancy requirement?
  numberOfBins <- max(1,floor(length(x)/BinSize))
  if(numberOfBins < 2){
    warning(
      paste(
        "min_occupancy_binning: Length of 'x' less than 2 *",
        BinSize,
        "- returning one bin."
      )
    )
  }
  if(BinType == "Optimal"){
    quantiles <- stats::quantile(
      -x,seq.int(0,1,1/numberOfBins),na.rm=TRUE,names=FALSE
    )
    # Which quantile does any given value fall into?
    splitMembership <- cut(-x,quantiles,labels=FALSE,include.lowest=TRUE)
    # Make output consistent
    splitMembership <- as.numeric(splitMembership)
  } else { #BinType == "Hard"
    # Build the split-assignment
    splits <- rep(seq_len(numberOfBins),each=BinSize)
    nx <- length(x)
    nsplits <- length(splits)
    if(nsplits < nx){
      # Not enough splits; make last bin bigger
      splits <- c(splits,rep.int(numberOfBins,(nx - nsplits)))
    } else if (nsplits > nx){
      # Too many splits; make last bin smaller
      splits <- splits[seq_len(nx)]
    }
    # Order it by decreasing value
    orderByDecreasingValue <- order(x, decreasing = TRUE)
    # Re-create the original order and extract the split assigments
    orderByIndex <- order(seq_len(nx)[orderByDecreasingValue])
    splitMembership <- splits[orderByIndex]
  }
  # Return
  return(splitMembership)
}

######################################################
# calc_sig_a
######################################################
#' @title Calculate sig A
#' @description A method to calculate significance A as described in the original MaxQuant
#' publication (with modifications)
#' @details The numerical values fed to this function should represent an approx.
#' symetrical distribution and might require prior logarithmization.
#' 
#' Details of the calculation may be looked up in the supplementary material of
#' the citation listed (page 19f).
#' The algorithm differs from that original - emulating later MaxQuant/Perseus
#' versions - in multiplying p-values by a factor of two, thus converting the
#' originally one-sided test into a two-sided one (J. Cox, personal communication).
#' 
#' \code{force=FALSE} implies that the data fed to the algorithm is checked for 
#' prior logarithmization. The function throws an error if \code{min(x,na.rm=TRUE) >= 0}. 
#' This behavior may be overridden by setting \code{force=TRUE}.
#' 
#' This implementation shows good correlation with the one performed by
#' Perseus as of version 1.1.1.17. See the example below for a
#' demonstration.
#' 
#' Note that no multiple hypothesis testing is done.
#' @param x Vector of numerical values to calculate the outlier significance for 
#' (see "Details" for more).
#' @param force Boolean indicating whether or not to go ahead with likely 
#' non-logarithmized data (see "Details" for more).
#' @param method Method for adjusting p-values for multiple comparisons. See
#' \code{\link[stats]{p.adjust}}.
#' @return Returns a vector of p-Values.
#' @author Johannes Graumann
#' @references Cox, J., and Mann, M. (2008). MaxQuant enables high peptide identification
#' rates, individualized p.p.b.-range mass accuracies and proteome-wide protein
#' quantification. Nat. Biotechnol 26, 1367-1372.
#' This function is borrowed from pdapmain; see that pkg for examples.
#' @seealso \code{\link{calc_sig_b}}, \code{\link[stats]{p.adjust}}
#' @keywords distribution methods manip
#' @export
calc_sig_a <- function(
  x,
  force=FALSE,
  method="none"
){
  # Prerequisites checks
  assertive.types::assert_is_a_bool(force)
  assertive.types::assert_is_numeric(x)
  # We next calculate an outlier significance score for log protein
  # ratios (significance A).
  if(!force && min(stats::na.omit(sign(x)))>=0){
    stop("Ratio distributions should be logarithmized before calculating Significance A. Override with \'force=TRUE\'.") 
  }
  method <- match.arg(method, stats::p.adjust.methods)
  xIsPresent <- !is.na(x)
  
  # To make a robust and asymmetrical estimate of the standard deviation of the main distribution we calculate
  # the 15.87, 50, and 84.13 percentiles r-1, r0, and r1.
  # 15.87 and 84.13 are pnorm(-1) and pnorm(1) respectively, rounded to 2 dp.  That is, the points 
  # on the cumulative distribution that would be one standard deviation away from the mean, if it
  # were normally distributed.
  allR <- stats::quantile(x[xIsPresent], stats::pnorm(-1:1), names = FALSE)
  rNeg1 <- allR[1]
  r0 <- allR[2]
  r1 <- allR[3]
  
  # We define r1-r0 and r0-r-1 as the right- and left-sided robust standard deviations.
  # For a normal distribution, these would be equal to each other and to the conventional definition
  # of a standard deviation. A suitable measure for a ratio r>r0 of being significantly far away from
  # the main distribution would be the distance to r0 measured in terms of the right standard deviation
  # z=(r-r0)/(r1-r0)
  # Similarly, for a ratio r<r0 one would take
  # z=(r0-r)/(r0-r-1)
  # for the same purpose.
  z <- numeric(length(x)) #Default to no significance for missing values
  z[xIsPresent] <- (x[xIsPresent]-r0) / ifelse(x[xIsPresent] > r0, r1-r0, rNeg1-r0)
  
  # Under the null hypothesis of normal tails of the log ratio distribution the probability of obtaining
  # a value this large or larger is
  # significance A = 1/2erfc(z/sqrt(2)=1/sqrt(2*Pi) * Int (z,Inf){e^(-t^(2)/2) dt
  # The assumption that the ratio distribution is normal in a case where no differential
  # regulation takes place is reasonable, since the protein ratios are obtained as averages of
  # many SILAC peptide ratios. In the limit of a large number of SILAC peptides per protein
  # a normal distribution can be assumed due to the central limit theorem.
  significance <- stats::pnorm(z, lower.tail = FALSE)
  
  # In the original paper a one-sided test was used here.
  # Another multiplication with 2 takes care of that for a two sided test
  # (according to personal communication with J. Cox).
  twoSidedSignificance <- 2 * significance
  
  # For our dataset we report significant ratios based on significance B with a Benjamini-
  # Hochberg corrected p-value threshold of 0.05.
  stats::p.adjust(twoSidedSignificance, method)
}

###########################################################################################
# calc_sig_b
###########################################################################################
#' @title Calculate sig B
#' @description A method to calculate significance B as described in the original MaxQuant
#' publication (with modifications)
#' @details Significance B equates to \link[=calc_sig_a]{Significance A}
#' on data that has bin binned according to a binning vector "y" (summed peptide intensity 
#' in the original application) and a minimum bin size.
#' 
#' The numerical values fed to this function as "x" should represent an approx. symetrical 
#' distribution and might require prior logarithmization.
#' 
#' Details of the calculation may be looked up in the supplementary material of the 
#' citation listed (page 19f).
#' 
#' The algorithm differs from that original - emulating later MaxQuant/Perseus versions - 
#' by using a default bin size of 500 (J. Cox, personal communication).
#' 
#' \code{BinType} indicates how values should be distributed if \code{length(x)} is not a 
#' multiple of \code{BinSize}. When using the MaxQuant/Perseus-analogous 
#' \code{BinType="Hard"}, all but the last bin contain \code{BinSize} members, while 
#' \code{BinType="Optimal"} calculates the maximum number of bins with at least 
#' \code{BinSize} members and evenly distributes the members to those bins.
#' 
#' See the example below for a comparative correlation with the one performed by Perseus 
#' as of version 1.1.1.17.
#' 
#' Note that no multiple hypothesis testing is done.
#' @param x A vector of numerical values to calculate the outlier significance for. See 
#' "Details" for more.
#' @param y A vector of numerical values according to which the binning of values in 
#' \code{x} #' will be done. See "Details"for more.
#' @param BinSize The minimum number of values per bin from \code{x} binned according to 
#' \code{y}. Defaulting to 500.
#' @param BinType Mode of bin positioning and size determination (see "Details"; 
#' defaulting to "Hard").
#' @param ... Passed to \code{\link{calc_sig_a}}.
#' @return Returns a vector of p-Values.
#' @author Johannes Graumann
#' @references Cox, J., and Mann, M. (2008). MaxQuant enables high peptide identification
#' rates, individualized p.p.b.-range mass accuracies and proteome-wide protein
#' quantification. Nat. Biotechnol 26, 1367-1372.
#' This function is borrowed from pdapmain; see that pkg for examples.
#' @seealso \code{\link{calc_sig_a}}
#' @keywords distribution methods manip
#' @export
calc_sig_b <- function(
  x,y,
  BinSize=500,
  BinType=c("Hard","Optimal"),# Defaults to HARD (match.arg)
  ...){
  # Check prerequisites
  if(length(x) != length(y)){
    stop("'x' and 'y' are not of same length.")
  }
  assertive.types::assert_is_numeric(x)
  assertive.types::assert_is_numeric(y)
  assertive.types::assert_is_a_number(BinSize)
  assertive.numbers::assert_all_are_whole_numbers(BinSize)
  BinType <- match.arg(BinType,c("Hard","Optimal"))
  # Assemble a data.frame to work on
  tmpDF <- data.frame(MergeID=seq_along(x),x=x,y=y)
  # Extract entries that have all required entries
  subsetter <- !is.na(x) & !is.na(y)
  goodOnes <- tmpDF[subsetter,]
  # ...[[ it can be seen that the width of the bulk distribution of logarithmic
  # ratios depends on the protein intensity. For highly abundant proteins the
  # statistical spread of unregulated proteins is much more focused than for low
  # abundant ones. Because of this, a protein that shows, for instance, a ratio
  # of two should be very significant when it is highly abundant, while at very
  # low abundance it should only be marginally significant. To capture this
  # effect we define another quantity, called significance B, which is
  # calculated in the same way as significance A, but on protein subsets
  # obtained by grouping them into intensity bins. We divide the proteins into
  # bins of equal occupancy such that each bin contains at least 300 proteins.
  # The above calculation for significance A is then repeated in each bin to
  # obtain significance B.
  
  # How do the entries have to be subdivided into groups with optim. "BinSize" members each?
  splitMembership <- min_occupancy_binning(goodOnes[["y"]],BinSize,BinType=BinType)
  # Split the dataframe accordingly
  goodOnesBinned <- split(goodOnes,splitMembership)
  # Process each bin
  goodOnesBinned <- lapply(
    goodOnesBinned,
    function(z){
      significanceAOfBin <- calc_sig_a(z[["x"]], ...)
      assembly <- cbind(z,Significance.B=significanceAOfBin)
      return(assembly)
    }
  )
  # Reconcatenate
  goodOnes <- unsplit(goodOnesBinned,splitMembership)
  # Remerge with parias
  tmpDF <- merge(tmpDF,goodOnes,by=c("MergeID","x","y"),all=TRUE)
  # Set default value to 1.0
  tmpDF[is.na(tmpDF[["Significance.B"]]),"Significance.B"] <- 1.0
  # Return
  return(tmpDF[["Significance.B"]])
}

################
# Run "sig B" analysis and add p values to assayData
# @param pset ProtSet
# @importFrom magrittr %>%
# @importFrom autonomics.import exprs
# @importFrom autonomics.import intensities_num
# @importFrom autonomics.import intensities_den
# @export
# calc_sig_b_on_pset <- function(pset){
#   L2R  <- exprs(pset)
#   L2I <- (2^intensities_num(pset) + 2^intensities_den(pset)) %>% log2()
#   runSigB <- function(i, L2R, L10I){
#     suppressWarnings(calc_sig_b(L2R[, i], L10I[, i]))
#   }
#   # warning('sig_b currently returns p = 1 for Inf and -Inf ratios')
#   # TODO: make warning conditional on existence of Inf ratios, or move this to documentation
#   sig_b_mat <- sapply(1:ncol(pset), runSigB, L2R, L2I)
#   rownames(sig_b_mat) <- autonomics.import::fnames(pset)
#   colnames(sig_b_mat) <- autonomics.import::snames(pset)
#   sig_b_mat
# }

# @importFrom magrittr   %>%
# @noRd
#    add_sigb_to_fdata_single_contrast <- function(pset, contrast){
#      
#      # Get ratios
#      L2R <- autonomics.import::fdata(pset) %>% magrittr::extract2(sprintf('coef.%s', names(contrast)))
#      
#      # Get intensities
#      contrast_mat <- limma::makeContrasts(contrasts = contrast, 
#                                           levels    = create_design_matrix(pset))
#      selector <- autonomics.import::sdata(pset)$subgroup %in% rownames(contrast_mat)[contrast_mat != 0]
#      num_intensities <- 2^(autonomics.import::intensities_num(pset[, selector]))
#      den_intensities <- 2^(autonomics.import::intensities_den(pset[, selector]))
#      total_intensities <- rowSums(num_intensities + den_intensities)
#      L2I <- log2(total_intensities)
#      
#      # Run sigb
#      autonomics.import::fdata(pset)[[paste0('sigb.', names(contrast))]] <- suppressWarnings(calc_sig_b(L2R, L2I))
#    
#      # Return
#      return(pset)
# }

# Perform sigb outlier analysis for each contrast
# 
# @param  pset   ProtSet
# @param  contrasts  a character vector of contrasts (of subgroup levels)
# @return updated eset with limma results in fdata
#         added column: sigb.xxx \cr
# @author Aditya Bhagwat
# @importFrom magrittr %<>%
# @export
# add_sigb_to_fdata <- function(pset, contrasts){
#   for (i in seq_along(contrasts)){
#     pset %<>% add_sigb_to_fdata_single_contrast(contrasts[i])
#   }
#   return(pset)
# }