#' Identification of biologically interesting features with p-value-free FDR control
#'
#' \code{Clipper} identifies biologically interesting features with p-value-free FDR control on high-throughput biological data.
#'  \code{Clipper} is applicable to both differential analysis, including
#'  \itemize{
#'  \item differentially expressed genes identification from RNA-seq data
#'  \item differentially chromosomal interation regions identifiacation from Hi-C data,
#'  }
#'  and enrichment analysis, including
#'  \itemize{
#'  \item peak calling from ChIP-seq data
#'  \item peptide identification from mass spectrometry data.
#'  }
#'
#'
#' @param score.exp a numeric matrix of measurements under the experimental condition with rows being features
#' and columns being replicates.
#' @param score.back a numeric matrix of measurements under the background condition
#' with rows being features and columns being replicates. Its features should match features from \code{score.exp}.
#' @param FDR a numeric value or vector indicating the target FDR threshold(s).
#' @param analysis a character string specifying the analysis goal, must be either "differential" ("d")
#'  or "enrichment" ("e"). See below for details.
#' @param procedure a character string specifying the FDR control procedure, must be "GZ", "BC" or "aBH". If not
#' specified, \code{Clipper} uses the GZ procedure when \code{analysis = "differential"} and the BC procedure when \code{analysis = "enrichment"}.
#'  See below for details.
#' @param contrast.score a character string specifying the contrast score, must be "max" or "diff". If not specified,
#'  \code{Clipper} uses the maximum (max) contrast score when \code{procedure = "GZ"} and the difference (diff) contrast score when \code{procedure = "BC"}.
#'
#' @param n.permutation an integer specifying the number of permutations.
#' Effective only when \code{procedure = "GZ"}. See below for details.
#' @param seed random seed, used in permutations.
#'
#' @details
#' \code{Clipper} aims to contrast two conditions to
#' reliably screen interesting features. The interesting means "differential" or "enriched". Differential
#' features are defined as those that have different expected measurements (without measurement errors) between
#' two conditions, and the detection of such differential features is called differential analysis. In contrast,
#' enriched features are defined as those
#'  that have higher expected measurements under the experimental/treatment condition than the background
#'  condition, i.e., the negative control. The detection of such enriched features is called enrichment analysis.
#'
#'
#' The specification of the \code{procedure} indicates the type of FDR control procedure implemented.
#' In our paper, we introduced three FDR control procedures: GZ, BC and aBH. The GZ procedure applies to both
#' differential and enrichment analysis. It requires total number of replicates under two conditions to be at least 3.
#' The BC procedure only applies to enrichment analysis with equal numbers of replicates. Both GZ and BC have
#' theoretical guarantee on FDR control. aBH only applies to enrichment analysis and does not have theoretical
#' FDR control. However, empirical evidence shows that even though aBH may slightly exceed the target FDR threshold,
#' it is more powerful than BC. Therefore, if BC makes no discovery, we recommend users to use aBH at
#' the current FDR threshold or better yet, to increase the FDR threshold.
#'
#'
#'
#'
#' @return \code{Clipper} returns a list containing the following components:
#' \describe{
#'   \item{\code{contrast.score:}}{the type of contrast scores, 'max' or 'diff';}
#'   \item{\code{contrast.score.value:}}{the values of contrast scores;}
#'   \item{\code{FDR:}}{a vector of the target FDR threshold(s);}
#'   \item{\code{contrast.score.thre:}}{a vector of the threshold(s) on contrast scores corresponding to the FDR threshold(s);}
#'   \item{\code{discovery:}}{a list of identified discoveries. Each component contains discovered
#'   features, coded as the row indices of \code{score.exp} and \code{score.back}, at a FDR threshold.}
#' }
#'
#' @export
#' @import parallel
#' @importFrom stats complete.cases dgamma dnorm prcomp quantile rgamma rnorm sd uniroot
#'
#' @references
#' @author Xinzhou Ge, \email{xinzhouge@ucla.edu}
#' @author Yiling Chen, \email{yiling0210@ucla.edu}
#' @author Jingyi Jessica Li, \email{jli@stat.ucla.edu}
#'
#' @examples
#' ## differential analysis
#' re1 <- Clipper(exp_d, back_d, analysis = "differential", FDR = c(0.01, 0.05, 0.1))
#' re1$discoveries
#' ## enrichment analysis with multiple
#' re2 <- Clipper(exp_e, back_e, analysis = "enrichment", FDR = c(0.01, 0.05, 0.1))
#' re2$discoveries
#'
Clipper <- function(score.exp, score.back, analysis, FDR = 0.05,
                    procedure = NULL, contrast.score = NULL,
         n.permutation = NULL, seed = 12345){
  analysis = match.arg(analysis, choices = c('differential', 'enrichment'), several.ok = F)
  if (analysis == 'differential') {
    if (is.null(contrast.score)) {
      contrast.score <- 'max'
    }else{
      contrast.score <- match.arg(contrast.score, choices = c('diff', 'max'), several.ok = F)
    }
    if (is.null(procedure)) {
      procedure <- 'GZ'
    }else{
      procedure = match.arg(procedure, choices = c('BC','aBH','GZ'), several.ok = F)
    }
    re <- clipper2sided(score_exp = score.exp, score_back = score.back, FDR = FDR,
                        nknockoff = n.permutation,
                        contrastScore_method = contrast.score, importanceScore_method = 'diff',
                        FDR_control_method = procedure, ifpowerful = F, seed = seed)
    FDR_nodisc = sapply(re$results, function(re_i){
      length(re_i$discovery) == 0
    })
    if( any(FDR_nodisc & contrast.score == 'max') ){
      warning(paste0('At FDR = ', paste0(FDR[FDR_nodisc], collapse = ', '), ', no discovery has been found using max contrast score. To make more discoveries, switch to diff contrast score or increase the FDR threshold. '))
    }
  }
  if (analysis == 'enrichment') {
    if (is.null(procedure)) {
      procedure <- 'BC'
    }else{
      procedure = match.arg(procedure, choices = c('BC','aBH','GZ'), several.ok = F)
    }
    if (ncol(score.exp)!=ncol(score.back)){
      procedure <- 'GZ'
    }
    if (is.null(contrast.score)) {
      if (procedure=='BC'){
        contrast.score <- 'diff'
      }
      if (procedure=='GZ'){
        contrast.score <- 'max'
      }
    }else{
      contrast.score <- match.arg(contrast.score, choices = c('diff', 'max'), several.ok = F)
    }
    if (procedure=='BC'){
      re <- clipper1sided(score_exp = score.exp, score_back = score.back, FDR = FDR,
                          nknockoff = n.permutation,
                          importanceScore_method = contrast.score,
                          FDR_control_method = procedure, ifpowerful = F, seed = seed)
    }
    if (procedure=='GZ'){
      re <- clipper1sided(score_exp = score.exp, score_back = score.back, FDR = FDR,
                          nknockoff = n.permutation,
                          contrastScore_method = contrast.score,
                          FDR_control_method = procedure, ifpowerful = F, seed = seed)
    }
    FDR_nodisc = sapply(re$results, function(re_i){
      length(re_i$discovery) == 0
    })
    if( any(FDR_nodisc & procedure != 'aBH') ){
      warning(paste0('At FDR = ', paste0(FDR[FDR_nodisc], collapse = ', '), ', no discovery has been found using current procedure. To make more discoveries, switch to aBH procedure or increase the FDR threshold.'))
    }
  }

  contrast.score.value <- re$contrastScore
  thre <- unlist(lapply(re$results,'[[','thre'))
  discoveries <- lapply(re$results,'[[','discovery')
  re <- list(contrast.score = contrast.score,
             contrast.score.value = contrast.score.value,
             FDR = FDR,
             contrast.score.thre = thre,
             discoveries = discoveries)
  return(re)

}
