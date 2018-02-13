#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @import methods
#' @importFrom rlang !!! !! .data abort inform sym syms warn
#' @importFrom utils globalVariables
"_PACKAGE"



globalVariables(".")
metadataPriorityCols <- c("sampleID", "sampleName", "description")
