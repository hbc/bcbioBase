#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @import methods
"_PACKAGE"

#' @importFrom utils globalVariables
globalVariables(".")

metadataPriorityCols <- c("sampleID", "sampleName", "description")
