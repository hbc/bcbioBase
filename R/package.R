#' bcbioBase
#'
#' Base functions for bcbio R packages.
#'
#' @keywords internal
"_PACKAGE"



## S4 generics and methods =====================================================

#' @importFrom AcidExperiment makeSampleData
#' @importFrom AcidPlyr rbindToDataFrame
#' @importFrom S4Vectors metadata<- na.omit tail
#' @importFrom pipette import removeNA sanitizeNA
#' @importFrom syntactic camelCase makeNames
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase printString realpath strMatch
#' @importFrom AcidCLI abort alert alertInfo alertWarning dl toInlineString ul
#' @importFrom S4Vectors DataFrame
#' @importFrom goalie allAreAtomic allAreFiles areDisjointSets
#' allAreMatchingRegex assert hasLength hasRownames isADirectory isAFile
#' isCharacter isMatchingRegex isInRange isInt isScalar isString isSubset
#' validNames
#' @importFrom methods as is new
#' @importFrom utils packageName packageVersion
NULL
