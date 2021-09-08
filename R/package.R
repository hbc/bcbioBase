#' bcbioBase
#'
#' Base functions for bcbio R packages.
#'
#' @keywords internal
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom basejump DataFrame abort alert alertInfo alertWarning camelCase
#'   dl import localOrRemoteFile makeNames makeSampleData metadata<- na.omit
#'   packageName packageVersion printString rbindToDataFrame realpath removeNA
#'   sanitizeNA str_match tail toInlineString ul
#' @importFrom goalie allAreAtomic allAreFiles areDisjointSets
#'   allAreMatchingRegex assert hasLength hasRownames isADirectory isAFile
#'   isCharacter isMatchingRegex isInRange isInt isScalar isString isSubset
#'   validNames
#' @importFrom methods as is new
"_PACKAGE"
