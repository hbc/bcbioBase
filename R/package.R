## FIXME RETHINK THE LDPLY APPROACH HERE.



#' bcbioBase
#'
#' Base functions for bcbio R packages.
#'
#' @keywords internal
#'
#' @importMethodsFrom S4Vectors coerce
#' @importMethodsFrom basejump coerce
#'
#' @importFrom S4Vectors DataFrame metadata<- na.omit tail
#' @importFrom basejump camelCase import localOrRemoteFile makeNames
#'   makeSampleData printString realpath removeNA sanitizeNA
#' @importFrom cli cli_alert cli_alert_info cli_alert_warning cli_div cli_dl
#'   cli_text cli_ul
#' @importFrom goalie allAreAtomic allAreFiles areDisjointSets
#'   allAreMatchingRegex assert hasLength hasRownames isADirectory isAFile
#'   isCharacter isMatchingRegex isInRange isInt isScalar isString isSubset
#'   validNames
#' @importFrom methods as is new
#' @importFrom plyr ldply
#' @importFrom stringr str_match
#' @importFrom utils packageName packageVersion
"_PACKAGE"
