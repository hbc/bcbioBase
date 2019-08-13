#' @importMethodsFrom S4Vectors coerce
#' @importMethodsFrom basejump coerce
#'
#' @importFrom S4Vectors DataFrame metadata<- na.omit tail
#' @importFrom basejump camelCase import localOrRemoteFile makeNames
#'   makeSampleData printString realpath removeNA sanitizeNA
#' @importFrom dplyr arrange mutate mutate_all mutate_if
#' @importFrom goalie allAreAtomic allAreFiles areDisjointSets assert
#'   hasRownames isADirectory isAFile isCharacter isMatchingRegex isInRange
#'   isInt isNonEmpty isScalar isString isSubset validNames
#' @importFrom magrittr %>% set_colnames
#' @importFrom methods as is new
#' @importFrom plyr ldply
#' @importFrom rdrop2 drop_acc drop_auth drop_create drop_delete drop_exists
#'   drop_get_metadata drop_share drop_upload
#' @importFrom readr read_csv
#' @importFrom rlang !!! !! sym syms
#' @importFrom stringr str_match str_replace
#' @importFrom tibble as_tibble tibble
#' @importFrom utils globalVariables
NULL
