#' @importMethodsFrom S4Vectors coerce
#' @importMethodsFrom basejump coerce
#'
#' @importFrom Biostrings reverseComplement
#' @importFrom S4Vectors DataFrame metadata<- na.omit tail
#' @importFrom basejump Tx2Gene camel import localOrRemoteFile makeNames
#'   makeSampleData printString realpath removeNA sanitizeNA
#' @importFrom dplyr arrange everything group_by left_join mutate mutate_all
#'   mutate_at mutate_if select ungroup
#' @importFrom goalie allAreAtomic allAreFiles allAreMatchingRegex
#'   areDisjointSets assert containsAURL false hasLength hasNoDuplicates
#'   hasRownames isADirectory isAFile isAny isCharacter isFile isMatchingRegex
#'   isInRange isInt isNonEmpty isNonNegative isScalar isString isSubset
#'   validNames
#' @importFrom magrittr %>% set_colnames
#' @importFrom methods as is new
#' @importFrom plyr ldply
#' @importFrom rdrop2 drop_acc drop_auth drop_create drop_delete drop_exists
#'   drop_get_metadata drop_share drop_upload
#' @importFrom readr read_csv read_lines
#' @importFrom rlang !!! !! sym syms
#' @importFrom stringr str_match str_pad str_replace str_trunc
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr expand
#' @importFrom utils globalVariables
NULL
