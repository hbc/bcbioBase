#' bcbioBase
#'
#' Base functions for bcbio R packages.
#'
#' @keywords internal
#'
#' @importMethodsFrom S4Vectors coerce
#' @importMethodsFrom basejump coerce
#'
#' @importFrom Biostrings reverseComplement
#' @importFrom S4Vectors DataFrame metadata<- na.omit tail
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_in_range
#'   assert_all_are_non_negative
#' @importFrom assertive.properties assert_has_no_duplicates assert_is_atomic
#'   assert_is_non_empty assert_is_scalar
#' @importFrom assertive.sets assert_are_disjoint_sets assert_is_subset
#' @importFrom assertive.strings assert_all_are_matching_regex
#'   assert_all_are_non_missing_nor_empty_character
#'   assert_any_are_matching_regex
#' @importFrom assertive.types assert_is_a_string assert_is_all_of
#'   assert_is_an_integer assert_is_any_of assert_is_character assert_is_list
#'   assert_is_matrix assert_is_tbl_df
#' @importFrom assertthat assert_that
#' @importFrom basejump Tx2Gene camel localOrRemoteFile makeNames makeSampleData
#'   printString realpath removeNA sanitizeNA
#' @importFrom dplyr arrange everything funs group_by left_join mutate
#'   mutate_all mutate_at mutate_if select ungroup
#' @importFrom goalie assertAreValidNames assertHasRownames
#'   assertIsImplicitInteger
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
"_PACKAGE"
