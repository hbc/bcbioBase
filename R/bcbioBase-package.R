#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @importClassesFrom basejump tx2gene
#'
#' @importMethodsFrom basejump coerce
#' @importMethodsFrom S4Vectors coerce
#'
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_in_range assert_all_are_positive
#' @importFrom assertive.properties assert_has_no_duplicates assert_is_atomic
#'   assert_is_non_empty assert_is_scalar
#' @importFrom assertive.sets assert_are_disjoint_sets assert_is_subset
#' @importFrom assertive.strings assert_all_are_matching_regex
#'   assert_all_are_non_missing_nor_empty_character
#' @importFrom assertive.types assert_is_a_string assert_is_any_of
#'   assert_is_character assert_is_list assert_is_matrix assert_is_tbl_df
#' @importFrom basejump assertAllAreValidNames assertHasRownames
#'   assertIsAnImplicitInteger camel localOrRemoteFile makeNames
#'   printString readFileByExtension readYAML removeNA sanitizeNA
#' @importFrom Biostrings reverseComplement
#' @importFrom dplyr arrange everything funs group_by left_join mutate
#'   mutate_all mutate_at mutate_if select ungroup
#' @importFrom magrittr %>% set_colnames
#' @importFrom methods as is new
#' @importFrom plyr ldply
#' @importFrom rdrop2 drop_acc drop_auth drop_create drop_delete drop_exists
#'   drop_get_metadata drop_share drop_upload
#' @importFrom readr read_csv read_lines
#' @importFrom rlang !!! !! sym syms
#' @importFrom S4Vectors DataFrame tail
#' @importFrom stringr str_pad str_trunc
#' @importFrom tibble tibble
#' @importFrom tidyr expand
#' @importFrom utils globalVariables
"_PACKAGE"
