#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_in_range assert_all_are_positive
#' @importFrom assertive.properties assert_has_dimnames assert_has_no_duplicates
#'   assert_is_non_empty
#' @importFrom assertive.sets assert_are_disjoint_sets assert_is_subset
#' @importFrom assertive.strings assert_all_are_matching_regex
#'   assert_all_are_non_missing_nor_empty_character
#' @importFrom assertive.types assert_is_a_string assert_is_an_integer
#'   assert_is_any_of assert_is_character assert_is_data.frame
#'   assert_is_function assert_is_list
#' @importFrom basejump assertAllAreValidNames assertIsAnImplicitInteger
#'   assertIsTx2gene camel fixNA localOrRemoteFile makeNames printString
#'   readFileByExtension readYAML removeNA
#' @importFrom Biostrings reverseComplement
#' @importFrom dplyr arrange everything funs group_by left_join mutate
#'   mutate_all mutate_at mutate_if select ungroup
#' @importFrom ggplot2 aes geom_hline geom_label geom_vline
#' @importFrom ggrepel geom_label_repel
#' @importFrom grid arrow unit
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom methods as is
#' @importFrom plyr ldply
#' @importFrom rdrop2 drop_acc drop_auth drop_create drop_delete drop_exists
#'   drop_get_metadata drop_share drop_upload
#' @importFrom readr read_csv read_lines
#' @importFrom rlang !!! !! sym syms
#' @importFrom S4Vectors aggregate as.data.frame merge tail
#' @importFrom stats as.formula
#' @importFrom stringr str_pad str_trunc
#' @importFrom tibble column_to_rownames rownames_to_column tibble
#' @importFrom tidyr expand
#' @importFrom utils globalVariables
"_PACKAGE"
