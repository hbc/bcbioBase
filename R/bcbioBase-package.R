#' bcbioBase
#'
#' Base functions and generics for bcbio R packages.
#'
#' @import methods
#' @importFrom assertive assert_all_are_dirs
#' @importFrom assertive assert_all_are_existing_files
#' @importFrom assertive assert_all_are_in_range
#' @importFrom assertive assert_all_are_non_missing_nor_empty_character
#' @importFrom assertive assert_are_disjoint_sets
#' @importFrom assertive assert_are_identical
#' @importFrom assertive assert_are_intersecting_sets
#' @importFrom assertive assert_has_colnames
#' @importFrom assertive assert_has_dimnames
#' @importFrom assertive assert_has_names
#' @importFrom assertive assert_has_rownames
#' @importFrom assertive assert_has_no_duplicates
#' @importFrom assertive assert_is_a_bool
#' @importFrom assertive assert_is_a_string
#' @importFrom assertive assert_is_all_of
#' @importFrom assertive assert_is_any_of
#' @importFrom assertive assert_is_character
#' @importFrom assertive assert_is_identical_to_na
#' @importFrom assertive assert_is_integer
#' @importFrom assertive assert_is_list
#' @importFrom assertive assert_is_non_empty
#' @importFrom assertive assert_is_subset
#' @importFrom assertive assert_is_tbl
#' @importFrom basejump assertIsAStringOrNULL
#' @importFrom rlang !!! !! .data abort inform sym syms warn
#' @importFrom utils globalVariables
"_PACKAGE"



globalVariables(".")
metadataPriorityCols <- c("sampleID", "sampleName", "description")



#' Project Directory Grep Pattern
#' @keywords internal
#' @export
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
