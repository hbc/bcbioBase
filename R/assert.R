#' Assert Checks
#'
#' @rdname assert
#' @name assert
#' @keywords internal
#'
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
#' @importFrom assertive assert_is_any_of
#' @importFrom assertive assert_is_character
#' @importFrom assertive assert_is_identical_to_na
#' @importFrom assertive assert_is_integer
#' @importFrom assertive assert_is_list
#' @importFrom assertive assert_is_non_empty
#' @importFrom assertive assert_is_subset
#' @importFrom assertive assert_is_tbl
#'
#' @param x Object.
#' @param severity How severe should the consequences of the assertion be?
#'   Either "`stop`", "`warning`", "`message`", or "`none`".
NULL



#' Interesting Groups Formal Assert Check
#'
#' Prevent unwanted downstream behavior when a missing interesting group
#' is requested by the user.
#'
#' @family Assert Checks
#' @inherit assert
#' @keywords internal
#'
#' @param x Object supporting [colnames()], typically a [data.frame].
#' @param interestingGroups Interesting groups character vector.
#'
#' @export
#'
#' @examples
#' demultiplexed <- file.path(
#'     "http://bcbiobase.seq.cloud",
#'     "sample_metadata",
#'     "demultiplexed.xlsx")
#' meta <- readSampleMetadataFile(demultiplexed)
#' assert_formal_interesting_groups(meta, "genotype")
assert_formal_interesting_groups <- function(
    x,
    interestingGroups,
    severity = "stop") {
    assert_has_colnames(x, severity = severity)
    assert_is_subset(
        x = interestingGroups,
        y = colnames(x),
        severity = severity
    )
}
