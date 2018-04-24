#' Assert Checks
#'
#' @name assert
#' @keywords internal
#'
#' @inheritParams general
#' @param severity How severe should the consequences of the assertion be?
#'   Either "`stop`", "`warning`", "`message`", or "`none`".
#'
#' @return Silent on success, stop on failure.
NULL



#' Assert Formal Annotation Column
#'
#' @family Assert Check Functions
#' @author Michael Steinbaugh
#' @inherit assert
#'
#' @param x Object supporting `dimnames`, such as a `matrix` or `data.frame`.
#' @param sampleData Sample data.
#'
#' @export
#'
#' @examples
#' x <- data.frame(
#'     "sample_1" = c(1L, 2L),
#'     "sample_2" = c(3L, 4L),
#'     row.names = c("gene_1", "gene_2"),
#'     stringsAsFactors = FALSE
#' )
#' sampleData <- data.frame(
#'     "genotype" = c("wt", "ko"),
#'     row.names = c("sample_1", "sample_2"),
#'     stringsAsFactors = TRUE
#' )
#' assertFormalAnnotationCol(x, sampleData)
assertFormalAnnotationCol <- function(
    x,
    sampleData,
    severity = getOption("assertive.severity", "stop")
) {
    assert_has_dimnames(x, severity = severity)
    assert_is_any_of(
        x = sampleData,
        classes = c("data.frame", "DataFrame", "logical", "NULL"),
        severity = severity
    )
    if (has_dims(sampleData)) {
        assert_has_dimnames(sampleData, severity = severity)
        assert_are_identical(
            x = colnames(x),
            y = rownames(sampleData),
            severity = severity
        )
        # All columns must be factors
        lapply(sampleData, function(x) {
            assert_is_factor(x, severity = severity)
        })
    }
    if (is.logical(sampleData)) {
        assert_is_identical_to_na(sampleData, severity = severity)
    }
}



#' Interesting Groups Formal Assert Check
#'
#' Prevent unwanted downstream behavior when a missing interesting group
#' is requested by the user.
#'
#' @family Assert Check Functions
#' @author Michael Steinbaugh
#'
#' @inherit assert
#'
#' @inheritParams general
#' @param x Object supporting [colnames()], typically a `data.frame`.
#'
#' @return Silent, stop on error.
#' @export
#'
#' @examples
#' assertFormalInterestingGroups(
#'     x = sampleData(rse_dds),
#'     interestingGroups = "condition"
#' )
assertFormalInterestingGroups <- function(
    x,
    interestingGroups,
    severity = "stop"
) {
    assert_has_colnames(x, severity = severity)
    assert_is_subset(
        x = interestingGroups,
        y = colnames(x),
        severity = severity
    )
}
