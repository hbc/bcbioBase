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
#'     x = colData(rse_bcb),
#'     interestingGroups = "treatment"
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
