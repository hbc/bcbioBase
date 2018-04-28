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
#' @inheritParams general
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
#' @inherit assert
#'
#' @inheritParams general
#'
#' @return Silent, stop on error.
#' @export
#'
#' @examples
#' assertFormalInterestingGroups(rse_bcb, "treatment")
#' assertFormalInterestingGroups(rse_dds, "condition")
assertFormalInterestingGroups <- function(
    x,
    interestingGroups,
    severity = getOption("assertive.severity", "stop")
) {
    fun <- get(severity)

    # Early return on `NULL` value (e.g. DESeqDataSet)
    if (is.null(interestingGroups)) {
        return(invisible())
    }

    assert_is_character(interestingGroups)

    # Obtain sampleData is S4 object is passed in
    if (!any(is(x, "DataFrame") || is(x, "data.frame"))) {
        # Don't want clean return, so we can check to see if interesting
        # groups are present but defined as non factor columns
        x <- sampleData(x, clean = FALSE, interestingGroups = NULL)
    }

    # Check that interesting groups are slotted into sampleData
    if (!all(interestingGroups %in% colnames(x))) {
        setdiff <- setdiff(interestingGroups, colnames(x))
        fun(paste(
            "The interesting groups",
            deparse(toString(setdiff)),
            "are not defined as columns in `sampleData()`"
        ))
    }

    # Check that interesting groups are factors
    isFactor <- vapply(
        X = x[, interestingGroups, drop = FALSE],
        FUN = is.factor,
        FUN.VALUE = logical(1L),
        USE.NAMES = TRUE
    )
    if (!all(isFactor)) {
        invalid <- names(isFactor)[which(!isFactor)]
        fun(paste(
            "The interesting groups",
            deparse(toString(invalid)),
            "are not factor"
        ))
    }

    # Don't allow interesting groups from our defined blacklist
    if (any(interestingGroups %in% metadataBlacklist)) {
        intersect <- intersect(interestingGroups, metadataBlacklist)
        fun(paste(
            paste(
                "The interesting groups",
                deparse(toString(intersect)),
                "are blacklisted."
            ),
            "Blacklist (bcbioBase::metadataBlacklist) :",
            basejump::printString(metadataBlacklist),
            sep = "\n"
        ))
    }
}
