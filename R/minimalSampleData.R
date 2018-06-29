#' Minimal Sample Data
#'
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @param sample Sample name (e.g. "description" in bcbio YAML).
#'
#' @return `data.frame`.
#' @export
#'
#' @examples
#' minimalSampleData(c("sample 1", "sample 2"))
minimalSampleData <- function(sample) {
    assert_is_character(sample)
    assert_has_no_duplicates(sample)
    data.frame(
        sampleName = sample,
        description = sample,
        row.names = makeNames(sample),
        stringsAsFactors = TRUE
    )
}
