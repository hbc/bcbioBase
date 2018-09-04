#' Minimal Sample Data
#'
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @param sample `character`. Sample names (e.g. "description" in bcbio YAML).
#'
#' @return `DataFrame`.
#' @export
#'
#' @examples
#' minimalSampleData(c("sample 1", "sample 2"))
minimalSampleData <- function(sample) {
    assert_is_character(sample)
    assert_has_no_duplicates(sample)
    DataFrame(
        sampleName = sample,
        description = sample,
        row.names = makeNames(sample)
    )
}
