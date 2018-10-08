#' Transcript-to-Gene Annotations
#'
#' Generates a `Tx2Gene` object containing `transcriptID` and `geneID` columns.
#' Rownames are set to the transcript IDs.
#'
#' @note Doesn't attempt to strip transcript versions.
#'
#' @family Import/Export
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `Tx2Gene`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "tx2gene.csv")
#' x <- readTx2Gene(file)
#' print(x)
readTx2Gene <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file)
    data <- read_csv(file, col_names = c("transcriptID", "geneID"))
    # Arrange by transcript.
    data <- arrange(data, !!sym("transcriptID"))
    # Coerce to DataFrame.
    data <- as(data, "DataFrame")
    rownames(data) <- data[["transcriptID"]]
    new(Class = "Tx2Gene", data)
}
