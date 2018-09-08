#' Transcript to Gene Annotations
#'
#' Generates a `tx2gene` object containing `transcriptID` and `geneID` columns.
#' Rownames are set to the transcript IDs.
#'
#' @note Doesn't attempt to strip transcript versions.
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `tx2gene`.
#' @export
#'
#' @examples
#' x <- readTx2gene("http://bcbiobase.seq.cloud/tx2gene.csv")
#' print(x)
readTx2gene <- function(file) {
    assert_is_a_string(file)
    file <- localOrRemoteFile(file)
    data <- read_csv(file, col_names = c("transcriptID", "geneID"))
    # Arrange by transcript.
    data <- arrange(data, !!sym("transcriptID"))
    # Coerce to DataFrame.
    data <- as(data, "DataFrame")
    rownames(data) <- data[["transcriptID"]]
    new("tx2gene", data)
}
