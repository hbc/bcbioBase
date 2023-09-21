#' Commands log parsing functions
#'
#' @name commands
#' @note Updated 2023-09-21.
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param log `character`.
#' Commands log.
#'
#' @return `atomic`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "surecell-commands.log")
#' log <- import(file)
#' getBarcodeCutoffFromCommands(log)
#' getLevelFromCommands(log)
#' getUMITypeFromCommands(log)
NULL



#' @describeIn commands `integer(1)`.
#' @export
getBarcodeCutoffFromCommands <- function(log) {
    assert(isCharacter(log))
    pattern <- "--cb_cutoff (\\d+)"
    assert(
        any(grepl(pattern, log)),
        msg = "Failed to detect cellular barcode cutoff."
    )
    x <- strMatch(x = log, pattern = pattern)
    x <- x[, 2L]
    x <- as.integer(unique(na.omit(x)))
    assert(isInt(x))
    alertInfo(sprintf(
        "%d %s per cellular barcode cutoff detected.",
        x, ngettext(n = x, msg1 = "read", msg2 = "reads")
    ))
    x
}



#' @describeIn commands `character(1)`. Return `"genes"` or `"transcripts"`.
#' @export
getLevelFromCommands <- function(log) {
    assert(isCharacter(log))
    pattern <- "--genemap (.+)-tx2gene.tsv"
    if (any(grepl(pattern, log))) {
        level <- "genes"
    } else {
        level <- "transcripts"
    }
    message(sprintf("Counts will imported as %s.", level))
    level
}



#' @describeIn commands `character(1)`.
#' @export
getUMITypeFromCommands <- function(log) {
    assert(isCharacter(log))
    pattern <- "fastqtransform.*/(.*)\\.json"
    assert(
        any(grepl(pattern, log)),
        msg = "Failed to detect UMI type."
    )
    x <- strMatch(x = log, pattern = pattern)
    x <- x[, 2L]
    x <- unique(na.omit(x))
    x <- sub(pattern = "-transform", replacement = "", x = x)
    assert(isString(x))
    dl(c(`UMI type` = x))
    x
}
