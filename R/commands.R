#' Commands log parsing functions
#'
#' @name commands
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param log `character`.
#'   Commands log.
#'
#' @return `atomic`.
#'
#' @examples
#' file <- file.path(bcbioBaseTestsURL, "surecell-commands.log")
#' log <- basejump::import(file)
#' getBarcodeCutoffFromCommands(log)
#' getLevelFromCommands(log)
#' getUMITypeFromCommands(log)
NULL



#' @describeIn commands `integer(1)`.
#' @export
getBarcodeCutoffFromCommands <- function(log) {
    assert(isCharacter(log))
    pattern <- "--cb_cutoff (\\d+)"
    if (!any(grepl(pattern, log))) {
        stop("Failed to detect cellular barcode cutoff.")  # nocov
    }
    cutoff <- log %>%
        str_match(pattern) %>%
        .[, 2L] %>%
        na.omit() %>%
        unique() %>%
        as.integer()
    assert(isInt(cutoff))
    message(paste(cutoff, "reads per cellular barcode cutoff detected."))
    cutoff
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
    message(paste0("Counts will imported as ", level, "."))
    level
}



#' @describeIn commands `character(1)`.
#' @export
getUMITypeFromCommands <- function(log) {
    assert(isCharacter(log))
    pattern <- "fastqtransform.*/(.*)\\.json"
    if (!any(grepl(pattern, log))) {
        stop("Failed to detect UMI type.")
    }
    type <- log %>%
        str_match(pattern = pattern) %>%
        .[, 2L] %>%
        na.omit() %>%
        unique() %>%
        str_replace(pattern = "-transform", replacement = "")
    assert(isString(type))
    message(paste("UMI type:", type))
    type
}
