#' Commands Log Parsing Functions
#'
#' @name commands
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param log `character`. Commands log.
#'
#' @return `atomic`.
#'
#' @examples
#' file <- file.path(bcbioBaseCacheURL, "surecell_commands.log")
#' log <- basejump::import(file)
#'
#' getBarcodeCutoffFromCommands(log)
#' getLevelFromCommands(log)
#' getUMITypeFromCommands(log)
NULL



#' @describeIn commands `scalar integer`.
#' @export
getBarcodeCutoffFromCommands <- function(log) {
    assert_is_character(log)
    pattern <- "--cb_cutoff (\\d+)"
    if (!any(grepl(pattern, log))) {
        stop("Failed to detect cellular barcode cutoff.")  # nocov
    }
    match <- str_match(
        string = log,
        pattern = pattern
    )
    cutoff <- match %>%
        .[, 2L] %>%
        na.omit() %>%
        unique() %>%
        as.integer()
    assert_is_an_integer(cutoff)
    message(paste(cutoff, "reads per cellular barcode cutoff detected."))
    cutoff
}



#' @describeIn commands `string`. Return `"genes"` or `"transcripts"`.
#' @export
getLevelFromCommands <- function(log) {
    assert_is_character(log)
    pattern <- "--genemap (.+)-tx2gene.tsv"
    if (any(grepl(pattern, log))) {
        level <- "genes"
    } else {
        level <- "transcripts"
    }
    message(paste0("Counts will imported as ", level, "."))
    level
}



#' @describeIn commands `string`.
#' @export
getUMITypeFromCommands <- function(log) {
    assert_is_character(log)
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
    assert_is_a_string(type)
    assert_is_non_empty(type)
    message(paste("UMI type:", type))
    type
}
