# FIXME Add working examples and code coverage.



#' Get from Commands Log
#'
#' @name commands
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param log `character`. Commands log.
#'
#' @return `integer`.
NULL



#' @rdname commands
#' @export
getBarcodeCutoffFromCommands <- function(log) {
    assert_is_character(log)
    pattern <- "--cb_cutoff (\\d+)"
    assert_any_are_matching_regex(x = log, pattern = pattern)
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



#' Get Gene or Transcript Level from Commands Log
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @export
#'
#' @param log `character`. Commands log.
#'
#' @return `string`.
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



#' Get UMI Type from Commands Log
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @export
#'
#' @param log `character`. Commands log.
#'
#' @return `string`.
getUMITypeFromCommands <- function(log) {
    assert_is_character(log)
    pattern <- "fastqtransform.*/(.*)\\.json"
    assert_any_are_matching_regex(x = log, pattern = pattern)
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
