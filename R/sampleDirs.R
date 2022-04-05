#' Sample directories
#'
#' @author Michael Steinbaugh
#' @note Updated 2022-03-07.
#' @export
#'
#' @inheritParams AcidRoxygen::params
#'
#' @note Function will intentionally error when no sample directories match.
#'
#' @return Named `character`.
#' Sample directory paths.
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' x <- sampleDirs(uploadDir)
#' basename(x)
sampleDirs <- function(uploadDir) {
    assert(isADirectory(uploadDir))
    uploadDir <- realpath(uploadDir)
    ## Get the subdirectories in the upload directory.
    dirs <- sort(list.dirs(
        path = uploadDir,
        full.names = TRUE,
        recursive = FALSE
    ))
    ## Detect and remove nested dated project directory.
    suppressMessages({
        projectDir <- projectDir(uploadDir)
    })
    dirs <- setdiff(x = dirs, y = projectDir)
    assert(hasLength(dirs))
    ## Double check that we're nuking any remaining dated directories, in case
    ## bcbio has been run multiple times.
    keep <- !grepl(pattern = projectDirPattern, x = basename(dirs))
    dirs <- dirs[keep]
    assert(hasLength(dirs))
    ## Exclude any nested directories from bcbio-nextgen pipelines, such as new
    ## `bcbioRNASeq/` subdirectory, which was added in 2021.
    denylist <- c(
        "bcbioRNASeq",
        "bcbioSingleCell"
    )
    keep <- !basename(dirs) %in% denylist
    dirs <- dirs[keep]
    assert(hasLength(dirs))
    ## Use the directory basenames for vector names.
    basenames <- basename(dirs)
    ## Return sample directory basenames that are valid names in R.
    ## We're using these to define the `sampleId` column in our object metadata.
    ## Refer to `make.names()` for the valid name conventions in R.
    ## We're using the `makeNames()` variant here instead, which sanitizes using
    ## an underscore instead of a period.
    ##
    ## In particular, these are potentially problematic:
    ## - Contains dashes/hyphens (very common).
    ## - Begins with a number (very common).
    ## - Contains non-alphanumerics.
    ##
    ## Here we are informing the user when bcbio samples need to get sanitized
    ## to be valid in R.
    check <- basenames
    ## Note that multiplexed single-cell samples are expected to contain a dash
    ## in the name (e.g. multiplexed-AAAAAAAA), so we're removing this from the
    ## validity check.
    check <- gsub("[-ACGT]+$", "", basenames)
    if (!isTRUE(validNames(check))) {
        invalid <- setdiff(check, makeNames(check))
        alertInfo(sprintf(
            "Sanitizing sample names: %s.",
            toInlineString(invalid, n = 5L, class = "val")
        ))
    }
    ## Our `makeNames` function coerces periods and dashes to underscores.
    basenames <- makeNames(basenames)
    ## Assign our valid names to the absolute file paths.
    names(dirs) <- basenames
    alertInfo(sprintf(
        fmt = "%d %s detected:",
        length(dirs),
        ngettext(n = length(dirs), msg1 = "sample", msg2 = "samples")
    ))
    ul(sort(names(dirs)))
    dirs
}
