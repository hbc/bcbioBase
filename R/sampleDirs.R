#' Sample directories
#'
#' @author Michael Steinbaugh
#' @note Updated 2019-08-05.
#' @export
#'
#' @inheritParams acidroxygen::params
#'
#' @note Function will [`stop()`][base::stop] if no sample directories match.
#'
#' @return Named `character`. Sample directory paths.
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' x <- sampleDirs(uploadDir)
#' basename(x)
sampleDirs <- function(uploadDir) {
    assert(isADirectory(uploadDir))
    uploadDir <- realpath(uploadDir)

    ## Get the subdirectories in the upload directory.
    dirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)

    ## Detect and remove nested dated project directory.
    projectDir <- suppressMessages(projectDir(uploadDir))
    dirs <- setdiff(dirs, projectDir)

    ## Double check that we're nuking any remaining dated directories, in case
    ## bcbio has been run multiple times.
    isSample <- !grepl(
        pattern = projectDirPattern,
        x = basename(dirs)
    )
    dirs <- dirs[isSample]

    ## Ensure there are sample directories in the upload.
    assert(isNonEmpty(dirs))

    ## Use the directory basenames for vector names.
    basenames <- basename(dirs)

    ## Return sample directory basenames that are valid names in R.
    ## We're using these to define the `sampleID` column in our object metadata.
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
        message(sprintf(
            "Sanitizing sample names: %s.",
            toString(invalid, width = 100L)
        ))
    }

    ## Our `makeNames` function coerces periods and dashes to underscores.
    basenames <- makeNames(basenames)

    ## Assign our valid names to the absolute file paths.
    names(dirs) <- basenames

    message(sprintf(
        fmt = "%d %s detected:\n%s",
        length(dirs),
        ngettext(n = length(dirs), msg1 = "sample", msg2 = "samples"),s
        printString(sort(names(dirs)))
    ))

    dirs
}
