#' Project Directory
#'
#' @note This will pick the latest dated directory and warn the user, if bcbio
#' has been run multiple times to the same upload directory.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams params
#'
#' @return `string`. Dated project directory (e.g. "2018-01-01_rnaseq").
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' x <- projectDir(uploadDir)
#' basename(x)
projectDir <- function(uploadDir) {
    dir <- list.files(
        path = uploadDir,
        pattern = projectDirPattern,
        full.names = FALSE,
        recursive = FALSE
    )
    assert_is_non_empty(dir)
    # Check to see if user has run bcbio multiple times to the same upload
    # directory, and warn when this is detected.
    if (length(dir) > 1L) {
        newest <- tail(dir, n = 1L)
        warning(paste(
            "Multiple project directories detected:",
            printString(dir),
            paste("Using most recent:", newest),
            sep = "\n"
        ), call. = FALSE)
        dir <- newest
    }
    assert_is_a_string(dir)
    message(paste("Dated project directory:", dir))
    realpath(file.path(uploadDir, dir))
}



# Alternatively, can parse YAML data, but this is faster and simpler.
#' Run Date
#'
#' Get the run date from the project directory.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams params
#'
#' @examples
#' runDate("2018-01-01_illumina_rnaseq")
runDate <- function(projectDir) {
    assert_is_a_string(projectDir)
    projectDir <- basename(projectDir)
    assert_that(grepl(projectDirPattern, projectDir))
    match <- str_match(
        string = projectDir,
        pattern = projectDirPattern
    )
    as.Date(match[[2L]])
}



#' Sample Directories
#'
#' @note Function will stop if no sample directories match.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams params
#'
#' @return Named `character`. Sample directory paths.
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' x <- sampleDirs(uploadDir)
#' basename(x)
sampleDirs <- function(uploadDir) {
    assert_all_are_dirs(uploadDir)
    uploadDir <- realpath(uploadDir)

    # Get the subdirectories in the upload directory.
    dirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)

    # Detect and remove nested dated project directory.
    projectDir <- suppressMessages(projectDir(uploadDir))
    dirs <- setdiff(dirs, projectDir)

    # Double check that we're nuking any remaining dated directories, in case
    # bcbio has been run multiple times.
    isSample <- !grepl(
        pattern = projectDirPattern,
        x = basename(dirs)
    )
    dirs <- dirs[isSample]

    # Ensure there are sample directories in the upload.
    assert_is_non_empty(dirs)

    # Use the directory basenames for vector names.
    basenames <- basename(dirs)

    # Require that the sample directories are valid names in R. They cannot
    # contain non-alphanumeric characters, spaces, dashes, or begin with a
    # number. Prefix samples that start with a number using "X". Note that
    # multiplexed single-cell samples are expected to contain a dash in the name
    # (e.g. multiplexed-AAAAAAAA).
    x <- basenames
    x <- gsub("[-ACGT]+$", "", basenames)
    assertAllAreValidNames(x)

    # Our `makeNames()` function coerces periods and dashes to underscores.
    basenames <- makeNames(basenames, unique = TRUE)

    # Assign our valid names to the absolute file paths.
    names(dirs) <- basenames

    message(paste(
        paste(length(dirs), "sample(s) detected:"),
        str_trunc(toString(names(dirs)), width = getOption("width")),
        sep = "\n"
    ))

    dirs
}
