#' Sample Directories
#'
#' @note Function will stop if no sample directories match.
#'
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return Named `character`. Sample directory paths.
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' x <- sampleDirs(uploadDir)
#' basename(x)
sampleDirs <- function(uploadDir) {
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)

    # Get the subdirectories in the upload directory.
    dirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)
    # Ensure the file paths are normalized (for Windows).
    dirs <- normalizePath(dirs, winslash = "/", mustWork = TRUE)

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
        str_trunc(toString(names(dirs)), width = 80L),
        sep = "\n"
    ))

    dirs
}
