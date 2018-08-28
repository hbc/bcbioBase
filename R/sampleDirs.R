#' Sample Directories
#'
#' @note Function will stop if no sample directories match.
#'
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return Named `character` containing sample directory paths.
#' @export
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' sampleDirs(uploadDir)
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
    assert_is_non_empty(dirs)

    # Require that the sample directories are valid names in R.
    # They cannot contain non-alphanumeric characters, spaces, dashes, or begin
    # with a number. Prefix samples that start with a number using "X".
    basenames <- basename(dirs)
    # Note that multiplexed single-cell data contains a dash in the name
    # (e.g. multiplexed-AAAAAAAA).
    # We're allowing this in our checks here.
    basenames <- gsub("[-ACGT]+$", "", basenames)
    assertAllAreValidNames(basenames)
    # Our `makeNames()` function coerces periods to underscores.
    basenames <- makeNames(basenames, unique = TRUE)

    # Assign our valid names to the absolute file paths.
    names(dirs) <- basenames

    message(paste(
        paste(length(dirs), "samples detected:"),
        str_trunc(toString(names(dirs)), width = 80L),
        sep = "\n"
    ))

    dirs
}
