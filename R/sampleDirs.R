#' Sample Directories
#'
#' @note Function will stop if no sample directories match.
#'
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @param uploadDir File path to bcbio run upload directory.
#'
#' @return Named `character` vector containing sample directory paths.
#' @export
#'
#' @examples
#' uploadDir <- system.file("extdata/bcbio", package = "bcbioBase")
#' sampleDirs(uploadDir)
sampleDirs <- function(uploadDir) {
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)

    # Get the subdirectories in the upload directory
    subdirs <- list.dirs(uploadDir, full.names = TRUE, recursive = FALSE)

    # Require detection and removal of nested `projectDir`
    projectDir <- grep(
        pattern = projectDirPattern,
        x = basename(subdirs)
    )
    assert_is_non_empty(projectDir)

    sampleDirs <- subdirs[-projectDir]
    assert_is_non_empty(sampleDirs)

    # Generate names from file paths and make valid
    names(sampleDirs) <- makeNames(basename(sampleDirs), unique = TRUE)

    message(paste(
        paste(length(sampleDirs), "samples detected:"),
        str_trunc(toString(names(sampleDirs)), width = 80L),
        sep = "\n"
    ))

    sampleDirs
}
