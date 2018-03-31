#' Sample Directories
#'
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @param uploadDir bcbio run upload directory.
#'
#' @return Named character vector containing sample directory paths. Function
#'   will abort if no sample directories match.
#' @export
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

    inform(paste(
        paste(length(sampleDirs), "samples detected:"),
        str_trunc(toString(names(sampleDirs)), width = 80L),
        sep = "\n"
    ))

    sampleDirs
}
