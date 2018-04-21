#' Prepare R Markdown Template File
#'
#' If the required template dependency files aren't present, download latest
#' versions from the package website. Existing files are never overwritten.
#'
#' By default, this function will create local copies of these files:
#'
#' - `_footer.Rmd`
#' - `_header.Rmd`
#' - `_output.yaml`
#' - `_setup.R`
#' - `bibliography.bib`
#'
#' @family Prepare Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param file File name.
#' @param sourceDir Source directory.
#'
#' @return No value.
#' @export
#'
#' @examples
#' defaultFiles <- c(
#'     "_footer.Rmd",
#'     "_header.Rmd",
#'     "_output.yaml",
#'     "_setup.R",
#'     "bibliography.bib"
#' )
#'
#' # Copy all of the default shared template files
#' prepareTemplate()
#' file.exists(defaultFiles)
#' unlink(defaultFiles)
#'
#' # Request individual files
#' prepareTemplate("bibliography.bib")
#' file.exists("bibliography.bib")
#' unlink("bibliography.bib")
prepareTemplate <- function(
    file = c(
        "_output.yaml",
        "_footer.Rmd",
        "_header.Rmd",
        "_setup.R",
        "bibliography.bib"
    ),
    sourceDir = system.file("rmarkdown/shared", package = "bcbioBase")
) {
    assert_is_character(file)
    assert_all_are_dirs(sourceDir)
    assert_is_a_string(sourceDir)
    assert_all_are_existing_files(file.path(sourceDir, file))
    invisible(mapply(
        FUN = function(file, sourceDir) {
            if (!file.exists(file)) {
                file.copy(
                    from = file.path(sourceDir, file),
                    to = file,
                    overwrite = FALSE
                )
            }
        },
        file = file,
        MoreArgs = list(sourceDir = sourceDir)
    ))
}
