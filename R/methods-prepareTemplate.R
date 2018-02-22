#' Prepare R Markdown Template File
#'
#' If the required template dependency files aren't present, download latest
#' versions from the package website. *Existing files are never overwritten.*
#'
#' By default, this function will create local copies of these files:
#'
#' - `_footer.Rmd`
#' - `_header.Rmd`
#' - `_output.yaml`
#' - `_setup.R`
#' - `bibliography.bib`
#'
#' @rdname prepareTemplate
#' @name prepareTemplate
#' @family Infrastructure Utilities
#'
#' @importFrom fs file_copy file_exists path
#'
#' @inheritParams general
#'
#' @param object *Optional*. File name. If `NULL` (default), download the
#'   default dependency files for a new experiment.
#' @param sourceDir Source directory, typically a URL, where the dependency
#'   files are located.
#'
#' @return No value.
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
#' file_exists(defaultFiles)
#' file_delete(defaultFiles)
#'
#' # Request individual files
#' prepareTemplate("bibliography.bib")
#' file_exists("bibliography.bib")
#' file_delete("bibliography.bib")
#'
#' # Load the shared files from bcbioSingleCell
#' \dontrun{
#' prepareTemplate(
#'     sourceDir = system.file(
#'         "rmarkdown/shared",
#'         package = "bcbioSingleCell")
#' )
#' }
NULL



# Constructors =================================================================
.copyTemplateFile <- function(object, sourceDir = NULL) {
    assert_is_character(object)
    if (is.null(sourceDir)) {
        sourceDir <- system.file("rmarkdown/shared", package = "bcbioBase")
    }
    assert_is_a_string(sourceDir)
    assert_all_are_existing_files(path(sourceDir, object))
    invisible(lapply(object, function(file) {
        if (!file_exists(file)) {
            file_copy(
                from = path(sourceDir, file),
                to = file,
                overwrite = FALSE)
        }
    }))
}



# Methods ======================================================================
#' @rdname prepareTemplate
#' @export
setMethod(
    "prepareTemplate",
    signature("missing"),
    function(
        object,
        sourceDir = NULL) {
        .copyTemplateFile(
            c(
                "_output.yaml",
                "_footer.Rmd",
                "_header.Rmd",
                "_setup.R",
                "bibliography.bib"
            ),
            sourceDir = sourceDir)
    })



#' @rdname prepareTemplate
#' @export
setMethod(
    "prepareTemplate",
    signature("character"),
    .copyTemplateFile)
