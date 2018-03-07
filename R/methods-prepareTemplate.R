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
#' @name prepareTemplate
#' @family R Markdown Functions
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
#' file.exists(defaultFiles)
#' unlink(defaultFiles)
#'
#' # Request individual files
#' prepareTemplate("bibliography.bib")
#' file.exists("bibliography.bib")
#' unlink("bibliography.bib")
#'
#' # Load the shared files from bcbioSingleCell
#' \dontrun{
#' prepareTemplate(
#'     sourceDir = system.file(
#'         "rmarkdown/shared",
#'         package = "bcbioSingleCell"
#'     )
#' )
#' }
NULL



# Constructors =================================================================
.copyTemplateFile <- function(object, sourceDir = NULL) {
    assert_is_character(object)
    assertIsAStringOrNULL(sourceDir)
    if (is.null(sourceDir)) {
        sourceDir <- system.file("rmarkdown/shared", package = "bcbioBase")
    }
    assert_all_are_existing_files(file.path(sourceDir, object))
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
        file = object,
        MoreArgs = list(sourceDir = sourceDir)
    ))
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
    }
)



#' @rdname prepareTemplate
#' @export
setMethod(
    "prepareTemplate",
    signature("character"),
    .copyTemplateFile
)
