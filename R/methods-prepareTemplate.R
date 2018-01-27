#' Prepare R Markdown Template File
#'
#' If the required template dependency files aren't present, download latest
#' versions from the package website.
#'
#' By default, this function will create local copies of these files:
#' 
#' - `_footer.Rmd`
#' - `_header.Rmd`
#' - `_output.yaml`
#' - `bibliography.bib`
#' - `setup.R`
#'
#' @rdname prepareTemplate
#' @name prepareTemplate
#' @family Infrastructure Utilities
#'
#' @inheritParams AllGenerics
#'
#' @param object *Optional*. File name. If `NULL` (default), download the
#'   default dependency files for a new experiment.
#' @param sourceDir Source directory, typically a URL, where the dependency
#'   files are located.
#'
#' @return No value.
#'
#' @examples
#' # Copy all of the default shared template files
#' prepareTemplate()
#' unlink(c(
#'     "_footer.Rmd",
#'     "_header.Rmd",
#'     "_output.yaml",
#'     "bibliography.bib",
#'     "setup.R"
#' ))
#'
#' # Request individual files
#' prepareTemplate("setup.R")
#' unlink("setup.R")
#'
#' # Load the shared files from bcbioSingleCell
#' \dontrun{
#' prepareTemplate(
#'     sourceDir = system.file("rmarkdown/shared", 
#'                             package = "bcbioSingleCell")
#' )
#' }
#' 

NULL



# Constructors =================================================================
.copyTemplateFile <- function(object, sourceDir = NULL) {
    if (is.null(sourceDir)) {
        sourceDir <- system.file("rmarkdown/shared", package = "bcbioBase")
    }
    # Check that all source files exist
    if (!all(file.exists(file.path(sourceDir, object)))) {
        abort("Not all source template files exist")
    }
    # Note that we're not allowing accidental overwrite of locally modified
    # shared template files
    invisible(lapply(object, function(file) {
        if (!file.exists(file)) {
            file.copy(
                from = file.path(sourceDir, file),
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
            c("_output.yaml",
              "_footer.Rmd",
              "_header.Rmd",
              "bibliography.bib",
              "setup.R"),
            sourceDir = sourceDir)
    })



#' @rdname prepareTemplate
#' @export
setMethod(
    "prepareTemplate",
    signature("character"),
    .copyTemplateFile)
