# v0.2.0
# nocov start

#' @rdname deprecated
#' @export
setGeneric("flatFiles", function(object, ...) {
    standardGeneric("flatFiles")
})

#' @rdname deprecated
#' @export
setMethod(
    "flatFiles",
    signature("SummarizedExperiment"),
    function(object) {
        .Deprecated("as(object, \"list\")")
        as(object, "list")
    }
)

# nocov end
