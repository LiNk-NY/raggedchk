#' Utility function to measure object size
#'
#' @description `object_size` works on an object and provides the size of the
#'   object using `utils::object.size`. `as.object_size` converts a numeric
#'   value to the units indicated in the `unit` argument.
#'
#' @inheritParams utils::object.size
#'
#' @return A scalar character vector of the input object's size in the units
#'   requested
#'
#' @examples
#' var <- 1e10
#' object_size(var, unit = "KB")
#'
#' as.object_size(2500000)
#' @export
object_size <- function(x, unit = "MB", standard = "SI")
{
    x <- utils::object.size(x)
    class(x) <- "object_size"
    format(x, units = unit, standard = standard)
}

#' @rdname object_size
#' @export
as.object_size <- function(x, unit = "MB", standard = "SI")
{
    if (!is.numeric(x))
        stop("'x' must be numeric")
    class(x) <- "object_size"
    format(x, units = unit, standard = standard)
}
