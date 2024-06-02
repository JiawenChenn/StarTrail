#' Calculate distance between two arrays
#'
#' @param x1 An array, each row is a data point.
#' @param x2 Another array, each row is a data point.
#' @returns A distance matrix.
#' @examples
#' cdist_r(matrix(c(1,1)), matrix(c(1,2)))
#' @import stats

cdist_r <- function(x1, x2) {
  as.matrix(stats::dist(rbind(x1, x2)))[1:nrow(x1), (nrow(x1)+1):(nrow(x1)+nrow(x2))]
}
