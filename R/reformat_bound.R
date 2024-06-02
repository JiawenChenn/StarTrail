#' Format the boundary matrix to a matrix of segments
#'
#' @param bound coordinate matrix of the boundary, each row is a data point (x,y), # * 2 matrix
#' @returns a matrix of the boundary, each row is a segment of the boundary, # * 4 matrix

reformat_bound = function(bound){
  bound = bound %>% as.data.frame(check.names=F)
  bound_new = cbind(bound[-nrow(bound),],bound[-1,])
  colnames(bound_new) = c('x_start','y_start','x_end','y_end')
  return(bound_new)
}
