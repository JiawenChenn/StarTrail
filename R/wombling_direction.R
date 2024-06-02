#' For checking the direction of wombling
#'
#' @param boundary boundary coordinate n*2 matrix
#' @returns a data.table containing the direction of the boundary

wombling_direction = function(boundary){
  boundary = as.matrix(boundary)
  start_point = boundary[-nrow(boundary),]
  end_point = boundary[-1,]
  mid_point = (start_point+end_point)/2
  u = end_point - start_point
  u_len = sqrt(u[,1]^2+u[,2]^2)
  u = u/u_len
  v = cbind(-u[,2],u[,1],start_point) %>% as.data.table()
  colnames(v) = c('v1','v2','start_x','start_y')
  return(v)
}