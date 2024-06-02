#' Wombling method to calculate the gradient of the boundary
#'
#' @param model a fitted NNGP model
#' @param boundary boundary(a curve), each row is a point on the curve, each segment is defined by two neighboring points
#' @param h step size used in finite difference, default is 0.8* minimum distance between data points
#' @param threads number of threads used in NNGP model
#' @param path path to save the results
#' @param prefix prefix of the file name
#' @param save_file whether to save the results
#' @param verbose whether to print out the progress
#' @returns a data.table containing the mean prediction and gradient
#'
#' @import data.table
#' @import spNNGP
#' @import stats
#'

wombling = function(model,boundary,h,threads,path="./",prefix = 'boundary0',save_file=FALSE,verbose=FALSE){
    boundary <- as.matrix(boundary)

    start_point = boundary[-nrow(boundary),]
    end_point = boundary[-1,]
    mid_point = (start_point+end_point)/2
    u = end_point - start_point
    u_len = sqrt(u[,1]^2+u[,2]^2)
    u = u/u_len
    v = cbind(-u[,2],u[,1])

    f_a = predict(model, X.0 = matrix(1.0,nr=nrow(start_point)), coords.0 = start_point, n.omp.threads=threads,verbose=verbose)
    f_b = predict(model, X.0 = matrix(1.0,nr=nrow(end_point)), coords.0 = end_point, n.omp.threads=threads,verbose=verbose)
    f_mid = predict(model, X.0 = matrix(1.0,nr=nrow(mid_point)), coords.0 = mid_point, n.omp.threads=threads,verbose=verbose)
    f_ahv = predict(model, X.0 = matrix(1.0,nr=nrow(start_point)), coords.0 = start_point+h*v, n.omp.threads=threads,verbose=verbose)
    f_bhv = predict(model, X.0 = matrix(1.0,nr=nrow(end_point)), coords.0 = end_point+h*v, n.omp.threads=threads,verbose=verbose)
    f_midhv = predict(model, X.0 = matrix(1.0,nr=nrow(mid_point)), coords.0 = mid_point+h*v, n.omp.threads=threads,verbose=verbose)

    f_a = f_a$p.y.0
    f_b = f_b$p.y.0
    f_mid = f_mid$p.y.0
    f_ahv = f_ahv$p.y.0
    f_bhv = f_bhv$p.y.0
    f_midhv = f_midhv$p.y.0

    g_a = (f_ahv-f_a)/h
    g_b = (f_bhv-f_b)/h
    g_mid = (f_midhv-f_mid)/h
    g_simpson = u_len * (g_a+4*g_mid+g_b)/6
    if(save_file){data.table::fwrite(data.table(g_simpson), paste0(path,"/",prefix,'_wombling_all.txt'), row.names=F, col.names=F, sep="\t", quote=F)}

    g_a = apply(g_a, 1, mean)
    g_b = apply(g_b, 1, mean)
    g_mid = apply(g_mid, 1, mean)
    g_simpson = apply(g_simpson, 1, mean)


    boundary_result <- cbind(g_simpson, u_len,start_point,end_point)
    colnames(boundary_result) = c('g_simpson','u_len','start_x','start_y','end_x','end_y')
    if(save_file){data.table::fwrite(boundary_result, paste0(path,"/",prefix,'_wombling.txt'), row.names=F, col.names=T, sep="\t", quote=F)}
    return(boundary_result)
}
