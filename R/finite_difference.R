#' Calculate gradient using finite difference
#'
#' @param coord coordinates of the data points, each row is a data point, sample * 2 matrix
#' @param h step size used in finite difference, default is 0.8* minimum distance between data points
#' @param model a fitted NNGP model
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

finite_difference = function(coord,h,model,threads=8,path='./',prefix='ori',save_file=FALSE,verbose=FALSE){
    f =    predict(model, X.0 = matrix(1.0,nr=nrow(coord)), coords.0 = coord, n.omp.threads=threads,verbose=verbose)
    f_h0 = predict(model, X.0 = matrix(1.0,nr=nrow(coord)), 
                   coords.0 = coord+h*cbind(rep(1,nrow(coord)),rep(0,nrow(coord))), n.omp.threads=threads,verbose=verbose)
    f_0h = predict(model, X.0 = matrix(1.0,nr=nrow(coord)),
                     coords.0 = coord+h*cbind(rep(0,nrow(coord)),rep(1,nrow(coord))), n.omp.threads=threads,verbose=verbose)

    f = f$p.y.0
    f_h0 = f_h0$p.y.0
    f_0h = f_0h$p.y.0

    g1 = (f_h0-f)/h
    g2 = (f_0h-f)/h
    
    if(save_file){
    data.table::fwrite(f,paste0(path,"/",prefix,'_mean.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    data.table::fwrite(g1,paste0(path,"/",prefix,'_g1.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    data.table::fwrite(g2,paste0(path,"/",prefix,'_g2.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    }

    g1_min = apply(g1, 1, min)
    g2_min = apply(g2, 1, min)
    g1_max = apply(g1, 1, max)
    g2_max = apply(g2, 1, max)

    f = apply(f, 1, mean)
    g1 = apply(g1, 1, mean)
    g2 = apply(g2, 1, mean)
    
    result = data.table::data.table(mean=f,g1=g1,g2=g2,g1_min=g1_min,g1_max=g1_max,g2_min=g2_min,g2_max=g2_max)
    colnames(result) = c('pred','g1','g2','g1_min','g1_max','g2_min','g2_max')
    return(result)
}
