#.libPaths('/proj/yunligrp/users/jiawen/software/Rlibs')
library(spNNGP)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggsci)
library(egg)
library(paletteer)


cdist_r <- function(x1, x2) {
  as.matrix(dist(rbind(x1, x2)))[1:nrow(x1), (nrow(x1)+1):(nrow(x1)+nrow(x2))]
}

finite_difference = function(coord,h,model,threads=8,path='./',prefix='ori'){
    f =    predict(model, X.0 = matrix(1.0,nr=nrow(coord)), coords.0 = coord, n.omp.threads=threads)
    f_h0 = predict(model, X.0 = matrix(1.0,nr=nrow(coord)), 
                   coords.0 = coord+h*cbind(rep(1,nrow(coord)),rep(0,nrow(coord))), n.omp.threads=threads)
    f_0h = predict(model, X.0 = matrix(1.0,nr=nrow(coord)),
                     coords.0 = coord+h*cbind(rep(0,nrow(coord)),rep(1,nrow(coord))), n.omp.threads=threads)

    f = f$p.y.0
    f_h0 = f_h0$p.y.0
    f_0h = f_0h$p.y.0

    g1 = (f_h0-f)/h
    g2 = (f_0h-f)/h

    fwrite(f,paste0(path,"/",prefix,'_mean.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    fwrite(g1,paste0(path,"/",prefix,'_g1.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    fwrite(g2,paste0(path,"/",prefix,'_g2.txt'),row.names=F,col.names=F,sep="\t",quote=F)

    g1_min = apply(g1, 1, min)
    g2_min = apply(g2, 1, min)
    g1_max = apply(g1, 1, max)
    g2_max = apply(g2, 1, max)

    f = apply(f, 1, mean)
    g1 = apply(g1, 1, mean)
    g2 = apply(g2, 1, mean)


    return(data.table(mean=f,g1=g1,g2=g2,g1_min=g1_min,g1_max=g1_max,g2_min=g2_min,g2_max=g2_max))
}

finite_difference_curvature = function(coord,h,model,threads=8,path='./',prefix='ori'){
    f =    predict(model, X.0 = matrix(1.0,nr=nrow(coord)), coords.0 = coord, n.omp.threads=threads)
    f_hh = predict(model, X.0 = matrix(1.0,nr=nrow(coord)), 
                   coords.0 = coord+h*cbind(rep(1,nrow(coord)),rep(1,nrow(coord))), n.omp.threads=threads)
    f_h0 = predict(model, X.0 = matrix(1.0,nr=nrow(coord)), 
                   coords.0 = coord+h*cbind(rep(1,nrow(coord)),rep(0,nrow(coord))), n.omp.threads=threads)
    f_0h = predict(model, X.0 = matrix(1.0,nr=nrow(coord)),
                     coords.0 = coord+h*cbind(rep(0,nrow(coord)),rep(1,nrow(coord))), n.omp.threads=threads)
    f_2h0= predict(model, X.0 = matrix(1.0,nr=nrow(coord)), 
                    coords.0 = coord+h*cbind(rep(2,nrow(coord)),rep(0,nrow(coord))), n.omp.threads=threads)
    f_02h= predict(model, X.0 = matrix(1.0,nr=nrow(coord)),
                    coords.0 = coord+h*cbind(rep(0,nrow(coord)),rep(2,nrow(coord))), n.omp.threads=threads)

    f = f$p.y.0
    f_hh = f_hh$p.y.0
    f_h0 = f_h0$p.y.0
    f_0h = f_0h$p.y.0
    f_2h0 = f_2h0$p.y.0
    f_02h = f_02h$p.y.0

    g1 = (f_h0-f)/h
    g2 = (f_0h-f)/h
    c11 = (f_2h0-2*f_h0+f)/(h^2)
    c22 = (f_02h-2*f_0h+f)/(h^2)
    c12 = (f_hh-f_h0-f_0h+f)/(h^2)
    #fwrite(f,paste0(path,"/",prefix,'_mean.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    #fwrite(g1,paste0(path,"/",prefix,'_g1.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    #fwrite(g2,paste0(path,"/",prefix,'_g2.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    #fwrite(c11,paste0(path,"/",prefix,'_c11.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    #fwrite(c22,paste0(path,"/",prefix,'_c22.txt'),row.names=F,col.names=F,sep="\t",quote=F)
    #fwrite(c12,paste0(path,"/",prefix,'_c12.txt'),row.names=F,col.names=F,sep="\t",quote=F)

    g1 = apply(g1, 1, mean)
    g2 = apply(g2, 1, mean)
    c11 = apply(c11, 1, mean)
    c22 = apply(c22, 1, mean)
    c12 = apply(c12, 1, mean)
    f = apply(f, 1, mean)

    return(data.table(mean=f,g1=g1,g2=g2,c11=c11,c22=c22,c12=c12))
}

coord_filter <- function(coord_mesh, coords, close_dis = 0.2, distance_to_boundary = 0) {
  
  # Ensure the inputs are matrices
  coord_mesh <- as.matrix(coord_mesh)
  coords <- as.matrix(coords)
  
  n <- nrow(coord_mesh)
  indices <- logical(n)
  
  for (i in 1:n) {
    close_point_indices_x <- which(abs(coords[, 1] - coord_mesh[i, 1]) < close_dis)
    close_point_indices_y <- which(abs(coords[, 2] - coord_mesh[i, 2]) < close_dis)
    
    if (length(close_point_indices_x) < 2 || length(close_point_indices_y) < 2) {
      indices[i] <- FALSE
    } else {
      close_point_x <- coords[close_point_indices_x, ]
      close_point_y <- coords[close_point_indices_y, ]
      
      in_y_range <- (coord_mesh[i, 2] < max(close_point_x[, 2]) - distance_to_boundary) & 
                    (coord_mesh[i, 2] > min(close_point_x[, 2]) + distance_to_boundary)
      
      in_x_range <- (coord_mesh[i, 1] < max(close_point_y[, 1]) - distance_to_boundary) & 
                    (coord_mesh[i, 1] > min(close_point_y[, 1]) + distance_to_boundary)
      
      indices[i] <- in_y_range & in_x_range
    }
  }
  
  return(indices)
}

gradient_plot = function(gradient_ori,gradient_mesh,quantile=0.9,size_ori=1,size_mesh=1){
    gradient_ori = gradient_ori %>% as.data.frame(check.names=F)%>% mutate(L2 = (g1^2+g2^2)^0.5)
    gradient_mesh = gradient_mesh %>% as.data.frame(check.names=F)%>% mutate(L2 = (g1^2+g2^2)^0.5)

k1=ggplot(gradient_ori)+geom_point(aes(x=s1,y=s2,color=y))+
    coord_fixed()+#scale_color_gradient2(low='#1F92BF',mid='white',high='#F24C3D',name='')+
    scale_color_gradientn(colours = c('#1F92BF','white','#F24C3D'),name='') + 
    #scale_color_viridis_c(name='')+
     scale_y_reverse()+ theme_article()+
    theme(legend.position = 'bottom',axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())

k2=ggplot(gradient_mesh)+geom_point(aes(x=s1,y=s2,color=pred))+
    coord_fixed()+#scale_color_gradient2(low='#1F92BF',mid='white',high='#F24C3D',name='')+
    scale_color_gradientn(colours = c('#1F92BF','white','#F24C3D'),name='') + 
    #scale_color_viridis_c(name='')+
     scale_y_reverse()+ theme_article()+
    theme(legend.position = 'bottom',axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())

pal =paletteer_d("khroma::sunset")
k3=ggplot()+geom_point(aes(x=s1,y=s2),data=gradient_mesh,color='grey95',size=1)+
    geom_point(aes(x=s1,y=s2,color=L2),data=gradient_mesh %>% filter(L2 > quantile(L2,quantile)),size=1)+
    coord_fixed()+#scale_color_gradient2(low='#1F92BF',mid='white',high='#F24C3D',name='')+
    scale_color_gradientn(colours = pal,name='',limits=c(min(gradient_mesh$L2),max(gradient_mesh$L2))) + 
    #scale_colour_ghibli_c("LaputaMedium")+
    #scale_color_viridis_c(name='',option='plasma')+
     scale_y_reverse()+ theme_article()+
    theme(legend.position = 'bottom',axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())

return(plot_grid(k1,k2,k3,nrow=1))

}

curvature_plot = function(gradient_mesh,quantile=0.9,size_ori=1,size_mesh=1){

    gradient_mesh = gradient_mesh %>% as.data.frame(check.names=F)%>% mutate(L2 = (g1^2+g2^2)^0.5,Laplacian=c11+c22)

k2=ggplot(gradient_mesh)+geom_point(aes(x=s1,y=s2,color=pred))+
    coord_fixed()+#scale_color_gradient2(low='#1F92BF',mid='white',high='#F24C3D',name='')+
    scale_color_gradientn(colours = c('#1F92BF','white','#F24C3D'),name='') + 
    #scale_color_viridis_c(name='')+
     scale_y_reverse()+ theme_article()+
    theme(legend.position = 'bottom',axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())

pal =paletteer_d("khroma::sunset")
k3=ggplot()+geom_point(aes(x=s1,y=s2),data=gradient_mesh,color='grey95',size=1)+
    geom_point(aes(x=s1,y=s2,color=L2),data=gradient_mesh %>% filter(L2 > quantile(L2,quantile)),size=1)+
    coord_fixed()+#scale_color_gradient2(low='#1F92BF',mid='white',high='#F24C3D',name='')+
    scale_color_gradientn(colours = pal,name='',limits=c(min(gradient_mesh$L2),max(gradient_mesh$L2))) + 
    #scale_colour_ghibli_c("LaputaMedium")+
    #scale_color_viridis_c(name='',option='plasma')+
     scale_y_reverse()+ theme_article()+
    theme(legend.position = 'bottom',axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())

k4=ggplot()+geom_point(aes(x=s1,y=s2),data=gradient_mesh,color='grey95',size=1)+
    geom_point(aes(x=s1,y=s2,color=Laplacian),data=gradient_mesh %>% filter(L2 > quantile(L2,quantile)),size=1)+
    coord_fixed()+#scale_color_gradient2(low='#1F92BF',mid='white',high='#F24C3D',name='')+
    scale_color_gradientn(colours = pal,name='',limits=c(min(gradient_mesh$L2),max(gradient_mesh$L2))) + 
    #scale_colour_ghibli_c("LaputaMedium")+
    #scale_color_viridis_c(name='',option='plasma')+
     scale_y_reverse()+ theme_article()+
    theme(legend.position = 'bottom',axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.x = element_blank(),axis.ticks.y = element_blank())

return(plot_grid(k2,k3,k4,nrow=1))

}

wombling = function(model,boundary,h,threads,path="./",prefix = 'boundary0'){
    boundary <- as.matrix(boundary)

    start_point = boundary[-nrow(boundary),]
    end_point = boundary[-1,]
    mid_point = (start_point+end_point)/2
    u = end_point - start_point
    u_len = sqrt(u[,1]^2+u[,2]^2)
    u = u/u_len
    v = cbind(-u[,2],u[,1])

    f_a = predict(model, X.0 = matrix(1.0,nr=nrow(start_point)), coords.0 = start_point, n.omp.threads=threads)
    f_b = predict(model, X.0 = matrix(1.0,nr=nrow(end_point)), coords.0 = end_point, n.omp.threads=threads)
    f_mid = predict(model, X.0 = matrix(1.0,nr=nrow(mid_point)), coords.0 = mid_point, n.omp.threads=threads)
    f_ahv = predict(model, X.0 = matrix(1.0,nr=nrow(start_point)), coords.0 = start_point+h*v, n.omp.threads=threads)
    #f_nahv = predict(model, X.0 = matrix(1.0,nr=nrow(start_point)), coords.0 = start_point-h*v, n.omp.threads=threads)
    f_bhv = predict(model, X.0 = matrix(1.0,nr=nrow(end_point)), coords.0 = end_point+h*v, n.omp.threads=threads)
    #f_nbhv = predict(model, X.0 = matrix(1.0,nr=nrow(end_point)), coords.0 = end_point-h*v, n.omp.threads=threads)
    f_midhv = predict(model, X.0 = matrix(1.0,nr=nrow(mid_point)), coords.0 = mid_point+h*v, n.omp.threads=threads)
    #f_nmidhv = predict(model, X.0 = matrix(1.0,nr=nrow(mid_point)), coords.0 = mid_point-h*v, n.omp.threads=threads)

    f_a = f_a$p.y.0
    f_b = f_b$p.y.0
    f_mid = f_mid$p.y.0
    f_ahv = f_ahv$p.y.0
    #f_nahv = f_nahv$p.y.0
    f_bhv = f_bhv$p.y.0
    #f_nbhv = f_nbhv$p.y.0
    f_midhv = f_midhv$p.y.0
    #f_nmidhv = f_nmidhv$p.y.0

    g_a = (f_ahv-f_a)/h
    g_b = (f_bhv-f_b)/h
    g_mid = (f_midhv-f_mid)/h
    g_simpson = u_len * (g_a+4*g_mid+g_b)/6
    fwrite(data.table(g_simpson), paste0(path,"/",prefix,'_wombling_all.txt'), row.names=F, col.names=F, sep="\t", quote=F)

    g_a = apply(g_a, 1, mean)
    g_b = apply(g_b, 1, mean)
    g_mid = apply(g_mid, 1, mean)
    g_simpson = apply(g_simpson, 1, mean)


    boundary_result <- cbind(g_simpson, u_len,start_point,end_point)
    colnames(boundary_result) = c('g_simpson','u_len','start_x','start_y','end_x','end_y')
    fwrite(boundary_result, paste0(path,"/",prefix,'_wombling.txt'), row.names=F, col.names=T, sep="\t", quote=F)
    return(boundary_result)
}


format_principle_curve = function(boundary){
  colnames(boundary) = c('start_x','end_x','start_y','end_y')
  boundary$start = paste0(round(boundary$start_x,7),',',round(boundary$start_y,7))
  boundary$end = paste0(round(boundary$end_x,7),',',round(boundary$end_y,7))


  if(sum(table(c(boundary$start,boundary$end))==1)!=2){
    print('error')
    return(NA)
    }else{

        start_point = names(which(table(c(boundary$start,boundary$end))==1))[1]
        organize_line = start_point
        while(start_point!=(names(which(table(c(boundary$start,boundary$end))==1))[2])){
        next_point = c(boundary$end[which(boundary$start==start_point)],boundary$start[which(boundary$end==start_point)])
        if(length(next_point)==1){
            organize_line = c(organize_line,next_point)
        }else{
            next_point = next_point[which(! next_point%in% organize_line)]
            organize_line = c(organize_line,next_point)
        }
        start_point = next_point
        }

        organize_line=strsplit(organize_line,split=',') %>% unlist() %>% as.numeric() %>% matrix(ncol=2,byrow=T)
        return(organize_line)
    }

}

refine_boundary = function(boundary,min_dis){
    boundary = as.matrix(boundary)
    refined_line = boundary[1,]
    for(i in 1:(nrow(boundary)-1)){
        dis_pair = sqrt((boundary[i,1]-boundary[i+1,1])^2+(boundary[i,2]-boundary[i+1,2])^2)
        if(dis_pair>min_dis){
            add_point_num = floor(dis_pair/min_dis)
            add_point_x = seq(boundary[i,1],boundary[i+1,1],length.out = add_point_num+2)[2:(add_point_num+2)]
            add_point_y = seq(boundary[i,2],boundary[i+1,2],length.out = add_point_num+2)[2:(add_point_num+2)]
            refined_line = rbind(refined_line,cbind(add_point_x,add_point_y))
        }else{
            refined_line = rbind(refined_line,cbind(boundary[i+1,1],boundary[i+1,2]))
        }
    }
    rownames(refined_line) = NULL
    colnames(refined_line) = NULL
    return(refined_line)
}