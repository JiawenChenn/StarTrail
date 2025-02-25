
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
