#' Reformat the principle curve output
#'
#' @param boundary boundary coordinate n*4 matrix, output from elpigraph python package
#' @returns formatted boundary coordinate n*2 matrix


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