#.libPaths('/proj/yunligrp/users/jiawen/software/Rlibs')
library(spNNGP)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggsci)
library(egg)
library(paletteer)

source('./Laputa.R')

rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

## simulate data
set.seed(1)
N <- 100
tau <- 1
# synthetic location
#coords <- matrix(runif(2*N),nr=N,nc=2) 

coords_x = seq(0,1,0.05)
coords_y = seq(0,1,0.05)
coords = expand.grid(coords_x,coords_y) %>% as.matrix()

y <- rnorm(nrow(coords),10*(sin(3*pi*coords[,1])+cos(3*pi*coords[,2])),tau)
# y <- rnorm(N,10*(sin(3*pi*coords[,1])*cos(3*pi*coords[,2])),tau)
y=y # response
x=matrix(1,nr=nrow(coords)) # intercept just use 1 here


## start fit NNGP model: m.r
starting <- list("phi"=6, "sigma.sq"=5, "tau.sq"=1,"nu" = 2.5)
tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5,"nu"=0.5)
priors <- list("phi.Unif"=c(0.01,300), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1),"nu.Unif"=c(2,3))
cov.model <- "matern"
n.report <- 100
n.samples = 500

m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads=10, n.report=n.report)


# create a grid with higher resolution
coord_mesh_x = seq(min(coords[,1]),max(coords[,1]),length.out=50)
coord_mesh_y = seq(min(coords[,2]),max(coords[,2]),length.out=50)
coord_mesh = expand.grid(coord_mesh_x,coord_mesh_y) %>% as.matrix()

# define step to use in finite difference, here I use minimal sepatation * 0.8
min_sep = 0.05
path='./'

# start estimating gradients for points in coords
gradient_all = finite_difference(coords,min_sep*0.8,m.r,threads=thread,prefix = 'ori',path=path)
gradient_all = cbind(coords,y,gradient_all)
colnames(gradient_all) = c('s1','s2','y','pred','g1','g2','g1_min','g1_max','g2_min','g2_max')
fwrite(gradient_all,paste0(path,'/gradient_ori.txt'),row.names=F,col.names=T,sep="\t",quote=F)

# start estimating gradients for points in coords_mesh
gradient_mesh = finite_difference(coord_mesh,min_sep*0.8,m.r,threads=thread,prefix = 'mesh',path=path)
gradient_mesh = cbind(coord_mesh,gradient_mesh)
colnames(gradient_mesh) = c('s1','s2','pred','g1','g2','g1_min','g1_max','g2_min','g2_max')
fwrite(gradient_mesh,paste0(path,'/gradient_mesh.txt'),row.names=F,col.names=T,sep="\t",quote=F)

# plot
ggsave(gradient_plot(gradient_all,gradient_mesh),filename = paste0(path,'/gradient_plot.png'),width=10,height=5,bg='white',dpi=300)