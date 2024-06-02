#' Fit NNGP model
#'
#' @param coords coordinates of the data points, each row is a data point, sample * 2 matrix
#' @param y response variable (e.g. gene expression, cell type proportion), sample * 1 matrix
#' @param neighbor number of neighbors used in NNGP model
#' @param threads number of threads used in NNGP model
#' @param n.samples number of samples used in posterior sampling
#' @param verbose whether to print out the progress
#' @param cov.model covariance model used in NNGP model, default is "matern"
#' @param starting starting values for the parameters used in NNGP model
#' @param tuning tuning parameters for the parameters used in NNGP model
#' @param priors prior distributions for the parameters used in NNGP model
#' @param n.report report at every n.report iterations
#' @returns a fitted NNGP model
#'
#' @import spNNGP


fit_NNGP = function(coords,y,neighbor=10,threads=1,n.samples=500,verbose=FALSE,
                    cov.model = "matern", # here we use mater kernal
                    starting = list("phi"=6, "sigma.sq"=5, "tau.sq"=1,"nu" = 2.5),
                    tuning = list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5,"nu"=0.5),
                    priors = list("phi.Unif"=c(0.01,300), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1),"nu.Unif"=c(2.00001,5)),
                    n.report = 100
                    ){
  # fit NNGP
  x=matrix(1,nr=nrow(coords)) # intercept
  m.r <- spNNGP::spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=neighbor,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads= threads , n.report=n.report,verbose=verbose)
  return(m.r)
}
