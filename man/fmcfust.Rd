\name{fmcfust}
\alias{fmcfust}
\alias{summary.fmcfust}
\alias{print.fmcfust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Finite Mixtures of Multivariate Canonical Fundamental Skew t-Distributions
}
\description{
Computes maximum likelihood estimators (MLE) for finite mixtures of canonical fundamental 
multivariate skew t (FM-CFUST) model via the EM algorithm.
}
\usage{
fmcfust(g=1, dat, q, initial=NULL, known=NULL, clust=NULL, itmax=100, eps=1e-6, 
    nkmeans=20, verbose=T, method=c("moments","transformation","EMMIXskew","EMMIXuskew"), 
    convergence=c("Aitken","likelihood","parameters"))
\method{summary}{fmcfust}(object, ...)
\method{print}{fmcfust}(x, ...)                     
}
\arguments{
  \item{object, x}{
    an object class of class \code{"fmcfust"}, i.e. a fitted model.
}
  \item{g}{
    a scalar specifying the number of components in the mixture model
}
  \item{dat}{
    the data matrix giving the coordinates of the point(s) where the density is evaluated.
    This is either a vector of length \code{p} or a matrix with \code{p} columns.
}   
  \item{q}{
    a scalar specifying how many number of columns the skewness matrix \code{delta} has.
}
  \item{initial}{
  (optional) a list containing the initial parameters of the mixture model. 
  See the 'Details' section. The default is \code{NULL}. 
}
  \item{known}{
  (optional) a list containing parameters of the mixture model that are known
  and not required to be estimated. See the 'Details' section. The default is \code{NULL}.
}
  \item{itmax}{
  (optional) a positive integer specifying the maximum number of EM iterations 
  to perform. The default is \code{100}.
}
  \item{eps}{
  (optional) a numeric value used to control the termination criteria for the EM loops.
  It is the maximum tolerance for the absolute difference between the log-likelihood value
  and the asymptotic log likelihood value. The default is \code{1e-6}.
}   
  \item{clust}{
  (optional) a numeric value of length \code{nrow(dat)} containing the 
  initial labels for each data point in \code{dat}. The default is \code{NULL},
  indicating no initial clustering is known.
}
  \item{nkmeans}{
  (optional) a numeric value indicating how many k-means trials to be used 
  when searching for initial values. The default is \code{20}.
}
  \item{verbose}{
  (optional) a logical value. If \code{TRUE}, output for each iteration will be printed out.
  if \code{FALSE}, no output is printed. The default is \code{TRUE}. See the 'Details' section. 
}    
  \item{method}{
  (optional) a string indicating which method to use to generate initial values. 
  See \code{\link{init.cfust}}.
}    
  \item{convergence}{
  (optional) a string indicating which convergence criterion to use to terminate the iterations. 
  The default \code{"Aitken"} uses Aitken accelation, 
  whereas \code{"likelihood"} uses the relative difference in log likleihood value, 
  and \code{"parameters"} checks the changes in parameter estimates. 
}
  \item{\dots}{ not used. }
}
\details{
The arguments \code{init} and \code{known}, if specified, is a list structure containing
at least one of \code{mu}, \code{sigma}, \code{delta}, \code{dof}, \code{pro} 
(See \code{\link{dfmcfust}} for the structure of each of these elements). 
If \code{init=FALSE} (default), the program uses an automatic approach based on 
moments estimate and k-means clustering to generate an initial value for the model parameters.
Note that this may not provide the best results.

As the EM algorithm is sensitive to the starting value, 
it is highly recommended to apply a wide range different initializations.
Some strategies are implemented in \code{\link{init.cfust}}.  
}
\value{
  \item{mu}{
    a list of \code{g} numeric matrices containing the location parameter for each component.
}
  \item{sigma}{
    a list of \code{g} numeric matrices containing the scale parameter for each component.   
}
  \item{delta}{
    a list of \code{g} numeric matrices containing the skewness parameter for each component.
}
  \item{dof}{
    a numeric vector of length \code{g} representing the degrees of freedom for each component.
}
  \item{pro}{
    a vector of length of \code{g} specifying the mixing proportions for each component.
}
  \item{tau}{
   an \code{g} by \code{n} matrix of posterior probability of component membership. 
}
  \item{clusters}{
   a vector of length {n} of final partition.  
}
  \item{loglik}{
   the final log likelihood value.    
}
  \item{lk}{
   a vector of log likelihood values at each EM iteration.
}
  \item{iter}{
   number of iterations performed.
}
  \item{eps}{
   the final absolute difference between the log likelihood value and the asymptotic log likelihood value.
} 
  \item{aic, bic}{
   Akaike Information Criterion (AIC), Bayes Information Criterion (BIC)
}
}
\references{
Lee S.X. and McLachlan, G.J. (2016). Finite mixtures of canonical fundamental 
skew t-distributions: the unification of the restricted and unrestricted 
skew t-mixture models. Statistics and Computing 26, 573-589

Lee S.X. and McLachlan, G.J. (2017). EMMIXcskew: An R Package for the 
Fitting of a Mixture of Canonical Fundamental Skew t-Distributions. 
Journal of Statistical Software 83(3), 1-32. 
URL \url{10.18637/jss.v083.i03}.   
 
}
\seealso{
\code{\link{init.cfust}}, \code{\link{rfmcfust}}, \code{\link{dfmcfust}}, \code{\link{fmcfust.contour.2d}}
}
\examples{   
\donttest{                
#a short demo using geyser data
library(MASS)
Fit <- fmcfust(3, geyser)
summary(Fit)
print(Fit)   
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ multivariate skew t distribution}% __ONLY ONE__ keyword per line
\keyword{ maximum likelihood estimation}
