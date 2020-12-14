\name{init.cfust}
\alias{init.cfust}
\alias{init.fmcfust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialization for Fitting Finite Mixtures of Canonical Fundamental Skew t-Distributions
}
\description{
Computes different sets of initial values for finite mixtures of canonical fundamental skew t (FM-CFUST) model based on an initial clustering, transformation approiach, moment-based approach, or nested-model appraoch.
}
\usage{
init.cfust(g, dat, q=p, initial=NULL, known=NULL, clust=NULL, nkmeans=20, 
    method=c("moments","transformation","EMMIXskew","EMMIXuskew"))
init.fmcfust(g, dat, q=p, initial=NULL, known=NULL, clust=NULL, nkmeans=20, 
    method=c("moments","transformation","EMMIXskew","EMMIXuskew"))
}
\arguments{
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
  \item{clust}{
  (optional) a numeric value of length \code{nrow(dat)} containing the 
  initial labels for each data point in \code{dat}. The default is \code{NULL},
  indicating no initial clustering is known.
}
  \item{nkmeans}{
  (optional) a numeric value indicating how many k-means trials to be used 
  when searching for initial values. The default is \code{20}.
}
  \item{method}{
  (optional) a string indicating which method to use to generate initial values. 
  See Details.
}
}
\details{
As the EM algorithm is sensitive to the starting value, 
it is highly recommended to apply a wide range different initializations.
To obtain different sets of starting values using the strategy described in 
Section 5.1.3 of Lee and McLachlan (2014), \code{init.cfust()} can be used, 
which will return a list of objects with the same structure as \code{initial}.
An example is given in the examples section below.
                                                          
The argument \code{known}, if specified, is a list structure containing
at least one of \code{mu}, \code{sigma}, \code{delta}, \code{dof}, \code{pro} 
(See \code{\link{dfmcfust}} for the structure of each of these elements). 
Note that although not all parameters need to be provided in \code{known}, 
the parameters that are provided must be fully specified. 
They cannot be partially specified, e.g. only some elements or some components are specified. 
}
\value{
a list object containing the following parameters:
  \item{mu}{
    a list of \code{g} numeric matrices containing the location parameter for each component.
}
  \item{sigma}{
    a list of \code{g} numeric matrices containing the scale parameter for each component.   
}
  \item{delta}{
    a list of \code{g} numeric matrices containing the skewness matrix for each component.
}
  \item{dof}{
    a numeric vector of length \code{g} representing the degrees of freedom for each component.
}
  \item{pro}{
    a vector of length of \code{g} specifying the mixing proportions for each component.
}
  \item{tau}{
   an \code{g} by \code{n} matrix of initial probability of component membership. 
}
  \item{clusters}{
   a vector of length {n} of initial partition.  
}
  \item{loglik}{
   the initial log likelihood value.    
}
}
\references{
Lee S.X. and McLachlan, G.J. (2016). Finite mixtures of canonical fundamental 
skew t-distributions: the unification of the restricted and unrestricted 
skew t-mixture models. Statistics and Computing 26, 573-589.

Lee S.X. and McLachlan, G.J. (2017). EMMIXcskew: An R Package for the 
Fitting of a Mixture of Canonical Fundamental Skew t-Distributions. 
Journal of Statistical Software 83(3), 1-32. 
URL \url{10.18637/jss.v083.i03}.  
   
}
\seealso{
\code{\link{rfmcfust}}, \code{\link{dfmcfust}}, \code{\link{fmcfust.contour.2d}}
}
\examples{

	#a short demo using geyser data
	library(MASS)
	data(geyser)
	initial.transformation <- init.cfust(3, geyser, method="transformation") 
	initial.transformation$loglik
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM algorithm }
\keyword{ multivariate skew t distribution}% __ONLY ONE__ keyword per line
\keyword{ maximum likelihood estimation}
