\name{rfmcfust}
\alias{rfmcfust}
\alias{rcfust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of Mixture Data
}
\description{
Generate random sample from a specified mixture of 
multivariate  canonical fundamental skew t distribution
}
\usage{
rfmcfust(g, n, mu, sigma, delta, dof=rep(10,g), pro=rep(1/g,g), known=NULL) 
rcfust(n=1, mu = NULL, sigma=NULL, delta=NULL, dof=1, known=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
    a scalar specifying the number of components in the mixture model
}
  \item{n}{
    either a positive integer specifying the total number of points to be generated
    or a vector (of length \code{g}) of positive integers specifying the number of points to be generated in each component.
}
  \item{mu}{
    for \code{rcfust}, this is a numeric vector of length \code{p} representing the location parameter;
    for \code{rfmcfust}, this is list of \code{g} numeric matrices each having \code{p} rows and \code{1} column
    containing the location parameter for each component.
}
  \item{sigma}{
    for \code{rcfust}, this is a numeric positive definite matrix with dimension \code{(p,p)} representing the scale parameter;
    for \code{rfmcfust}, this is list of \code{g} numeric matrices containing the scale parameter for each component.
}
  \item{delta}{
    for \code{rcfust}, this is a numeric matrix of size \code{p} by \code{q} representing the skewness matrix;
    for \code{rfmcfust}, this is list of \code{g} numeric matrices each having \code{p} rows and \code{q} column
    containing the skewness parameter for each component.
}
  \item{dof}{
    for \code{rcfust}, this is a positive integer specifying the degrees of freedom;
    for \code{rfmcfust}, this is numeric vector of length \code{g} representing the degrees of freedom for each component.
}
  \item{pro}{
    the mixing proportions; for \code{rcfust}, this is equal to \code{1};
    for \code{rfmcfust}, this is vector of length of \code{g}
    specifying the mixing proportions for each component.
}
  \item{known}{
    a list containing the parameters of the model. 
    If specified, it overwrites the values of \code{mu}, \code{sigma}, 
    \code{delta}, \code{dof} and \code{pro}.
}
}
\details{
\code{rcfust} generates a sample \code{n} multivariate CFUST observations.
\code{rfmcfust} generates a mixture of CFUST observation. 
Note that model parameters can be passed to \code{rcfust} and \code{rfmcfust} 
through the argument \code{known} or listed as individual arguments.
If both methods of input were used, the parameters specified in \code{known}
will be used.    
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{rcfust} returns an \code{n} by \code{p} numeric matrix of generated data.
\code{rfmcfust} returns an \code{n} by \code{p+1} numeric matrix of genertaed data.
The first \code{p} gives the coordinates of the generated data. 
The last column specifies which component each data point is generated from.
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
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dcfust}}, \code{\link{dfmcfust}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
rcfust(10,c(1,2),diag(2),matrix(c(2,1,1,2),2,2),4)

obj <- list()
obj$mu <- list(c(17,19), c(5,22), c(6,10))
obj$sigma <- list(diag(2), matrix(c(2,0,0,1),2), matrix(c(3,7,7,24),2))
obj$delta <- list(matrix(c(3,0,2,1.5),2,2), matrix(c(5,0,0,10),2,2), matrix(c(2,0,5,0),2,2))
obj$dof <- c(1, 2, 3)
obj$pro <- c(0.25, 0.25, 0.5)
rfmcfust(3, 100, known=obj)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ random number }
\keyword{ multivariate skew t distribution}% __ONLY ONE__ keyword per line
\keyword{ multivariate distribution }
