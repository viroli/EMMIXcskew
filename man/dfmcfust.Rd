\name{dfmcfust}
\alias{dfmcfust}
\alias{dcfust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate canonical fundamental skew t-distribution
}
\description{
The probability density function for the multivariate canonical fundamental 
skew t (CFUST) distribution and finite mixture of CFUST distributions
}
\usage{
dfmcfust(dat, mu=NULL, sigma=NULL, delta=NULL, dof=NULL, pro=NULL, known=NULL)
dcfust(dat, mu=NULL, sigma=NULL, delta=NULL, dof=1, known=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    the data matrix giving the coordinates of the point(s) where the density is evaluated.
    This is either a vector of length \code{p} or a matrix with \code{p} columns.
}
  \item{mu}{
    for \code{dcfust}, this is a numeric vector of length \code{p} representing the location parameter;
    for \code{dfmcfust}, this is list of \code{g} numeric matrices each having \code{p} rows and \code{1} column
    containing the location parameter for each component.
}
  \item{sigma}{
    for \code{dcfust}, this is a numeric positive definite matrix with dimension \code{(p,p)} representing the scale parameter;
    for \code{dfmcfust}, this is list of \code{g} numeric matrices containing the scale parameter for each component.   
}
  \item{delta}{
    for \code{dcfust}, this is a numeric matrix of size \code{p} by \code{q} representing the skewness matrix;
    for \code{rfmcfust}, this is list of \code{g} numeric matrices each having \code{p} rows and \code{q} column
    containing the skewness parameter for each component.
}
  \item{dof}{
    for \code{dcfust}, this is a positive integer specifying the degrees of freedom;
    for \code{dfmcfust}, this is numeric vector of length \code{g} representing the degrees of freedom for each component.
}
  \item{pro}{
    the mixing proportions; for\code{dcfust}, this is equal to \code{1};
    for \code{dfmcfust}, this is vector of length of \code{g}
    specifying the mixing proportions for each component.
}
  \item{known}{
    a list containing the parameters of the model. 
    If specified, it overwrites the values of \code{mu}, \code{sigma}, 
    \code{delta}, \code{dof} and \code{pro}.
} 
}
\details{
The function \code{dcfust} computes the density value of a specified multivariate 
canonical fundamental skew t (CFUST) distribution. 
If any model parameters are not specified, their default values are used:
\code{mu} and \code{delta} are zero vectors, \code{sigma} is the identity matrix,
and \code{dof} is \code{1}.
 
The function \code{dfmcfust} computes the density value for a specified mixture of
MST distribution. Note that \code{dfmcfust} expects at least \code{dof} is specified.
Other missing parameters will take the default value described above. 
When \code{g=1}, \code{dfmcfust} passes the call to \code{dcfust}. 
Model parameters can be passed to \code{dcfust} and \code{dfmcfust} 
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
\code{dcfust} and \code{dfmcfust} returns a numeric vector of density values
}
\references{
%% ~put references to the literature/web site here ~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{rcfust}}, \code{\link{rfmcfust}}
}
\examples{           
dcfust(c(1,2), mu=c(1,5), sigma=diag(2), delta=matrix(c(-3,1,1,1),2,2), dof=4)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mixture density }
\keyword{ multivariate skew t distribution}% __ONLY ONE__ keyword per line
\keyword{ multivariate distribution }
