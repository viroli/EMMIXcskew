\name{fmcfust.contour.2d}
\alias{fmcfust.contour.2d}
\alias{fmcfust.contour.3d}  
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
2D and 3D Visualisation of Fitted Contours
}
\description{
Create 2D or 3D contour plot.
}
\usage{   
fmcfust.contour.2d(dat, model, grid=50, drawpoints=TRUE, clusters=NULL,
    nlevels=10, map=c("scatter", "heat", "cluster"), 
    component=NULL, xlim, ylim, xlab, ylab, main, pcol=NULL, ccol=NULL, ...)     
fmcfust.contour.3d(dat, model, grid=20, drawpoints=TRUE, levels=0.9, 
    clusters=NULL, xlim, ylim, zlim, xlab, ylab, zlab, main, component=NULL, 
    pcol=NULL, ccol=NULL, ...)
}
\arguments{
  \item{dat}{
    the data matrix giving the coordinates of the point(s) where the density is evaluated.
    This must be a matrix with at least \code{2} columns for \code{fmcfust.contour.2d} 
    or \code{3} columns for \code{fmcfust.contour.3d}.
    If \code{dat} is not provided, then \code{xlim}, \code{ylim} and \code{zlim}
    must be provided, and \code{drawpoints} must be set to \code{FALSE}.
}
  \item{model}{
    a list containing the parameters of the model and also a vector of cluster labels for \code{dat}. 
    This is typically an output from \code{\link{fmcfust}}, 
    containing \code{mu}, \code{sigma}, \code{delta}, \code{dof}, 
    \code{pro} and \code{clusters}; see \code{\link{fmcfust}} for structure of \code{model}.
}
  \item{grid}{
  a positive integer specifying the grid size used to calculate the density map.
}
  \item{drawpoints}{
    logical. Points are plotted if \code{TRUE}.
}
  \item{clusters}{
  a vector of cluster labels to be applied when colouring the points. 
  This only applies when \code{drawpoints} is \code{TRUE}.
}
  \item{nlevels}{
  a positive integer specifying the number of contour levels 
}  
  \item{levels}{
  either a positive integer or a numeric vector specifying the contour levels to be drawn for each component
}
  \item{map}{
   character string specifying how to plot the points if \code{drawpoints=TRUE}.
   Possible values are \code{"scatter"} (default), \code{"heat"} and \code{"cluster"}.
   See the 'Details' section.
}
  \item{component}{
    the index of the components to be plotted. See the 'Details' section.
}
  \item{xlim, ylim, zlim}{
    x-, y- and z- limits for the plot
}
  \item{xlab, ylab, zlab}{
    labels for x-, y-  and z- axis
}
  \item{main}{
    title of the plot
}  
  \item{pcol}{
    the color(s) to be used for plotted points
}
  \item{ccol}{
    the color(s) to be used for plotted contours
}
  \item{\dots}{
    additional arguments to \code{\link{plot.default}}
}
}
\details{
\code{fmcfust.contour.2d} draw contour plots for bivariate densities. 
The argument \code{dat} must be provided and must contain at least \code{2} columns.
Note that only the first two columns of \code{dat} will be used if \code{dat} have more than \code{2} columns.
For bivariate datset, the data points can be drawn as a scatter plot by specifying \code{map="scatter"} (default),
or as an intensity plot (\code{map="heat"}). 
Alternatively, a cluster map can be drawn instead (\code{map="cluster"}). 
Note that if an intensity plot is used, the data points will not be drawn, that is,
\code{drawpoints} will be set to \code{FALSE}.

The argument \code{component} specifies which individual component is drawn.
When \code{component=FALSE}, the mixture contour is drawn. 
If specified, \code{component} is a integer vector of the index of the components
to be drawn. It can only take values between \code{1} an \code{g} inclusive.
For example, \code{component=c(1,3)} will draw the first and third component contours.

If the argument \code{model} contains the cluster labels (\code{model$clusters}), 
the data point will be coloured according to their cluster.  
}
\seealso{
  \code{\link{fmcfust}}, \code{\link{contour}}
}
\examples{          
\donttest{
#2D plots
data(iris)
iris.versicolor <- iris[iris$Species=="versicolor",2:3]
Fit.versicolor <- fmcfust(1, iris.versicolor)
fmcfust.contour.2d(iris.versicolor, Fit.versicolor, drawpoints=FALSE, main="versicolor")

#3D plot
obj <- list()
obj$mu <- list(matrix(c(0,0,0),3), matrix(c(5,5,5),3))
obj$sigma <- list(matrix(c(5,2,1,2,5,1,1,1,1),3,3), 2*diag(3))
obj$delta <- list(matrix(c(1,0,0,1,0,0,1,0,0),3,3), matrix(c(5,0,0,0,10,0,0,0,15),3,3))
obj$dof <- c(3,3)
obj$pro <- c(0.2, 0.8)
fmcfust.contour.3d(model=obj, level=0.98, drawpoints=FALSE, xlab="X", ylab="Y", zlab="Z", 
    xlim = c(-20, 50), ylim = c(-20, 50), zlim = c(-20, 80), ccol = "green")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ contour }
\keyword{ 3d }% __ONLY ONE__ keyword per line
