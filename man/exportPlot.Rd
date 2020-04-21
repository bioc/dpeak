\name{exportPlot}
\alias{exportPlot}
\title{
Export plots to pdf files.
}
\description{
Exports the plots of estimated binding sites (\code{plotType="fit"})
or the goodness of fit (GOF) plots (\code{plotType="GOF"}) to a PDF file.
}
\usage{
exportPlot(x, y, ...)
}
\arguments{
    \item{x}{Object of class \code{DpeakFit}},
    \item{y}{Name of file to export to.}
    \item{...}{ Other parameters to be passed through to generic \code{exportPlot}.}
}
\details{
Exports the plots of estimated binding sites (\code{plotType="fit"})
or the goodness of fit (GOF) plots (\code{plotType="GOF"}) to a PDF file.
Its file name needs to be specified in the \code{filename} argument.
In both of these plots, estimated binding sites or simulated fragments are superimposed on
the plots of reads (or fragments) aligned to each position.
For SET data, if \code{plotType="fit"} and \code{strand=TRUE},
reads will be plotted in a strand-specific manner,
where each read is extended to \code{extension} from its 5' end.
If \code{smoothing=TRUE}, a smoothed plot (using the smoothing spline) is provided.
Unsmoothed plot is provided by default.
}
\value{Export plots to files}
\author{ Dongjun Chung }
\examples{
data(exampleData)
exampleFit <- dpeakFit( exampleData, maxComp=5)
exportPlot( exampleFit, filename="exampleResult_combined.pdf" )
}
\keyword{models}
\keyword{methods}
