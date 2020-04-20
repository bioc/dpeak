\name{exportPeakList}
\alias{exportPeakList}
\alias{exportPeakList,DpeakFit-method}
\title{
Export deconvolution results to text files
}
\description{
Export deconvolution results to text files in TXT, BED, or GFF file formats.
}
\usage{
exportPeakList(object, ...)
\S4method{exportPeakList}{DpeakFit}( object, type=NA, filename=NA, ... )
}
\arguments{
    \item{object}{Object of class \code{DpeakFit},
        deconvolution model fits obtained using the method \code{dpeakFit}. }
    \item{type}{Format of the exported file.
        Possible values are \code{"txt"}, \code{"bed"}, and \code{"gff"}. See Details. }
    \item{filename}{Name of the exported file. }
    \item{...}{ Other parameters to be passed through to generic \code{exportPeakList}.}
}
\details{
Columns of TXT file format (\code{type="txt"}) include
chromosome, binding site, relative binding strength in each peak region,
and the peak region that each binding event belongs to.
\code{type="bed"} and \code{type="gff"} export deconvolution results
in standard BED and GFF file formats, respectively,
where score is the relative binding strength multiplied by 1000.
The feature of GFF file and the name of BED file are
the peak region that each binding event belongs to.
}
\value{Export deconvolution results to text files}
\author{ Dongjun Chung }
\seealso{
\code{\link{dpeakFit}}, \code{\linkS4class{DpeakFit}}.
}
\examples{
data(exampleData)
fit <- dpeakFit(exampleData)
exportPeakList( fit, type="txt", filename="result.txt" )
exportPeakList( fit, type="bed", filename="result.bed" )
exportPeakList( fit, type="gff", filename="result.gff" )
}
\keyword{models}
\keyword{methods}
