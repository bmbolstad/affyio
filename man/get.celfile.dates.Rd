\name{get.celfile.dates}
\alias{get.celfile.dates}
\title{Extract Dates from CEL files}
\description{This function reads the header information for a series of CEL files then extracts and returns the dates.}
\usage{
get.celfile.dates(filenames, ...)
}
\arguments{
  \item{filenames}{a vector of characters with the CEL filenames. May be fully pathed.}
  \item{\dots}{further arguments passed on to \code{\link{read.celfile.header}}.}
}
\details{
The function uses \code{\link{read.celfile.header}} to read in the header of each file. The \code{ScanDate} component is then parsed to extract the date.
Note that an assumption is made about the format. Namely, that dates are in the Y-m-d or m/d/y format.}
\value{
A vector of class \code{\link{Date}} with one date for each celfile.}
\author{
Rafael A. Irizarry <rafa@jimmy.harvard.edu>
}
\seealso{
See Also as \code{\link{read.celfile.header}}.
}
\keyword{IO}

