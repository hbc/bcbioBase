% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runDate.R
\name{runDate}
\alias{runDate}
\title{Run date}
\usage{
runDate(projectDir)
}
\arguments{
\item{projectDir}{\code{character(1)}.
Project directory path.}
}
\value{
\code{Date}.
}
\description{
Get the run date from the project directory.
}
\details{
Alternatively, can parse YAML data, but this approach is faster and simpler.
}
\note{
Updated 2023-09-21.
}
\examples{
runDate(projectDir = "2018-01-01_illumina_rnaseq")
}
\author{
Michael Steinbaugh
}
