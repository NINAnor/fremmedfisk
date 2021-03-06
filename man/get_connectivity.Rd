% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataIO.R
\name{get_connectivity}
\alias{get_connectivity}
\title{Download and store connectivity matrix from the NOFA database.}
\usage{
get_connectivity(serveradress, datafolder, tablename)
}
\arguments{
\item{serveradress}{URL of the server wher the database is located.}

\item{datafolder}{Relative or absolute path of the folder where the data should be stored locally.}

\item{tablename}{the name of the db table holding the connectivity matrix subset in quesion.}
}
\value{
A data.frame stored as .rds object named view_occurrence_by_event.rds
}
\description{
Download and store connectivity matrix from the NOFA database.
}
\examples{
get_connectivity(serveradress="my-server.no",datafolder="./data",tablename="temporary.connectivity_troendelag")

}
