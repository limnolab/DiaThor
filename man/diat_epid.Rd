% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diat_epid.R
\encoding{UTF-8}
\name{diat_epid}
\alias{diat_epid}
\title{Calculates the EPID index (EPID)}
\usage{
diat_epid(resultLoad, maxDistTaxa = 2)
}
\arguments{
\item{resultLoad}{The resulting list obtained from the diat_loadData() function}

\item{maxDistTaxa}{Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2}
}
\description{
The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
so the acronyms and species' names are recognized
References for the index:
\itemize{
\item Dell'Uomo, A. (1996). Assessment of water quality of an Apennine river as a pilot study for diatom-based monitoring of Italian watercourses. Use of algae for monitoring rivers, 65-72.
}
Sample data in the examples is taken from:
\itemize{
\item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. <https://doi:10.1016/j.ecolind.2019.105951>
}
}
\examples{
\donttest{
# Example using sample data included in the package (sampleData):
data("diat_sampleData")
# First, the diat_loadData() function has to be called to read the data
# The data will be stored into a list (loadedData)
# And an output folder will be selected through a dialog box if resultsPath is empty
# In the example, a temporary directory will be used in resultsPath
df <- diat_loadData(diat_sampleData, resultsPath = tempdir())
epidResults <- diat_epid(df)
}
}
\keyword{bioindicator}
\keyword{biotic}
\keyword{diatom}
\keyword{ecology}
