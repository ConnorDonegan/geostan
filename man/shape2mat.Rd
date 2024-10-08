% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convenience-functions.R
\name{shape2mat}
\alias{shape2mat}
\title{Create spatial and space-time connectivity matrices}
\source{
Bivand, Roger S. and Pebesma, Edzer and Gomez-Rubio, Virgilio (2013). Applied spatial data analysis with R, Second edition. Springer, NY. https://asdar-book.org/

Griffith, Daniel A. (2012). Space, time, and space-time eigenvector filter specifications that account for autocorrelation. Estadística Espanola, 54(177), 7-34.

Haining, Robert P. and Li, Guangquan (2020). Modelling Spatial and Spatial-Temporal Data: A Bayesian Approach. CRC Press.
}
\usage{
shape2mat(
  shape,
  style = c("B", "W"),
  queen,
  method = c("queen", "rook", "knn"),
  k = 1,
  longlat = NULL,
  snap = sqrt(.Machine$double.eps),
  t = 1,
  st.style = c("contemp", "lag"),
  quiet = FALSE
)
}
\arguments{
\item{shape}{An object of class \code{sf}, \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}.}

\item{style}{What kind of coding scheme should be used to create the spatial connectivity matrix? Defaults to "B" for binary; use "W" for row-standardized weights.}

\item{queen}{Deprecated: use the `method' argument instead. This option is passed to \code{\link[spdep]{poly2nb}} to set the contiguity condition. Defaults to \code{TRUE} so that a single shared boundary point (rather than a shared border/line) between polygons is sufficient for them to be considered neighbors.}

\item{method}{Method for determining neighbors: queen, rook, or k-nearest neighbors. See Details for more information.}

\item{k}{Number of neighbors to select for k-nearest neighbor method. Passed to \code{spdep::knearneigh}.}

\item{longlat}{If longlat = TRUE, Great Circle (rather than Euclidean) distances are used; great circle circle distances account for curvature of the Earth.}

\item{snap}{Passed to \code{spdep::poly2nb}; "boundary points less than ‘snap’ distance apart are considered to indicate contiguity."}

\item{t}{Number of time periods. Only the binary coding scheme is available for space-time connectivity matrices.}

\item{st.style}{For space-time data, what type of space-time connectivity structure should be used? Options are "lag" for the lagged specification and "contemp" (the default) for contemporaneous specification (see Details).}

\item{quiet}{If \code{TRUE}, messages will be silenced.}
}
\value{
A spatial connectivity matrix in sparse matrix format. Binary matrices are of class \code{ngCMatrix}, row-standardized are of class \code{dgCMatrix}, created by \code{\link[Matrix]{sparseMatrix}}.
}
\description{
Creates sparse matrix representations of spatial connectivity structures
}
\details{
The method argument currently has three options. The queen contiguity condition defines neighbors as polygons that share at least one point with one another. The rook condition requires that they share a line or border with one another. K-nearest neighbors is based on distance between centroids. All methods are implemented using the spdep package and then converted to sparse matrix format.

Haining and Li (Ch. 4) provide a helpful discussion of spatial connectivity matrices (Ch. 4).

The space-time connectivity matrix can be used for eigenvector space-time filtering (\code{\link[geostan]{stan_esf}}. The 'lagged' space-time structure connects each observation to its own past (one period lagged) value and the past value of its neighbors. The 'contemporaneous' specification links each observation to its neighbors and to its own in situ past (one period lagged) value (Griffith 2012, p. 23).
}
\examples{
data(georgia)

## binary adjacency matrix
C <- shape2mat(georgia, "B", method = 'rook')

## number of neighbors per observation
summary( n_nbs(C) )
head(Matrix::summary(C))

## row-standardized matrix 
W <- shape2mat(georgia, "W", method = 'rook')

## summary of weights
E <- edges(W, unique_pairs_only = FALSE)
summary(E$weight)

## space-time matricies 
## for eigenvector space-time filtering
## if you have multiple years with same geometry/geography,
## provide the geometry (for a single year!) and number of years \code{t}
Cst <- shape2mat(georgia, t = 5)
dim(Cst)
EVst <- make_EV(Cst)
dim(EVst)

}
\seealso{
\code{\link[geostan]{edges}} \code{\link[geostan]{row_standardize}} \code{\link[geostan]{n_nbs}}
}
