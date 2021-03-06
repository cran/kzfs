% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kzmd.R
\name{kzmd}
\alias{kzmd}
\title{Yet Another Multi-dimensional Kolmogorov-Zurbenko Filter}
\usage{
kzmd(ss, window, scale, k = 1, edges = TRUE)
}
\arguments{
\item{ss}{Data frame with value column behind time/space coordinates.}

\item{window}{Vector for window size of each dimension.}

\item{scale}{Vector for scale of each dimension.}

\item{k}{Iteration times of KZ filter. Defaults to 1.}

\item{edges}{Logic. Defaults to TRUE. FLASE means clear the data that
are located outside the time-space range of input data.}
}
\value{
Data framework with value column behind time/space coordinates.
}
\description{
This implement of spatial KZ-filter works for any dimensions. It is designed 
for cases with sparse data in large time-space.
}
\examples{
zs <- rbind(c(0,5,1,40),c(12,6,1,10),c(6,7,1,20),c(15,15,4,80))
colnames(zs) <- c("x","y","z","v")
zs <- kzmd(data.frame(zs), scale=c(1,1,1), window=c(3,5,3), k=4)
u <- zs[zs$z==1, -3]
x = sort(unique(u$x))
y = sort(unique(u$y))
z=df2mt(u, scale=c(1,1))	# Transfer from data frame to matrix.
image(x=x, y=y, z=z)

}
\seealso{
\code{\link[kza]{kz}}
}
\concept{Kolmogorov-Zurbenko filter}
\keyword{KZ-filter}
