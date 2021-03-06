% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{quantile_spike_and_slab}
\alias{quantile_spike_and_slab}
\title{Function for computing spike-and-slab quantiles.}
\usage{
quantile_spike_and_slab(p, w = 0.5, slab_precision = 1, spike_precision = 100)
}
\arguments{
\item{p}{Numeric probability value between 0 and 1.}

\item{w}{Numeric weight of spike-and-slab mixture between 0 (spike) and 1 (slab).}

\item{slab_precision}{Numeric precision of slab component.}

\item{spike_precision}{Numeric precision of spike component.}
}
\value{
Spike-and-slab quantile for given probability value.
}
\description{
Function for computing spike-and-slab quantiles.
}
\examples{
quantile_spike_and_slab(0.9)
quantile_spike_and_slab(0.5, w = 0.3)
}
