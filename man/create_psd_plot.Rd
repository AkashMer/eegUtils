% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frequency_plotting.R
\name{create_psd_plot}
\alias{create_psd_plot}
\title{Create a PSD plot}
\usage{
create_psd_plot(psd_out, freq_range, topo_freq, chan_info)
}
\arguments{
\item{psd_out}{PSD to plot.}

\item{freq_range}{Frequency range to plot.}

\item{topo_freq}{Display a topographic distribution of power spectral density
at particular frequency/frequencies.}

\item{chan_info}{Channel info of the data}
}
\value{
ggplot showing power spectral density.
}
\description{
Create a PSD plot
}
\keyword{internal}
