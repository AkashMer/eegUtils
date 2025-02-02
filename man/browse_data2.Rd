% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_viewers2.R
\name{browse_data2}
\alias{browse_data2}
\alias{browse_data2.eeg_data}
\title{Browse continuous EEG data.}
\usage{
browse_data2(data, ...)

\method{browse_data2}{eeg_data}(data, sig_length = 5, n_elecs = NULL, downsample = TRUE, ...)
}
\arguments{
\item{data}{\code{eeg_data} object to be plotted.}

\item{...}{Other parameters passed to browsing functions.}

\item{sig_length}{Length of signal to be plotted initially. Default is 5 sec}

\item{n_elecs}{Number of electrodes to be plotted on a single screen. Can be
an integer indicating the number of electrodes to display or a character
vector specifying the electrodes. Electrodes not present in the data will
be ignored}

\item{downsample}{Reduces size of data by only plotting every 4th point,
speeding up plotting considerably. Defaults to TRUE.}
}
\description{
A Shiny app for scrolling through continuous EEG data interactively.
With EEG data (continuous), data can be viewed as a butterfly plot
(all electrodes overlaid) or as individual traces (electrodes "stacked").
The arrow buttons help to navigate the time-domain with single arrow allowing
1 sec steps and double arrows allowing step length equal to \code{sig_length}
}
\section{Methods (by class)}{
\itemize{
\item \code{browse_data2(eeg_data)}: Browse continuous EEG data.

}}
\author{
Akash Mer
}
