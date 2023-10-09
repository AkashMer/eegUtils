#' Identify globally bad channels
#'
#' Identify bad channels that should be excluded from ICA decomposition.
#' Channel autocorrelation values, variance and Hurst exponents are used to
#' determine bad channels.
#' *See Reference for more details*
#'
#' @author Akash Mer
#'
#' @param data An object of class `eeg_data`
#' @param sds Number of standard deviations of normalized data to consider as
#' acceptable range. Defaults to `3`.
#' @param exclude Channels to be ignored
#' @param ... Further parameters (tbd)
#' @return A character vector with channel names
#' @family Pre-ICA artifact rejection
#' @seealso [bad_epochs()]
#' @references
#' Nolan, Whelan & Reilly (2010). FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods.
#' @export

bad_chans <- function(data,
                      sds = 3,
                      exclude = NULL,
                      ...) {

  orig_chan_info <- channels(data)

  orig_names <- channel_names(data)

  data_chans <- orig_names[!(orig_names %in% data$reference$ref_chans)]
  if (!is.null(exclude)) {
    if (!is.character(exclude)) {
      message("Channels to be excluded should be a character vector")
    }
    message("Excluding channel(s):",
            paste(exclude, ""))
    data_chans <- data_chans[!(data_chans %in% exclude)]
  }
  data_mat <- data$signals[, data_chans]

  # Calculate channel auto correlation and normalize it
  chan_corrs <- scale(colMeans(abs(stats::cor(data_mat))))

  # Calculate variation in the signal and normalize it
  chan_vars <- scale(apply(data_mat,
                           2,
                           stats::var))

  # Calculate the hurst exponent and normalize it
  chan_hurst <- scale(quick_hurst(data_mat))

  bad_chans <- matrix(c(abs(chan_hurst) > sds,
                        abs(chan_vars) > sds,
                        abs(chan_corrs) > sds),
                      nrow = 3,
                      byrow = TRUE)
  bad_chans <- apply(bad_chans,
                     2,
                     any)

  bad_chan_n <- data_chans[bad_chans]
  bad_chan_n
}

#' Identify globally bad epochs
#'
#' Identify bad epochs that should be excluded from ICA decomposition.
#' Uses standard deviation to identify bad epochs
#' *See Reference for more details*
#'
#' @author Akash Mer
#'
#' @param data An object of class `eeg_epochs`
#' @param sds Number of standard deviations of normalized data to consider as
#' acceptable range. Defaults to `3`.
#' @param exclude Channels to be ignored
#' @param ... Further parameters (tbd)
#' @return A numeric vector with epoch numbers
#' @family Pre-ICA artifact rejection
#' @seealso [bad_chans()]
#' @references
#' Nolan, Whelan & Reilly (2010). FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods.
#' @export

bad_epochs <- function(data,
                       sds = 3,
                       exclude = NULL,
                       ...) {

  # Exclude the channels specified by the user
  orig_chan_info <- channels(data)

  orig_names <- channel_names(data)

  data_chans <- orig_names[!(orig_names %in% data$reference$ref_chans)]
  if (!is.null(exclude)) {
    if (!is.character(exclude)) {
      message("Channels to be excluded should be a character vector")
    }
    message("Excluding channel(s):",
            paste(exclude, ""))
    data_chans <- data_chans[!(data_chans %in% exclude)]
  }
  data$signals <- data$signals[, data_chans]
  data_mat <- data.table::as.data.table(data)

  # Calculate overall channel means
  chan_means <- data_mat[, lapply(.SD, mean), .SDcols = data_chans]

  # Calculate amplitude range within in each epoch
  epoch_range <- data_mat[, lapply(.SD, function(x) max(x) - min(x)),
                          .SDcols = data_chans,
                          by = epoch]
  epoch_range <- epoch_range[, .(Mean = rowMeans(.SD)), by = epoch]
  # Normalize and check which are 'sds' deviations away
  epoch_range <- abs(scale(epoch_range$Mean)) > sds

  # Calculate channel deviation within each epoch
  epoch_diffs <- data_mat[, lapply(.SD, mean),
                          .SDcols = data_chans,
                          by = epoch][, lapply(.SD, function(x) x - mean(x)),
                                      .SDcols = data_chans][ ,
                                                        .(Mean = rowMeans(.SD))]
  # Normalize and check which are 'sds' deviations away
  epoch_diffs <- abs(scale(epoch_diffs$Mean)) > sds

  # Calculate channel variance within each epoch
  epoch_vars <- data_mat[, lapply(.SD, var), .SDcols = data_chans,
                         by = epoch][, apply(.SD, 1, mean),
                                     .SDcols = data_chans]
  # Normalize and check which are 'sds' deviations away
  epoch_vars <- abs(scale(epoch_vars)) > sds

  # Combine the results of all 3 parameters
  bad_epochs <- matrix(c(rowSums(epoch_vars) > 0,
                         rowSums(epoch_range) > 0,
                         rowSums(epoch_diffs) > 0),
                       ncol = 3)
  bad_epochs <- apply(bad_epochs, 1, any)
  bad_epochs <- unique(data$timings$epoch)[bad_epochs]
  bad_epochs
}

#' Quickly calculate simple Hurst exponent for a matrix
#'
#' @param data matrix of EEG signals
#' @importFrom data.table data.table
#' @keywords internal

quick_hurst <- function(data) {
  n <- nrow(data)
  data <- data.table::data.table(data)
  dat_cumsum <- data[, lapply(.SD, cumsum)]
  rs <- dat_cumsum[, lapply(.SD, max)] - dat_cumsum[, lapply(.SD, min)]
  rs <- rs / data[, lapply(.SD, stats::sd)]#column_sd
  as.numeric(log(rs) / log(n))
}
