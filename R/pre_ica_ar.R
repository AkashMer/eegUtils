#' Identify globally bad channels
#'
#' Identify bad channels that should be excluded from ICA decomposition.
#' Uses Hurst scores, variance and correlation to identify bad channels.
#' *See Reference for more details*
#'
#' @author Akash Mer
#'
#' @param data An object of class `eeg_data`
#' @param sds Standard deviation thresholds
#' @param ... Further parameters (tbd)
#' @return A character vector with channel names
#' @family Pre-ICA artifact rejection
#' @seealso [bad_epochs()]
#' @references
#' Nolan, Whelan & Reilly (2010). FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods.
#' @export

bad_chans <- function(data,
                      sds = 2,
                      ...) {

  if (is.null(data$reference)) {
    orig_ref <- NULL
    excluded <- NULL
  } else {
    orig_ref <- data$reference$ref_chans
    excluded <- data$reference$excluded
  }

  orig_chan_info <- channels(data)

  orig_names <- channel_names(data)
  # Exclude ref chan from subsequent computations (may be better to alter
  # reref_eeg...)
  data_chans <- orig_names[!(orig_names %in% data$reference$ref_chans)]
  if (!is.null(excluded)) {
    if (is.numeric(excluded)) {
      exclude <- orig_names[excluded]
    }
    message("Excluding channel(s):",
            paste(exclude, ""))
    data_chans <- data_chans[!(data_chans %in% excluded)]
  }
  data_mat <- data$signals[, data_chans]

  chan_hurst <- scale(quick_hurst(data_mat))
  chan_vars <- scale(apply(data_mat,
                           2,
                           stats::var))
  chan_corrs <- scale(colMeans(abs(stats::cor(data_mat))))
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
#' @param sds Standard deviation thresholds
#' @param ... Further parameters (tbd)
#' @return A numeric vector with epoch numbers
#' @family Pre-ICA artifact rejection
#' @seealso [bad_chans()]
#' @references
#' Nolan, Whelan & Reilly (2010). FASTER: Fully Automated Statistical Thresholding for
#' EEG artifact Rejection. J Neurosci Methods.
#' @export

bad_epochs <- function(data, sds = 2, ...) {
  chans <- channel_names(data)
  data_mat <- data.table::as.data.table(data)
  chan_means <- data_mat[, lapply(.SD, mean), .SDcols = chans]
  epoch_range <- data_mat[, lapply(.SD, function(x) max(x) - min(x)),
                          .SDcols = chans,
                          by = epoch]
  epoch_range <- epoch_range[, .(Mean = rowMeans(.SD)), by = epoch]
  epoch_range <- abs(scale(epoch_range$Mean)) > sds

  epoch_diffs <- data_mat[, lapply(.SD, mean),
                          .SDcols = chans,
                          by = epoch][, lapply(.SD, function(x) x - mean(x)),
                                      .SDcols = chans][ ,
                                                        .(Mean = rowMeans(.SD))]
  epoch_diffs <- abs(scale(epoch_diffs$Mean)) > sds

  epoch_vars <- data_mat[, lapply(.SD, var), .SDcols = chans,
                         by = epoch][, apply(.SD, 1, mean),
                                     .SDcols = chans]
  epoch_vars <- abs(scale(epoch_vars)) > sds

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
