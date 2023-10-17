#' Filter EEG data 2
#'
#' Perform IIR or FIR filtering on input EEG data of class `eeg_data` or
#' `eeg_epochs`. WARNING: with epoched data, epoch boundaries are currently
#' ignored, which can result in minor edge artifacts.
#'
#' low_freq and high_freq are the low and high cutoff frequencies. Pass low freq
#' or high freq alone to perform high-pass or low-pass filtering respectively.
#' For band-pass or band-stop filters, pass both low_freq and high_freq.
#'
#' If low_freq < high_freq, bandpass filtering is performed.
#'
#' If low_freq > high_freq, bandstop filtering is performed.
#'
#' Note that it is recommended to first zero-mean the signal using either
#' channel means or by-channel epoch means.
#'
#' The function allows parallelization using the `future` package, e.g. using
#' `plan(multiprocess)`
#'
#' @section FIR versus IIR filtering: Finite Impulse Response (FIR) filtering is
#'   performed using an overlap-add FFT method. Note that this only performs a
#'   single-pass; the data is shifted back in time by the group delay of the
#'   filter to compensate for the phase delay imposed by the linear filtering
#'   process. Infinite Impulse Response (IIR) filtering is performed using a
#'   two-pass (once forwards, once reversed) method to correct for phase
#'   alignment. Note that the Butterworth filter designs used here can become
#'   numerically unstable with only a small increase in filter order. For most
#'   purposes, use FIR filters.
#'
#' @examples
#' plot_psd(eeg_filter(demo_epochs, low_freq = 1, high_freq = 30))
#' plot_psd(eeg_filter(demo_epochs, low_freq = 12, high_freq = 8))
#' plot_psd(eeg_filter(demo_epochs, low_freq = 12, high_freq = 8, method = "iir"))
#'
#' @param data An `eeg_data` or `eeg_epochs` object to be filtered.
#' @param ... Additional parameters.
#' @return An object of the original class with signals filtered according to
#'   the user's specifications
#' @export

eeg_filter2 <- function(data, ...) {
  UseMethod("eeg_filter2", data)
}

#' @param low_freq Low cutoff frequency.
#' @param high_freq High cutoff frequency.
#' @param filter_length Defaults to "auto", which automatically estimates filter order for the specified filter characteristics (defaults to 4 if method = "iir").
#' @param low_trans_bw Transition bandwidth of high pass (Hz). "auto" or an integer.
#'   "auto" attempts to determine a suitable transition bandwidth using the
#'   heuristic given below. Ignored if `method` = "iir".
#' @param high_trans_bw Transition bandwidth of low pass (Hz). "auto" or an integer.
#'   "auto" attempts to determine a suitable transition bandwidth using the
#'   heuristic given below. Ignored if `method` = "iir".
#' @param method "fir" (Finite Impulse Response) or "iir" (Infinite Impulse
#'   Response). Defaults to "fir".
#' @param window Windowing function to use (FIR filtering only). Defaults to
#'   "hamming"; currently only "hann", "hamming" and "blackman" available.
#' @param demean Remove DC component (i.e. channel/epoch mean) before filtering. Defaults to TRUE.
#' @rdname eeg_filter2
#' @export
eeg_filter2.eeg_data <- function(data,
                                 low_freq = NULL,
                                 high_freq = NULL,
                                 filter_length = "auto",
                                 low_trans_bw = "auto",
                                 high_trans_bw = "auto",
                                 method = "fir",
                                 window = "hamming",
                                 demean = TRUE,
                                 ...) {

  if (identical(method, "iir") | identical(method, "fir")) {

    # Get the filter parameters
    filt_pars <- parse_filt_params(data$srate,
                                   method,
                                   low_freq,
                                   high_freq,
                                   low_trans_bw,
                                   high_trans_bw,
                                   filter_length,
                                   window)

  } else {
    stop("Unknown method; valid methods are 'fir' and 'iir'.")
  }

  # Prepare filter parameters to prepare for coefficient calculation
  if(identical(method, "iir")) {

    # Get the normalized critical frequencies
    W <- c(low_freq, high_freq)
    W <- W/(data$srate/2)

    # Calculate the IIR filter coefficients using a butterworth filter
    filt_coef <- signal::butter(
      n = filt_pars$filter_length,
      W = W,
      type = filt_pars$filt_type
    )

  } else { # FIR filters

    if(identical(filt_pars$filt_type, "high")) {

      freq <- c(filt_pars$stop_lower, filt_pars$pass_lower, data$srate/2)
      gain <- c(0,1,1)
      if(freq[1] != 0) {
        freq <- c(0, freq)
        gain <- c(0, gain)
      }

    } else if(identical(filt_pars$filt_type, "low")) {

      freq <- c(0, filt_pars$pass_upper, filt_pars$stop_upper)
      gain <- c(1,1,0)
      if(freq[length(freq)] != data$srate/2) {
        freq <- c(freq, data$srate/2)
        gain <- c(gain, 0)
      }

    } else if(identical(filt_pars$filt_type, "pass")) {

      freq <- c(filt_pars$stop_lower, filt_pars$pass_lower,
                filt_pars$pass_upper, filt_pars$stop_upper)
      gain <- c(0,1,1,0)
      if(freq[1] != 0) {
        freq <- c(0, freq)
        gain <- c(0, gain)
      }
      if(freq[length(freq)] != data$srate/2) {
        freq <- c(freq, data$srate/2)
        gain <- c(gain, 0)
      }

    } else { #band-stop

      freq <- c(filt_pars$pass_lower, filt_pars$stop_lower,
                filt_pars$stop_upper, filt_pars$pass_upper)
      gain <- c(1,0,0,1)
      freq <- freq[order(freq)]
      gain <- gain[order(freq)]
      if(freq[0] != 0) {
        freq <- c(0, freq)
        gain <- c(1, gain)
      }
      if(freq[length(freq)] != data$srate/2) {
        freq <- c(freq, data$srate/2)
        gain <- c(gain, 1)
      }
      # Make sure if stop bands are sufficiently separated
      if(any(abs(diff(gain, differences = 2)) > 1)) {
        stop("Stop bands are not sufficiently separated")
      }

    }

    # Normalize the frequencies
    freq <- freq/(data$srate/2)

    # Run checks
    if(filt_pars$filter_length > nrow(data$signals)) {
      message(paste("filter_length", filter_length, "is longer than the",
                    "signal", nrow(data$signals), ", distortion is likely.",
                    "Reduce filter length or filter a longer signal"))
    }

    # Apply the requested smoothing window
    smoothing_window <- select_window(window,
                                      filt_pars$filter_length)

    # Calculate FIR filter coefficients required for FFT
    filt_coef <- fir2(n = filt_pars$filter_length,
                      f = freq,
                      m = gain,
                      window = smoothing_window)
  }

  # Remove DC component
  if (demean) {
    data <- rm_baseline(data)
  }

  # Filter the data using the requested method and other parameters
  if (identical(method, "iir")) {
    data <- run_iir_n(data,
                      filt_coef)
  } else {
    data <- run_fir(data,
                    filt_coef,
                    filt_pars$filter_length)
  }
  data$signals <- tibble::as_tibble(data$signals)
  data
}


#' @rdname eeg_filter2
#' @export
eeg_filter2.eeg_epochs <- function(data,
                                   low_freq = NULL,
                                   high_freq = NULL,
                                   filter_length = "auto",
                                   low_trans_bw = "auto",
                                   high_trans_bw = "auto",
                                   method = "fir",
                                   window = "hamming",
                                   demean = TRUE,
                                   ...) {

  if (identical(method, "iir") | identical(method, "fir")) {

    # Get the filter parameters
    filt_pars <- parse_filt_params(data$srate,
                                   method,
                                   low_freq,
                                   high_freq,
                                   low_trans_bw,
                                   high_trans_bw,
                                   filter_length,
                                   window)

  } else {
    stop("Unknown method; valid methods are 'fir' and 'iir'.")
  }

  # Prepare filter parameters to prepare for coefficient calculation
  if(identical(method, "iir")) {

    # Get the normalized critical frequencies
    W <- c(low_freq, high_freq)
    W <- W/(data$srate/2)

    # Calculate the IIR filter coefficients using a butterworth filter
    filt_coef <- signal::butter(
      n = filt_pars$filter_length,
      W = W,
      type = filt_pars$filt_type
    )

  } else { # FIR filters

    if(identical(filt_pars$filt_type, "high")) {

      freq <- c(stop_lower, pass_lower, data$srate/2)
      gain <- c(0,1,1)
      if(freq[1] != 0) {
        freq <- c(0, freq)
        gain <- c(0, gain)
      }

    } else if(identical(filt_pars$filt_type, "low")) {

      freq <- c(0, pass_upper, stop_upper)
      gain <- c(1,1,0)
      if(freq[length(freq)] != data$srate/2) {
        freq <- c(freq, data$srate/2)
        gain <- c(gain, 0)
      }

    } else if(identical(filt_pars$filt_type, "pass")) {

      freq <- c(stop_lower, pass_lower, pass_upper, stop_upper)
      gain <- c(0,1,1,0)
      if(freq[1] != 0) {
        freq <- c(0, freq)
        gain <- c(0, gain)
      }
      if(freq[length(freq)] != data$srate/2) {
        freq <- c(freq, data$srate/2)
        gain <- c(gain, 0)
      }

    } else { #band-stop

      freq <- c(pass_lower, stop_lower, stop_upper, pass_upper)
      gain <- c(1,0,0,1)
      freq <- freq[order(freq)]
      gain <- gain[order(freq)]
      if(freq[0] != 0) {
        freq <- c(0, freq)
        gain <- c(1, gain)
      }
      if(freq[length(freq)] != data$srate/2) {
        freq <- c(freq, data$srate/2)
        gain <- c(gain, 1)
      }
      # Make sure if stop bands are sufficiently separated
      if(any(abs(diff(gain, differences = 2)) > 1)) {
        stop("Stop bands are not sufficiently separated")
      }

    }

    # Normalize the frequencies
    freq <- freq/(data$srate/2)

    # Run checks
    if(filt_pars$filter_length > nrow(data$signals)) {
      message(paste("filter_length", filter_length, "is longer than the",
                    "signal", nrow(data$signals), ", distortion is likely.",
                    "Reduce filter length or filter a longer signal"))
    }

    # Apply the requested smoothing window
    smoothing_window <- select_window(window,
                                      filt_pars$filter_length)

    # Calculate FIR filter coefficients required for FFT
    filt_coef <- signal::fir2(n = filt_pars$filter_length,
                              f = freq,
                              m = gain,
                              window = smoothing_window)
  }

  # Remove DC component
  if (demean) {
    data <- rm_baseline(data)
  }

  # Filter the data using the requested method and other parameters
  if (identical(method, "iir")) {
    data <- run_iir_n(data,
                      filt_coef)
  } else {
    data <- run_fir(data,
                    filt_coef,
                    filt_pars$filter_length)
  }
  data$signals <- tibble::as_tibble(data$signals)
  data
}

#' @rdname eeg_filter2
#' @export
eeg_filter2.eeg_group <- function(data,
                                  low_freq = NULL,
                                  high_freq = NULL,
                                  filter_length = "auto",
                                  low_trans_bw = "auto",
                                  high_trans_bw = "auto",
                                  method = "fir",
                                  window = "hamming",
                                  demean = TRUE,
                                  ...) {

  stop("Filtering not supported on eeg_group objects.")
}


#' Run FIR filter using overlap-add FFT
#'
#' @param data Data to be filtered.
#' @param filt_coef Filter coefficients
#' @param filter_order Order of filter
#' @importFrom future.apply future_lapply
#' @importFrom purrr map_df
#' @importFrom tibble as_tibble
#' @keywords internal
run_fir <- function(data,
                    filt_coef,
                    filter_order) {

  fft_length <- length(filt_coef) * 2 - 1
  fft_length <- stats::nextn(fft_length) #length(filt_coef) * 2 - 1
  sig_length <- nrow(data$signals)
  # pad the signals with zeros to help with edge effects
  pad_zeros <- stats::nextn(sig_length + fft_length - 1) - sig_length
  pad_zeros <- 2 * round(pad_zeros / 2)
  data$signals <- purrr::map_df(data$signals,
                                ~pad(.,
                                     pad_zeros))
  data$signals <- future.apply::future_lapply(data$signals,
                                              signal::fftfilt,
                                              b = filt_coef,
                                              n = fft_length)


  # fftfilt filters once and thus shifts everything in time by the group delay
  # of the filter (half the filter order). Here we correct for both the
  # padding and the group delay
  data$signals <- purrr::map_df(data$signals,
                                ~fix_grpdelay(.,
                                              pad_zeros,
                                              filter_order / 2))
  data$signals <- tibble::as_tibble(data$signals)
  data
}

run_iir_n <- function(data,
                      filt_coef) {

  data$signals <- future.apply::future_lapply(data$signals,
                                              signal::filtfilt,
                                              filt = filt_coef,
                                              a = 1)
  data$signals <- tibble::as_tibble(data$signals)
  data
}


#' Parse filter frequency input
#'
#' Parses the frequencies input by the user, converting them to a fraction of
#' the sampling rate and setting the filter type (low-pass, high-pass,
#' band-pass, band-stop) appropriately.
#'
#' @param srate Sampling rate (Hz)
#' @param method "iir" or "fir" method.
#' @param l_freq low frequency cutoff (Hz)
#' @param h_freq High frequency cutoff (Hz)
#' @param l_trans_bw Transition bandwidth of high pass (Hz)
#' @param h_trans_bw Transition bandwidth of low pass (Hz)
#' @param filter_length Filter length of the filter representing the filter order
#' @param window Windowing function to use (FIR filtering only). Defaults to
#'   "hamming"; currently only "hann", "hamming" and "blackman" available.
#' @keywords internal
parse_filt_params <- function(srate,
                              method,
                              l_freq,
                              h_freq,
                              l_trans_bw,
                              h_trans_bw,
                              filter_length,
                              window) {

  if (length(l_freq) > 1 | length(h_freq) > 1) {
    stop("Only one number should be passed to low_freq or high_freq.")
  }

  if (is.null(l_freq) & is.null(h_freq)) {
    stop("At least one frequency must be specified.")
  }

  reverse <- FALSE
  lbw <- NA
  hbw <- NA
  max_tbw <- NA
  length_factors <- data.frame(
    window_method = c("hann", "hamming", "blackman"),
    factor = c(3.1, 3.3, 5.0)
  )
  l_stop = h_stop = NULL

  # Determine the type of filter
  if (is.null(h_freq)) {
    filt_type <- "high"
    message(paste("High-pass",
                  toupper(method),
                  "filter at",
                  l_freq,
                  "Hz"))

    if(method == "iir") {
      l_stop = l_freq
    } else { # FIR method

      message(paste("Lower passband edge: ", l_freq))
      # Determine transition bandwidth of high pass
      if (identical(l_trans_bw, "auto")) {
        lbw <- min(max(l_freq * 0.25, 2), l_freq)
      } else if (!is.character(l_trans_bw) & length(l_trans_bw) == 1) {
        lbw <- l_trans_bw
      } else {
        stop("low_trans_bandwidth has to be set to auto or a single number")
      }
      message(paste("Lower transition bandwidth: ", lbw, "Hz"))

      # Determine the stop band frequency
      l_stop <- low_freq - lbw

    }

  } else if (is.null(l_freq)) {
    filt_type <- "low"
    message(paste("Low-pass",
                  toupper(method),
                  "filter at",
                  h_freq, "Hz"))

    if(method == "iir") {
      h_stop = h_freq
    } else { # FIR method

      message(paste("Higher passband edge: ", h_freq))
      # Determine transition bandwidth of low pass
      if (identical(h_trans_bw, "auto")) {
        hbw <- min(max(h_freq * 0.25, 2), srate / 2 - h_freq)
      } else if (!is.character(h_trans_bw) & length(h_trans_bw) == 1) {
        hbw <- h_trans_bw
      } else {
        stop("high_trans_bandwidth has to be set to auto or a single number")
      }
      message(paste("Upper transition bandwidth: ", hbw, "Hz"))

      # Determine the stop band frequency
      h_stop <- high_freq + hbw

    }

  } else {
    if (l_freq > h_freq) {

      filt_type <- "stop"
      message(paste("Band-stop",
                    toupper(method),
                    "filter from",
                    h_freq,
                    "-",
                    l_freq, "Hz."))
      reverse <- TRUE

    } else if (l_freq < h_freq) {

      filt_type <- "pass"
      message(paste("Band-pass",
                    toupper(method),
                    "filter from",
                    l_freq, "-",
                    h_freq, "Hz"))

      reverse <- FALSE
    }

    if(method == "iir") {
      l_stop = l_freq
      h_stop = h_freq
    } else { # FIR method

      message(paste("Lower passband edge: ", l_freq))
      # Determine transition bandwidth of high pass
      if (identical(l_trans_bw, "auto")) {
        lbw <- min(max(l_freq * 0.25, 2), l_freq)
      } else if (!is.character(l_trans_bw) & length(l_trans_bw) == 1) {
        lbw <- l_trans_bw
      } else {
        stop("low_trans_bandwidth has to be set to auto or a single number")
      }
      message(paste("Lower transition bandwidth: ", lbw, "Hz"))

      message(paste("Higher passband edge: ", h_freq))
      # Determine transition bandwidth of low pass
      if (identical(h_trans_bw, "auto")) {
        hbw <- min(max(h_freq * 0.25, 2), srate / 2 - h_freq)
      } else if (!is.character(h_trans_bw) & length(h_trans_bw) == 1) {
        hbw <- h_trans_bw
      } else {
        stop("high_trans_bandwidth has to be set to auto or a single number")
      }
      message(paste("Upper transition bandwidth: ", hbw, "Hz"))

      # Determine the stop band frequency
      l_stop <- l_freq - lbw
      h_stop <- h_freq + hbw
      if(reverse) {
        l_stop = l_stop + lbw
        l_freq = l_freq + lbw
        h_stop = h_stop - hbw
        h_freq = h_freq - hbw
        max_tbw <- abs(l_freq - h_freq) / 2
      }

    }

  }

  # Determine the filter_length
  if(identical(filter_length, "auto")) {
    if(identical(method, "iir")) {
      filter_length <- 2
      message(paste("Effective filter order:",
                    filter_length * 2,
                    "(two-pass)"))
    } else { # FIR method
      if (!is.na(match(window, length_factors$window_method))) {
        filt_ord <- length_factors$factor[match(
          window,
          length_factors$window_method
        )]
      } else{
        stop("Only hamming, hann and blackman windows currently implemented.")
      }
      filt_ord <- filt_ord / min(c(lbw, hbw, max_tbw), na.rm = TRUE)
      # Convert to samples
      filter_length <- filt_ord * srate
      # Make sure it is even and a whole number
      filter_length <- max(round(filter_length, 1))
      filter_length <- ceiling(filter_length / 2) * 2
      message(paste("Filter Length:", filter_length, "samples", "(", filt_ord, "s)"))
    }
  }

  # Return the required value
  if(!reverse) {
    list("srate" = srate,
         "filt_type" = filt_type,
         "pass_lower" = l_freq,
         "pass_upper" = h_freq,
         "stop_lower" = l_stop,
         "stop_upper" = h_stop,
         "filter_length" = filter_length,
         "window" = window)
  } else {
    list("srate" = srate,
         "filt_type" = filt_type,
         "stop_lower" = l_freq,
         "stop_upper" = h_freq,
         "pass_lower" = l_stop,
         "pass_upper" = h_stop,
         "filter_length" = filter_length,
         "window" = window)
  }
}

#' Create windowing function
#'
#' Create a windowing function for use in creating a windowed-sinc kernel
#'
#' @param type Window function to apply
#' @param m Filter order
#' @param a alpha/beta to be used for some window functions
#' @keywords internal

select_window <- function(type,
                          m,
                          a = NULL) {

  m <- m + 1
  w <- switch(type,
              "bartlett" = signal::bartlett(m),
              "hann" = signal::hanning(m),
              "hamming" = signal::hamming(m),
              "blackman" = signal::blackman(m),
              "kaiser" = signal::kaiser(m, 5.653))
  w
}
