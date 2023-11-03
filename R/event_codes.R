#' Tag events
#'
#' Give trigger events meaningful labels. Existing labels will be overwritten.
#' Use hierarchical labelling to tag an event with multiple labels: separate
#' labels with a "/" symbol. (e.g. "cond1" for a trigger that belongs to one
#' condition, "cond1/cond2" for a trigger that could belong to more than one
#' condition).
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class `eeg_data` or `eeg_epochs`
#' @param ... Parameters passed to S3 methods
#' @family event handlers
#' @export

tag_events <- function(data, ...) {
  UseMethod("tag_events", data)
}

#' @param trigs Character vector of trigger numbers
#' @param event_label Labels for the events.
#' @importFrom tibble as_tibble
#' @export
#' @describeIn tag_events Tag events in an `eeg_data` object

tag_events.eeg_data <- function(data,
                                trigs,
                                event_label,
                                ...) {

  if (length(trigs) != length(event_label)) {
    stop("Trigs and event_label parameters must be the same length.")
  }

  if (!any(trigs %in% unlist(list_events(data)))) {
    stop(paste0("Trigger(s) not found. Check trigger values with list_events()."))
  }

  data$events <- dplyr::left_join(data$events,
                                  data.frame(event_type = trigs,
                                             event_label = as.character(event_label),
                                             stringsAsFactors = FALSE),
                                  by = "event_type")
  data$events <- tibble::as_tibble(data$events)
  data
}

#' @describeIn tag_events Tag events in an epoched dataset
#' @export
tag_events.eeg_epochs <- function(data,
                                  trigs,
                                  event_label,
                                  ...) {

  if (length(trigs) != length(event_label)) {
    stop("Trigs and event_label parameters must be the same length.")
  }

  data$events <- dplyr::left_join(data$events,
                                  data.frame(event_type = trigs,
                                             event_label = as.character(event_label),
                                             stringsAsFactors = FALSE),
                                  by = "event_type")
  data$events <- tibble::as_tibble(data$events)
  data
}

#' List events
#'
#' List trigger types and any labels found in an `eeg_data` object.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param data An object of class `eeg_data`
#' @examples
#' list_events(demo_epochs)
#' @export
#' @family event handlers
#' @seealso [tag_events()] and [list_epochs()]

list_events <- function(data) {
  if (!is.eeg_data(data)) {
    stop("For eeg_data objects only.")
  }

  if ("event_label" %in% names(data$events)) {
    data$events[!duplicated(data$events$event_type), c("event_type", "event_label")]
  } else {
    data.frame(event_type = unique(data$events$event_type))
  }

}

#' List epochs
#'
#' List trigger types and any labels found in an `eeg_epochs` object.
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#'
#' @param data An object of class `eeg_epochs`
#' @param ... Additional arguments
#' @export
#' @family event handlers
#' @seealso [tag_events()] and [list_events()]

list_epochs <- function(data, ...) {
  UseMethod("list_epochs", data)
}

#' @describeIn list_epochs List epochs and associated events from `eeg_epochs` objects
#' @export
list_epochs.eeg_epochs <- function(data, ...) {
  data$events[, c("epoch", "event_type", "event_label")]
}

#' @describeIn list_epochs List epochs and associated events from `eeg_ICA` objects
#' @export
list_epochs.eeg_ICA <- function(data, ...) {
  data$events[, c("epoch", "event_type", "event_label")]
}

#' Modify events structure
#'
#' Get or set the values in the `events` structure of an eegUtils object.
#'
#' @examples
#' events(demo_epochs)
#' events(demo_epochs) <- mutate(events(demo_epochs),
#'  sf = dplyr::case_when(
#'          event_type %% 2 == 0 ~ "HSF",
#'          event_type %% 2 == 1 ~ "LSF",
#'  ))
#' events(demo_epochs)
#'
#' @author Matt Craddock \email{matt@@mattcraddock.com}
#' @param .data `eegUtils` object to view
#' @family event handlers
#' @export
events <- function(.data) {
  UseMethod("events", .data)
}

#' @export
events.eeg_data <- function(.data) {
  .data$events
}

#' @export
events.eeg_epochs <- function(.data) {
  .data$events
}


#' @param value Value to replace `events` structure with.
#' @rdname events
#' @export
`events<-` <- function(.data, value) {
  UseMethod("events<-", .data)
}

#' @rdname events
#' @export
`events<-.eeg_epochs` <- function(.data, value) {

  .data$events <- value
  .data
}

#' @rdname events
#' @export
`events<-.eeg_data` <- function(.data, value) {
  .data$events <- value
  .data
}

#' Shift events by a fixed amount of time
#'
#' Shift event onset and timings to account for monitor/display device delay
#' measured by a photo sensor
#'
#' @param data An `eeg_data` object
#' @param event_codes Numerical event codes for which the delay has to be
#' applied
#' @param delay Time delay in seconds
#' @param ... Additional arguments.
#' @return An eeg_data object with the specified events corrected
#' @family event handlers
#' @export
shift_events <- function(data,
                         event_codes,
                         delay,
                         ...) {
  UseMethod("shift_events", data)
}

#'@describeIn shift_events Shift events in an `eeg_data` object.
#'@export
shift_events.eeg_data <- function(data,
                                  event_codes,
                                  delay,
                                  ...) {

  # Check if event_codes supply do match the data
  if (!any(event_codes %in% unique(data$events$event_type))) {
    stop("No events found - check event codes.")
  }

  # Throw a warning in case some event_codes do not match
  if (!all(event_codes %in% unique(data$events$event_type))) {
    warning("Some events not found - check event codes.")
  }

  # Add the delay to the onset of the specified events
  new_events_table <- data$events %>%
    filter(event_type %in% event_codes) %>%
    # Multiply by sampling rate to convert to number of samples
    mutate(event_onset = event_onset + delay*data$srate) %>%
    bind_rows(
      data$events %>%
        filter(!event_type %in% event_codes)
    ) %>%
    arrange(event_onset)

  # Find the nearest samples to the updated event onset samples and replace
  data_samps <- sort(unique(data$timings$sample))
  nearest_samps <- findInterval(new_events_table$event_onset,
                                data_samps)
  new_events_table$event_onset <- data_samps[nearest_samps]

  # Update the event time to match the samples
  new_events_table$event_time <- 1 / (data$srate) * (new_events_table$event_onset - 1)

  # Replace the old event table with the new one
  data$events <- new_events_table

  return(data)
}


#' Tag epochs with labels
#'
#' Tag epochs with labels indicating details such as experimental condition,
#' based on the occurrence of event triggers from the events() structure. This
#' adds a new column to the epochs structure in an `eeg_epochs` object.
#'
#' @param .data An `eegUtils` object
#' @param ... Additional arguments.
#' @export
tag_epochs <- function(.data,
                       ...) {
  UseMethod("tag_epochs", .data)
}

#'@rdname tag_epochs
#'@export
tag_epochs.default <- function(.data,
                               ...) {
  warning("Not implemented for objects of class ", class(.data))

}

#'@describeIn tag_epochs Tag epochs in an `eeg_epochs` object.
#'@param event_type Label epochs according to specific event_types (typically a
#'  trigger)
#'@param event_label Label epochs according to specific event_labels
#'@export
tag_epochs.eeg_epochs <- function(.data,
                                  event_type = NULL,
                                  event_label = NULL,
                                  ...) {

  # need to work out how to deal with multiple labels in single epochs.
  # 1) error and demand the labels be specified?
  # 2) request specific event_types?
  # 3) allow specification of new column name in the epochs structure?

  if (!is.null(event_type) && !is.null(event_label)) {
    stop("Only event_type or event_label should be supplied, not both.")
  } else if (is.null(event_type) && is.null(event_label)) {
   epochs(.data) <- dplyr::left_join(epochs(.data),
                                     dplyr::select(events(.data),
                                                   epoch,
                                                   event_type,
                                                   event_label),
                                     by = "epoch")
  }

  if (!is.null(event_type)) {
    epochs(.data) <- dplyr::left_join(epochs(.data),
                                    dplyr::select(events(.data),
                                                  epoch,
                                                  event_type),
                                    by = "epoch")
    return(.data)
  }

  if (!is.null(event_label)) {
    epochs(.data) <- dplyr::left_join(epochs(.data),
                                      dplyr::select(events(.data),
                                                    epoch,
                                                    event_label),
                                      by = "epoch")
    return(.data)
  }
  .data
}
