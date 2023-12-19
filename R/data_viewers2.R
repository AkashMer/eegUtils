#' Browse continuous EEG data.
#'
#' A Shiny app for scrolling through continuous EEG data interactively.
#' With EEG data (continuous), data can be viewed as a butterfly plot
#' (all electrodes overlaid) or as individual traces (electrodes "stacked").
#' The arrow buttons help to navigate the time-domain with single arrow allowing
#' 1 sec steps and double arrows allowing step length equal to `sig_length`
#'
#' @author Akash Mer
#' @importFrom dplyr %>%
#' @import ggnewscale
#' @import cowplot
#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import ggiraph
#' @param data `eeg_data` object to be plotted.
#' @param ... Other parameters passed to browsing functions.
#' @export

browse_data2 <- function(data, ...) {
  UseMethod("browse_data2", data)
}

#' @param sig_length Length of signal to be plotted initially. Default is 5 sec
#' @param n_elecs Number of electrodes to be plotted on a single screen. Can be
#' an integer indicating the number of electrodes to display or a character
#' vector specifying the electrodes. Electrodes not present in the data will
#' be ignored
#' @param downsample Reduces size of data by only plotting every 4th point,
#'  speeding up plotting considerably. Defaults to TRUE.
#'
#' @export
#' @describeIn browse_data2 Browse continuous EEG data.

browse_data2.eeg_data <- function(data,
                                  sig_length = 5,
                                  n_elecs = NULL,
                                  downsample = TRUE,
                                  ...) {

  # Verbose message
  message("Please wait! Building continuous data scroller")

  # Calculate raw sampling rate
  srate <- data$srate
  samp_diff <- min(diff(data$timings$sample))
  raw_srate <- srate * samp_diff

  # Downsample the data if user needs to
  if (downsample) {
    data <- eeg_downsample(data,
                           q = 4)
  }

  # Try to find standard locations of the electrodes
  data <- electrode_locations(data, overwrite = TRUE)

  # Build a palette for the butterfly plot based on channel types info
  colors_dict <- list(ECG = "#E41A1C",
                      Photo_Sens = "#F781BF",
                      HEOG_left = "#377EB8",
                      HEOG_right = "#984EA3",
                      HEOG = "#377EB8",
                      VEOG_upper = "#FF7F00",
                      VEOG_lower = "#4DAF4A",
                      VEOG = "#4DAF4A") %>%
    unlist()
  eog_colors <- c("#377EB8", "#984EA3", "#4DAF4A", "#FF7F00")
  if(!any(names(data$chan_info) == "type")) {
    data_channels <- data$chan_info %>%
      mutate(type = ifelse(grepl("eog", electrode, ignore.case = TRUE),
                           "EOG",
                           ifelse(grepl("ecg", electrode, ignore.case = TRUE),
                                  "ECG", "EEG"))) %>%
      select(electrode, type)
  } else {
    data_channels <- data$chan_info %>%
      select(electrode, type)
  }
  eog_channels <- data_channels %>%
    filter(type %in% c("EOG", "eog")) %>%
    mutate(colors = colors_dict[which(names(colors_dict) %in% electrode)])
  if(any(is.na(eog_channels$colors))) {
    eog_channels <- eog_channels %>%
      mutate(colors = eog_colors[nrow(eog_channels)])
  }
  colors <- vector(mode = "character", length = length(channel_names(data)))
  names(colors) <- channel_names(data)
  colors[which(names(colors) %in% eog_channels$electrode)] <-
    eog_channels$colors[which(eog_channels$electrode %in% names(colors))]
  colors[which(names(colors[which(colors == "")]) %in% names(colors_dict))] <-
    colors_dict[which(names(colors_dict) %in% names(colors[which(colors == "")]))]
  colors[which(colors == "")] <-
    scales::viridis_pal(option = "rocket")(length(colors[which(colors == "")]))

  # Build a topographic legend for butterfly plot
  topo_color_legend <- ggplot(channels(data)) +
    geom_point(aes(x, y, color = electrode), size = 2.5, na.rm = TRUE) +
    scale_color_manual(values = colors) +
    coord_equal() +
    geom_head() +
    labs(x = "",
         y = "") +
    theme_void() +
    theme(legend.position = "none")
  topo_color_legend <- suppressWarnings(plot_grid(topo_color_legend, NULL,
                                                  nrow = 1,
                                                  rel_widths = c(1,5)))

  # Calculate the initial optimized scaling based on data to display
  # Subset the electrodes
  temp <- data
  if(is.numeric(n_elecs) | is.integer(n_elecs)) {
    temp$signals <- temp$signals[,1:n_elecs]
  }
  if(is.character(n_elecs)) {
    suppressWarnings({
      temp <- temp %>%
        select_elecs(electrode = n_elecs)
    })
  }
  # Subset the data to include data according to starting point and display
  # length
  start_time <- min(data$timings$time)
  end_time <- start_time + sig_length
  if(end_time > max(data$timings$time)) {
    end_time = max(data$timings$time)
  }
  temp <- select_times(temp,
                       time_lim = c(start_time, end_time))
  # Convert into data frame for plotting
  stds <- as.data.frame(temp,
                        long = TRUE,
                        coords = FALSE,
                        events = TRUE) %>%
    summarise(std = sd(amplitude), .by = electrode) %>%
    arrange(std) %>%
    select(std) %>%
    unlist()
  if(length(stds) > 2) {
    stds <- mean(stds[2:length(stds) - 1])
  }
  else {
    stds <- mean(stds)
  }
  spacing_var <- round(stds*3)
  rm(temp)

  # App UI definition
  ui <- dashboardPage(

    # Adding a skin
    skin = "green",

    # Adding the app title
    dashboardHeader(title = "Continuous EEG data browser"),

    # Define the side bar
    dashboardSidebar(
      sliderInput(
        "time_range",
        label = "Display start time (seconds)",
        step = 1,
        min = 0,
        max = max(unique(data$timings$time)),
        value = min(unique(data$timings$time))
      ), # Closed slider input start time

      numericInput("sig_time",
                   "Display length (seconds)",
                   value = sig_length,
                   min = 1, max = 60
      ), # Closed numeric input for sig length

      checkboxInput("dc_offset",
                    "Remove DC offset",
                    value = FALSE
      ), # Closed checkbox for DC offset

      # Render the scale UI
      numericInput("spacing",
                   "Scale (μV)",
                   value = spacing_var,
                   min = 5, max = 500
      ) # Closed scale input

    ), # Closed sidebar

    # Defining the body of the app
    dashboardBody(

      # First row
      fluidRow(

        # Tab box for output of butterfly or stacked EEG plot
        tabBox(
          width = 12,
          side = "right",
          selected = "Individual",
          title = "",
          id = "outputBox",
          # Panel for butterfly plot
          tabPanel(
            "Butterfly",
            plotOutput("color_topo",
                       height = 140),
            br(),
            girafeOutput("butterfly_plot",
                         width = "100%",
                         height = "500px")
          ), # Closed butterfly tab panel
          # Panel for individual stacked plot
          tabPanel(
            "Individual",
            girafeOutput("time_plot",
                         width = "100%",
                         height = "780px")
          ) # Closed stacked tab panel
        ) # Closed output tabBox

      ), # Closed first row

      # Define second row for input
      fluidRow(align = "center",

               # Box for navigation buttons
               box(
                 background = "black",
                 width = 12,
                 radioGroupButtons("nav_buttons",
                                   choiceNames = c("<<", "<", ">", ">>"),
                                   choiceValues = c(1,2,3,4),
                                   selected = character(0),
                                   individual = TRUE,
                                   size = "lg")
               ) # Closed box for navigation buttons

      ) # Closed 2nd row

    ) # Closed dashboard body

  ) # Closed dashboard page

  # Define the server logic of the app
  server <- function(input,
                     output,
                     session) {

    output$display <- renderPrint({
      reject_data()
    })

    # Observe for navigation button clicks and change the starting point
    # and reset the buttons
    observeEvent(input$nav_buttons, {
      if(input$nav_buttons == 1) {
        updateSliderInput(
          session = session,
          "time_range",
          value = max(input$time_range - input$sig_time,
                      min(unique(data$timings$time)))
        )
        updateRadioGroupButtons(
          session = session,
          "nav_buttons",
          selected = character(0),
        )
      }
      if(input$nav_buttons == 2) {
        updateSliderInput(
          session = session,
          "time_range",
          value = max(input$time_range - 1,
                      min(unique(data$timings$time)))
        )
        updateRadioGroupButtons(
          session = session,
          "nav_buttons",
          selected = character(0),
        )
      }
      if(input$nav_buttons == 3) {
        updateSliderInput(
          session = session,
          "time_range",
          value = min(input$time_range + 1,
                      max(unique(data$timings$time)) - input$sig_time)
        )
        updateRadioGroupButtons(
          session = session,
          "nav_buttons",
          selected = character(0),
        )
      }
      if(input$nav_buttons == 4) {
        updateSliderInput(
          session = session,
          "time_range",
          value = min(input$time_range + input$sig_time,
                      max(unique(data$timings$time)) - input$sig_time)
        )
        updateRadioGroupButtons(
          session = session,
          "nav_buttons",
          selected = character(0),
        )
      }
    }) # Closed observer

    # Data processing for plotting
    tmp_data <- reactive({

      req(input$sig_time, input$spacing)

      temp <- data

      # Subset the electrodes
      if(is.numeric(n_elecs) | is.integer(n_elecs)) {
        temp$signals <- temp$signals[,1:n_elecs]
      }
      if(is.character(n_elecs)) {
        suppressWarnings({
          temp <- temp %>%
            select_elecs(electrode = n_elecs)
        })
      }

      # Subset the data to include data according to starting point and display
      # length
      start_time <- input$time_range
      end_time <- input$time_range + input$sig_time
      if(start_time < min(data$timings$time)) {
        start_time = min(data$timings$time)
      }
      if(start_time == max(data$timings$time)) {
        start_time = max(data$timings$time) - input$sig_time
      }
      if(end_time > max(data$timings$time)) {
        end_time = max(data$timings$time)
      }
      temp <- select_times(temp,
                           time_lim = c(start_time, end_time))

      # # Remove baseline if DC offset was not required
      # if (input$dc_offset) {
      #     temp <- rm_baseline(temp,
      #                         verbose = FALSE)
      # }

      # Convert into data frame for plotting
      temp <- as.data.frame(temp,
                            long = TRUE,
                            coords = FALSE,
                            events = TRUE)
    })

    if(!is.null(data$reject)) {
      disp_reject_data <- reactive({

        req(tmp_data())

        # Check the displayed data's time range
        start_time <- min(tmp_data()["time"])
        end_time <- max(tmp_data()["time"])

        # Subset the rejection data accordingly
        data$reject %>%
          filter(between(start/raw_srate, start_time, end_time) |
                   between(stop/raw_srate, start_time, end_time)|
                   between(start_time, start/raw_srate, stop/raw_srate) |
                   between(end_time, start/raw_srate, stop/raw_srate))
      })
    }


    # Render the interactive individual stacked plot
    output$time_plot <- renderGirafe({

      channels <- channel_names(data)

      # Prepare the data to plot
      if(input$dc_offset) {
        data_to_plot <- tmp_data() %>%
          mutate(electrode = fct_rev(electrode)) %>%
          mutate(amplitude = amplitude - mean(amplitude) +
                   (as.numeric(electrode)*input$spacing),
                 .by = electrode)
      } else {
        data_to_plot <- tmp_data() %>%
          mutate(electrode = fct_rev(electrode)) %>%
          mutate(amplitude = amplitude +
                   (as.numeric(electrode)*input$spacing),
                 .by = electrode)
      }

      # Build the plot
      ind_plot <- ggplot() +
        geom_line(data = data_to_plot,
                  aes(x = time, y = amplitude, color = electrode)) +
        scale_color_manual(values = rep("blue3", length(channels))) +
        new_scale_color() +
        geom_vline(data = data_to_plot,
                   aes(xintercept = event_time,
                       color = event_type), na.rm = TRUE) +
        geom_text(data = data_to_plot,
                  aes(x = event_time,
                      label = event_type,
                      color = event_type),
                  y = (length(channels) + 1.5)*input$spacing,
                  angle = 90,
                  size = 3,
                  na.rm = TRUE) +
        scale_color_distiller(palette = "Dark2") +
        scale_x_continuous(name = NULL,
                           expand = c(0,0)) +
        scale_y_continuous(name = NULL,
                           breaks = seq(input$spacing,
                                        length(channels)*input$spacing,
                                        by = input$spacing),
                           labels = rev(channels),
                           expand = c(0,0)) +
        coord_cartesian(xlim = c(min(data_to_plot$time),
                                 max(data_to_plot$time)),
                        ylim = c(0, (length(channels) + 1)*input$spacing),
                        clip = "off") +
        theme_bw() +
        theme(legend.position = "none",
              plot.margin = margin(0.8, 0.3, 0.3, 0.5, unit = "cm"))

      # Add information on rejected data
      if(!is.null(data$reject)) {
        ind_plot <- ind_plot +
          geom_rect_interactive(data = disp_reject_data(),
                                aes(xmin = start/raw_srate,
                                    xmax = stop/raw_srate,
                                    fill = status,
                                    tooltip = reason),
                                ymin = -Inf,
                                ymax = Inf,
                                alpha = 0.2)
      }

      # Return the plot
      girafe(ggobj = ind_plot,
             height_svg = 8,
             width_svg = 13,
             options = list(
               opts_hover(css = "fill:red;")
             ))

    }) # Closed render of individual plot

    output$color_topo <- renderPlot({

      topo_color_legend
    })

    output$butterfly_plot <- renderGirafe({

      # Prepare the data to plot
      if(input$dc_offset) {
        data_to_plot <- tmp_data() %>%
          mutate(amplitude = amplitude - mean(amplitude),
                 .by = electrode)
      } else {
        data_to_plot <- tmp_data()
      }

      # Build the plot
      but_plot <- ggplot() +
        geom_line(data = data_to_plot,
                  aes(x = time, y = amplitude, color = electrode)) +
        scale_color_manual(values = colors) +
        new_scale_color() +
        geom_vline(data = data_to_plot,
                   aes(xintercept = event_time,
                       color = event_type), na.rm = TRUE) +
        geom_text(data = data_to_plot,
                  aes(x = event_time,
                      label = event_type,
                      color = event_type),
                  y = max(abs(min(data_to_plot$amplitude)),
                          max(data_to_plot$amplitude))*1.1,
                  angle = 90,
                  size = 3,
                  na.rm = TRUE) +
        scale_color_distiller(palette = "Dark2")  +
        scale_x_continuous(name = NULL,
                           expand = c(0,0)) +
        scale_y_continuous(name = NULL,
                           expand = c(0,0)) +
        coord_cartesian(xlim = c(min(data_to_plot$time),
                                 max(data_to_plot$time)),
                        ylim = c(
                          -max(abs(min(data_to_plot$amplitude)),
                               max(data_to_plot$amplitude)),
                          max(abs(min(data_to_plot$amplitude)),
                              max(data_to_plot$amplitude))
                        ), clip = "off") +
        theme_bw() +
        theme(legend.position = "none",
              plot.margin = margin(0.8, 0.3, 0.3, 0.3, unit = "cm"))

      # Add information on rejected data if any
      if(!is.null(data$reject)) {
        but_plot <- but_plot +
          geom_rect_interactive(data = disp_reject_data(),
                                aes(xmin = start/raw_srate,
                                    xmax = stop/raw_srate,
                                    fill = status,
                                    tooltip = reason),
                                ymin = -Inf,
                                ymax = Inf,
                                alpha = 0.2)
      }

      # Return the plot
      girafe(ggobj = but_plot,
             height_svg = 5,
             width_svg = 13)
    })
  }
  shiny::shinyApp(ui,
                  server)
}
