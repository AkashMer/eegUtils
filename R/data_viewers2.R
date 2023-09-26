#'Browse continuous EEG data.
#'
#'A Shiny gadget for browsing continuous EEG data interactively.
#'With EEG data (continuous), data can be viewed as a butterfly plot
#'(all electrodes overlaid) or as individual traces (electrodes "stacked").
#'The arrow buttons help to navigate the time-domain with single arrow allowing
#'1 sec steps and double arrows allowing step length equal to `sig_length`
#'
#'@author Akash Mer
#'@import shiny
#'@import shinydashboard
#'@import shinyWidgets
#'@param data `eeg_data` object to be plotted.
#'@param ... Other parameters passed to browsing functions.
#'@export

browse_data2 <- function(data, ...) {
  UseMethod("browse_data2", data)
}

#'@param sig_length Length of signal to be plotted initially. Default is 5 sec
#'@param n_elecs Number of electrodes to be plotted on a single screen. Can be
#'an integer indicating the number of electrodes to display or a character
#'vector specifying the electrodes. Electrodes not present in the data will
#'be ignored
#'@param downsample Reduces size of data by only plotting every 4th point,
#'  speeding up plotting considerably. Defaults to TRUE.
#'@param scale Apply amplitude minmax scaling. Defaults to TRUE.
#'  If TRUE, Scaling factor is 50 μV by default, which can be changed in the app.
#'  Range accepted from 5 to 500 μV
#'
#'@export
#'@describeIn browse_data2 Browse continuous EEG data.

browse_data2.eeg_data <- function(data,
                                 sig_length = 5,
                                 n_elecs = NULL,
                                 downsample = TRUE,
                                 scale = TRUE,
                                 ...) {

  library(dplyr)
  library(plotly)
  srate <- data$srate

  # Downsample the data if user needs to
  if (downsample) {
    data <- eeg_downsample(data,
                           q = 4)
  }

  # App UI definition
  ui <- dashboardPage(

    # Adding a skin
    skin = "green",

    # Adding the app title
    dashboardHeader(title = "Continuous EEG data browser"),

    # Disable the side bar
    dashboardSidebar(
      sliderInput(
        "time_range",
        label = "Display start time",
        step = 1,
        min = 0,
        max = max(unique(data$timings$time)),
        value = min(unique(data$timings$time))
      ), # Closed slider input start time

      numericInput("sig_time",
                   "Display length",
                   value = sig_length,
                   min = 1, max = 60
      ), # Closed numeric input for sig length

      checkboxInput("dc_offset",
                    "Remove DC offset",
                    value = TRUE
      ), # Closed checkbox for DC offset

      checkboxInput("scale_check",
                    "Apply scaling to amplitude",
                    value = scale
      ), # Closed check box for scale manipulation

      # UI for scale factor control
      numericInput("scale_value",
                   "Scale",
                   value = 50,
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
              plotlyOutput("butterfly_plot")
            ), # Closed butterfly tab panel
            # Panel for individual stacked plot
            tabPanel(
              "Individual",
              plotlyOutput("time_plot", height = 780)
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

      # Remove baseline if DC offset was not required
      if (input$dc_offset) {
        temp <- rm_baseline(temp,
                            verbose = FALSE)
      }

      # Convert into data frame for plotting
      temp <- as.data.frame(temp,
                            long = TRUE,
                            coords = FALSE,
                            events = TRUE)

      # Scale EEG signal by amplitude range customization
      if(input$scale_check) {
        scale_min <- -input$scale_value
        scale_max <- input$scale_value
        min <- min(temp$amplitude)
        max <- max(temp$amplitude)
        scaling_factor <- (scale_max - scale_min)/(max - min)
        temp$amplitude <- (temp$amplitude - min)*scaling_factor + scale_min
      }

      return(temp)
    })

    # Render the interactive individual stacked plot
    output$time_plot <- renderPlotly({

      # Get event timings for plotting
      event_times <- unique(tmp_data()$event_time) %>%
        na.omit() %>%
        as.vector()
      line_list <- list()
      if(length(event_times) > 0) {
        for(i in 1:length(event_times)) {
          line_list[[i]] <- list(type = "line",
                                 fillcolor = "darkgreen",
                                 line = list(color = "darkgreen"),
                                 x0 = event_times[i],
                                 x1 = event_times[i],
                                 xref = "x",
                                 y0 = 0,
                                 y1 = 1,
                                 yref = "paper")
        }
      }

      # Define plot margins
      m <- list(
        l = 90,
        r = 20,
        b = 60,
        t = 20
      )

      # Define the interactive plot
      init_plot <- tmp_data() %>%
        group_by(electrode) %>%
        do(p=plot_ly(., x = ~time, y = ~amplitude,
                     type = "scatter", mode = "lines",
                     split = ~electrode, color = I("black")) %>%
             add_annotations(
               text = ~unique(electrode),
               x = 0,
               y = 0.8,
               yref = "paper",
               xref = "paper",
               xanchor = "right",
               yanchor = "top",
               showarrow = FALSE,
               font = list(size = 12)
             ) %>%
             layout(shapes = line_list,
                    xaxis = list(title = "Time(s)"),
                    yaxis = list(showticklabels = FALSE),
                    margin = m)) %>%
        subplot(nrows = length(unique(tmp_data()$electrode)),
                shareX = TRUE, margin = 0) %>%
        hide_legend()

      # Return the plot output
      init_plot
    }) # Closed render of individual plot

    output$butterfly_plot <- renderPlotly({

      # Get event timings for plotting
      event_times <- unique(tmp_data()$event_time) %>%
        na.omit() %>%
        as.vector()
      line_list <- list()
      if(length(event_times) > 0) {
        for(i in 1:length(event_times)) {
          line_list[[i]] <- list(type = "line",
                                 fillcolor = "black",
                                 line = list(color = "black"),
                                 x0 = event_times[i],
                                 x1 = event_times[i],
                                 xref = "x",
                                 y0 = 0,
                                 y1 = 1,
                                 yref = "paper")
        }
      }


      # Define the interactive plot
      plot <- tmp_data() %>%
        plot_ly(., x =~time, y = ~amplitude,
                     type = "scatter", mode = "lines",
                     split = ~electrode) %>%
        layout(shapes = line_list,
               xaxis = list(title = "Time(s)"),
               yaxis = list(title = "Amplitude(μV)")
        ) %>%
        hide_legend()

      # Return the plot output
      plot
    })
  }
  shiny::shinyApp(ui,
                  server)
}
