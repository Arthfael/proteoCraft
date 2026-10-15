#' .DIANN_to_MQ_server1
#'
#' @description 
#' Server for app used to decide whether to update PTM marks.
#' 
#' @param input The server function input, cf. https://shiny.posit.co/r/getstarted/shiny-basics
#' @param output The server function's output, cf. https://shiny.posit.co/r/getstarted/shiny-basics
#' @param session Shiny session.

.DIANN_to_MQ_server1 <- function(input,
                                 output,
                                 session) {
  REMAP <- reactiveVal(tstMap)
  output$PTMs <- shiny::renderUI({
    lst <- list()
    for (i in mrks) {
      lst <- append(lst,
                    list(shiny::fluidRow(shiny::column(1L, ""),
                                         shiny::column(4L, shiny::em(i)))))
    }
    return(lst)
  })
  shiny::observeEvent(input$reMap, {
    REMAP(input$reMap)
  }, ignoreInit = TRUE)
  shiny::observeEvent(input$saveBtn, {
    tstMap <<- REMAP()
    assign("tstMap", tstMap, pos = 1)
    shiny::stopApp()
  })
  session$onSessionEnded(function() {
    assign("tstMap", tstMap, pos = 1)
    shiny::stopApp()
  })
}
