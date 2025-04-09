#' DIANN_to_MQ_server2
#'
#' @description 
#' Server for app used to decide whether to update PTM marks.
#' 
#' @export

DIANN_to_MQ_server1 <- function(input, output, session) {
  REMAP <- reactiveVal(tstMap)
  output$PTMs <- shiny::renderUI({
    lst <- list()
    for (i in mrks) {
      lst <- append(lst,
                    list(shiny::fluidRow(shiny::column(1, ""),
                                         shiny::column(4, shiny::em(i)))))
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
