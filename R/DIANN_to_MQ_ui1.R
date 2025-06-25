#' DIANN_to_MQ_ui1
#'
#' @description 
#' Function to return UI for app used to deciding whether to update PTM marks.
#' 
#' @export

DIANN_to_MQ_ui1 <- function(dflt = tstMap) {
  if (!dflt %in% c("yes", "no")) { dlft <- "yes" }
  return(shiny::shinyUI(shiny::fluidPage(
    shinyjs::useShinyjs(),
    shinyWidgets::setBackgroundColor( # Doesn't work
      color = c(#"#F8F8FF",
        "#F5ECE6"),
      gradient = "linear",
      direction = "bottom"
    ),
    shinyjs::extendShinyjs(text = paste(readLines(system.file("extdata",
                                                              "jsToggleFS.txt", package = "proteoCraft")),
                                        collapse = "\n"),
                           functions = c("toggleFullScreen")),
    shiny::titlePanel(""),
    shiny::mainPanel(shiny::br(),
                     shiny::fluidRow(shiny::strong("Non UniMod-based PTM marks detected:")),
                     shiny::uiOutput("PTMs"),
                     shiny::strong(" -> open shiny-app to re-map them interactively to different PTM names?"),
                     shiny::fluidRow(shiny::column(8,
                                                   shiny::radioButtons("reMap",
                                                                       NULL,
                                                                       c("yes", "no"),
                                                                       dflt))),
                     shinyWidgets::actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill")))))
}
