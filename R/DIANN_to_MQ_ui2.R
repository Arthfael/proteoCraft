#' DIANN_to_MQ_ui2
#'
#' @description 
#' UI for app used to update PTM marks.
#' 
#' @export

DIANN_to_MQ_ui2 <- function() {
  return(shiny::shinyUI(shiny::fluidPage(
    shinyjs::useShinyjs(),
    shinyWidgets::setBackgroundColor( # Doesn't work
      color = c(#"#F8F8FF",
        "#F5EDE4"),
      gradient = "linear",
      direction = "bottom"
    ),
    shinyjs::extendShinyjs(text = paste(readLines(system.file("extdata",
                                                              "jsToggleFS.txt", package = "proteoCraft")),
                                        collapse = "\n"),
                           functions = c("toggleFullScreen")),
    shiny::titlePanel(shiny::tag("u", "Map unknown PTM(s) to UniMod name(s):")),
    shiny::mainPanel(shiny::br(),
                     shiny::fluidRow(shiny::actionButton("aaRemove", "Remove all AA- and position-specific suffixes")),
                     shiny::br(),
                     shiny::fluidRow(shiny::column(3, shiny::strong(" ")),
                                     shiny::column(4, shiny::strong("Enter name here...")),
                                     shiny::column(4, shiny::strong("... or pick from UniMod matches here:"))),
                     shiny::uiOutput("PTMs"),
                     shiny::span(shiny::uiOutput("Msg"), style = "color:blue", .noWS = "outside"),
                     shinyWidgets::actionBttn("saveBtn", "Save", icon = icon("save"), color = "success", style = "pill")))))
}
